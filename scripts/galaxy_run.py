#!/usr/bin/env python3
"""Run the DNA and RNA analysis chains on Galaxy from the command line.

Everything the checklist asks you to click is encoded here: tool ids, versions and every parameter. Each step's
inputs are wired to the previous step's outputs, outputs are renamed to the canonical names
`scripts/import_galaxy.py` expects, and the script resumes — a step whose output already exists in the history
is skipped, so an interrupted run can simply be restarted.

Before anything is submitted, every parameter path is checked against the tool definition that Galaxy currently
serves, so a renamed parameter fails loudly here instead of silently producing a wrong result.

  export GALAXY_API_KEY=...
  python3 scripts/galaxy_run.py --history "TUMOR01 tumour profiling" --arm dna --dry-run
  python3 scripts/galaxy_run.py --history "TUMOR01 tumour profiling" --arm dna
  python3 scripts/galaxy_run.py --history "TUMOR01 tumour profiling" --arm rna
"""
import argparse, os, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
try:
    from bioblend.galaxy import GalaxyInstance
except ImportError:
    import subprocess
    for cand in [os.environ.get("CONDA_PREFIX", "") + "/bin/python",
                 "/opt/miniconda3/envs/tp-py/bin/python"]:
        if cand and Path(cand).exists() and cand != sys.executable:
            if subprocess.run([cand, "-c", "import bioblend"], capture_output=True).returncode == 0:
                os.execv(cand, [cand, __file__] + sys.argv[1:])
    sys.exit("needs BioBlend: /opt/miniconda3/envs/tp-py/bin/python -m pip install bioblend")

T = {  # tool ids verified against usegalaxy.eu on 2026-08-30
    "umi_extract": "toolshed.g2.bx.psu.edu/repos/iuc/umi_tools_extract/umi_tools_extract/1.1.6+galaxy0",
    "bwa_mem": "toolshed.g2.bx.psu.edu/repos/devteam/bwa/bwa_mem/0.7.19+galaxy1",
    "umi_dedup": "toolshed.g2.bx.psu.edu/repos/iuc/umi_tools_dedup/umi_tools_dedup/1.1.6+galaxy0",
    "ampliconclip": "toolshed.g2.bx.psu.edu/repos/iuc/samtools_ampliconclip/samtools_ampliconclip/1.22+galaxy2",
    "mosdepth": "toolshed.g2.bx.psu.edu/repos/iuc/mosdepth/mosdepth/0.3.8+galaxy0",
    "mutect2": "toolshed.g2.bx.psu.edu/repos/iuc/gatk4_mutect2/gatk4_mutect2/4.6.2.0+galaxy0",
    "vardict": "toolshed.g2.bx.psu.edu/repos/iuc/vardict_java/vardict_java/1.8.4+galaxy0",
    "lofreq": "toolshed.g2.bx.psu.edu/repos/iuc/lofreq_call/lofreq_call/2.1.5+galaxy3",
    "rna_star": "toolshed.g2.bx.psu.edu/repos/iuc/rgrnastar/rna_star/2.7.8a+galaxy1",
    "arriba": "toolshed.g2.bx.psu.edu/repos/iuc/arriba/arriba/2.5.1+galaxy1",
    "arriba_filters": "toolshed.g2.bx.psu.edu/repos/iuc/arriba_get_filters/arriba_get_filters/2.5.1+galaxy1",
}
DNA_UMI = r"^(?P<umi_1>.{12})(?P<discard_1>AGTCGTCTCGAAGT?){s<=2}"
RNA_UMI = r"^(?P<umi_1>.{12})(?P<discard_1>CTGGATAGTACGCT){s<=2}"
GTF_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/"
           "gencode.v44.primary_assembly.annotation.gtf.gz")


def flatten(inputs, prefix=""):
    """Every legal parameter path of a tool, so a payload can be validated before it is submitted."""
    paths = set()
    for i in inputs or []:
        nm = i.get("name", "")
        full = f"{prefix}{nm}"
        t = i.get("type")
        if t == "conditional":
            paths.add(full)
            for c in i.get("cases", []) or []:
                paths |= flatten(c.get("inputs"), f"{full}|")
        elif t in ("section", "repeat"):
            paths |= flatten(i.get("inputs"), f"{full}|")
        else:
            paths.add(full)
    return paths


class Runner:
    def __init__(self, gi, history_id, dry_run=False):
        self.gi, self.hid, self.dry = gi, history_id, dry_run
        self._tools = {}
        self.refresh()

    def refresh(self):
        self.datasets = {}
        for d in self.gi.histories.show_history(self.hid, contents=True, deleted=False):
            if d.get("visible") is False:
                continue
            self.datasets.setdefault(d["name"], []).append(d)

    def find(self, name, states=("ok",)):
        for d in self.datasets.get(name, []):
            if d.get("state") in states:
                return d
        return None

    def need(self, name):
        d = self.find(name)
        if not d:
            raise SystemExit(f"missing input dataset {name!r} in the history — upload it first "
                             f"(scripts/galaxy_upload.py)")
        return {"src": "hda", "id": d["id"]}

    def tool_paths(self, key):
        if key not in self._tools:
            self._tools[key] = flatten(self.gi.tools.show_tool(T[key], io_details=True).get("inputs"))
        return self._tools[key]

    def validate(self, key, inputs):
        legal = self.tool_paths(key)
        bad = [k for k in inputs if k not in legal]
        if bad:
            raise SystemExit(f"{key}: Galaxy does not know these parameters: {bad}\n"
                             f"  (the tool was probably updated; check with "
                             f"`show_tool('{T[key]}', io_details=True)`)")

    def run(self, key, label, inputs, rename):
        """Run one tool unless its renamed outputs are already present."""
        done = all(self.find(n) for n in rename.values())
        if done:
            print(f"skip  {label} — outputs already in the history")
            return
        self.validate(key, inputs)
        if self.dry:
            print(f"would run {label} ({T[key].split('/')[-2]}) -> {list(rename.values())}")
            return
        print(f"run   {label} ...", flush=True)
        out = self.gi.tools.run_tool(self.hid, T[key], inputs)
        ids = {o["output_name"]: o["id"] for o in out["outputs"]} if "output_name" in (out["outputs"][0] or {}) \
            else {}
        if not ids:                       # older API shape: outputs come back in tool order
            ids = {}
            info = self.gi.tools.show_tool(T[key], io_details=True)
            for spec, o in zip(info.get("outputs", []), out["outputs"]):
                ids[spec["name"]] = o["id"]
        for out_name, new_name in rename.items():
            if out_name in ids:
                self.gi.histories.update_dataset(self.hid, ids[out_name], name=new_name)
        self.wait([i for i in ids.values()], label)
        self.refresh()

    def wait(self, dataset_ids, label, poll=30):
        while True:
            states = []
            for i in dataset_ids:
                try:
                    states.append(self.gi.datasets.show_dataset(i).get("state"))
                except Exception:
                    states.append("unknown")
            if any(s == "error" for s in states):
                raise SystemExit(f"{label}: a job failed — open the history and read the job's stderr")
            if all(s in ("ok", "empty") for s in states):
                print(f"      {label} finished", flush=True)
                return
            print(f"      {label}: {', '.join(sorted(set(states)))}", flush=True)
            time.sleep(poll)


def dna_arm(r, alias):
    r.run("umi_extract", "1 UMI extract (DNA)", {
        "input_type_cond": "paired",
        "input_type_cond|input_read1": r.need(f"{alias}_DNA_R1.fastq.gz"),
        "input_type_cond|input_read2": r.need(f"{alias}_DNA_R2.fastq.gz"),
        "input_type_cond|bc_pattern": DNA_UMI,
        "extract_method_cond": "regex",
        "log": "true",
    }, {"out": "dna_umi_R1.fastq.gz", "out2": "dna_umi_R2.fastq.gz", "out_log": "dna_umi_extract.log"})

    r.run("bwa_mem", "2 BWA-MEM (DNA)", {
        "reference_source": "cached",
        "reference_source|ref_file": "hg38",
        "fastq_input": "paired",
        "fastq_input|fastq_input1": r.need("dna_umi_R1.fastq.gz"),
        "fastq_input|fastq_input2": r.need("dna_umi_R2.fastq.gz"),
        "rg": "set",
        "rg|read_group_id_conditional": "true",
        "rg|read_group_id_conditional|ID": alias,
        "rg|read_group_sm_conditional": "true",
        "rg|read_group_sm_conditional|SM": alias,
        "rg|read_group_lb_conditional": "true",
        "rg|read_group_lb_conditional|LB": "VariantPlex",
        "rg|PL": "ILLUMINA",
        "analysis_type": "illumina",
    }, {"bam_output": "dna_aligned.bam"})

    r.run("umi_dedup", "3 UMI dedup (DNA)", {
        "input": r.need("dna_aligned.bam"),
        "bc": "read_id",                       # the UMI sits in the read name, put there by UMI-tools extract
        "umi|method": "directional",
        "sambam|paired": "true",
        "log": "true",
    }, {"output": "dna_dedup.bam", "out_log": "dna_dedup.log"})

    r.run("ampliconclip", "4 primer clip (DNA)", {
        "input_bam": r.need("dna_dedup.bam"),
        "input_bed": r.need("panel_dna.bed"),
        "hard_clip_mode": "true",
        "both_ends": "true",
    }, {"output_bam": "dna_clipped.bam"})

    r.run("mosdepth", "5 coverage", {
        "input_alignment": r.need("dna_clipped.bam"),
        "per_base_coverage": "false",
        "window": "bed",
        "window|region_file": r.need("panel_dna.bed"),
    }, {"output_regions_bed": "coverage.regions.bed", "output_summary": "coverage.summary.txt"})

    r.run("mutect2", "6a Mutect2 (tumour-only)", {
        "mode": "tumor_only",
        "mode|tumor": r.need("dna_clipped.bam"),
        "reference_source": "cached",
        "reference_source|reference_sequence": "hg38",
        "optional": "yes",
        "optional|germline_resource": r.need("af-only-gnomad.hg38.vcf.gz"),
        "optional|ival_type": "ival_file",
        "optional|ival_type|intervals": r.need("panel_dna.bed"),
        "gzipped_output": "false",
    }, {"output_vcf": "mutect2.vcf"})

    r.run("vardict", "6b VarDict", {
        "select_mode": "single",
        "select_mode|tumor": r.need("dna_clipped.bam"),
        "select_mode|interval_file": r.need("panel_dna.bed"),
        "reference_source": "cached",
        "reference_source|ref_file": "hg38",
        "advancedsettings|f": "0.01",
    }, {"all_variants": "vardict.vcf"})

    r.run("lofreq", "6c LoFreq", {
        "reads": r.need("dna_clipped.bam"),
        "reference_source": "cached",
        "reference_source|ref": "hg38",
        "regions": "regions_from_file",
        "regions|bed": r.need("panel_dna.bed"),
        "variant_types": "--call-indels",
        "call_control": "no",
    }, {"variants": "lofreq.vcf"})


def rna_arm(r, alias):
    if not r.find("gencode.gtf.gz") and not r.dry:
        print("fetching the GENCODE annotation into the history ...", flush=True)
        r.gi.tools.put_url(GTF_URL, r.hid, file_name="gencode.gtf.gz", file_type="gtf.gz")
        time.sleep(20); r.refresh()

    r.run("umi_extract", "1 UMI extract (RNA)", {
        "input_type_cond": "paired",
        "input_type_cond|input_read1": r.need(f"{alias}_RNA_R1.fastq.gz"),
        "input_type_cond|input_read2": r.need(f"{alias}_RNA_R2.fastq.gz"),
        "input_type_cond|bc_pattern": RNA_UMI,
        "extract_method_cond": "regex",
        "log": "true",
    }, {"out": "rna_umi_R1.fastq.gz", "out2": "rna_umi_R2.fastq.gz", "out_log": "rna_umi_extract.log"})

    r.run("rna_star", "2 STAR (arriba preset)", {
        "singlePaired": "paired",
        "singlePaired|input1": r.need("rna_umi_R1.fastq.gz"),
        "singlePaired|input2": r.need("rna_umi_R2.fastq.gz"),
        "refGenomeSource": "indexed",
        "refGenomeSource|GTFconditional": "without-gtf",
        "refGenomeSource|GTFconditional|genomeDir": "hg38",
        "refGenomeSource|GTFconditional|sjdbGTFfile": r.need("gencode.gtf.gz"),
        "refGenomeSource|GTFconditional|sjdbOverhang": "150",
        "chimOutType": "WithinBAM SoftClip",
        "algo|params": "arriba",               # STAR preset that sets every chimeric parameter Arriba needs
        "quantmode_output": "GeneCounts",
    }, {"mapped_reads": "rna_star.bam", "splice_junctions": "SJ.out.tab",
        "reads_per_gene": "rna_gene_counts.tsv", "output_log": "rna_star.log"})

    r.run("arriba_filters", "3 Arriba filter files", {
        "arriba_reference_name": "hg38_GRCh38",
    }, {"blacklist": "arriba_blacklist.tsv.gz", "known_fusions": "arriba_known_fusions.tsv.gz",
        "protein_domains": "arriba_protein_domains.gff3", "cytobands": "arriba_cytobands.tsv"})

    r.run("arriba", "4 Arriba", {
        "input": r.need("rna_star.bam"),
        "genome": "cached",
        "genome|ref_file": "hg38",
        "genome_gtf": "history",
        "genome_gtf|annotation": r.need("gencode.gtf.gz"),
        "blacklist": r.need("arriba_blacklist.tsv.gz"),
        "known_fusions": r.need("arriba_known_fusions.tsv.gz"),
        "protein_domains": r.need("arriba_protein_domains.gff3"),
        "wgs_cond": "no",
    }, {"fusions_tsv": "fusions.tsv", "discarded_fusions_tsv": "fusions.discarded.tsv"})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--server", default=os.environ.get("GALAXY_SERVER", "https://usegalaxy.eu"))
    ap.add_argument("--history", required=True, help="history name or id")
    ap.add_argument("--arm", choices=["dna", "rna", "both"], default="both")
    ap.add_argument("--alias", default=None, help="dataset name prefix (default: sample_id from config)")
    ap.add_argument("--dry-run", action="store_true", help="validate every parameter, submit nothing")
    a = ap.parse_args()
    key = os.environ.get("GALAXY_API_KEY", "")
    if not key:
        sys.exit("export GALAXY_API_KEY=<your key> first")
    gi = GalaxyInstance(url=a.server, key=key)
    hs = [h for h in gi.histories.get_histories() if a.history in (h["name"], h["id"])]
    if not hs:
        sys.exit(f"no history called {a.history!r}. Available: "
                 + ", ".join(h["name"] for h in gi.histories.get_histories()))
    hid = hs[0]["id"]
    import yaml
    alias = a.alias or yaml.safe_load(open(REPO / "config" / "config.yaml")).get("sample_id", "TUMOR01")
    print(f"history {hs[0]['name']!r} ({hid}); dataset prefix {alias}")
    r = Runner(gi, hid, dry_run=a.dry_run)
    if a.arm in ("dna", "both"):
        dna_arm(r, alias)
    if a.arm in ("rna", "both"):
        rna_arm(r, alias)
    print("\ndone. Next: python3 scripts/import_galaxy.py --history "
          f"'{hs[0]['name']}' --out results/galaxy_import")
    return 0


if __name__ == "__main__":
    sys.exit(main())
