#!/usr/bin/env python3
"""Download the small result files from a Galaxy history so the local half of the pipeline can continue.

Deliberately does NOT download BAM files: they are large, they stay on the Galaxy server, and every local step
works from VCFs, coverage tables and fusion calls. Requires GALAXY_API_KEY (Galaxy → User → Preferences →
Manage API Key).

Usage:
  export GALAXY_API_KEY=...
  python3 scripts/import_galaxy.py --list
  python3 scripts/import_galaxy.py --history "Tumour profiling — DNA" --out results/galaxy_import
"""
import argparse, json, os, re, sys, urllib.parse, urllib.request
from pathlib import Path

DEFAULT_SERVER = "https://usegalaxy.eu"
# what the local pipeline needs, matched case-insensitively against dataset names
WANTED = [
    (r"mutect", "mutect2.vcf"), (r"vardict", "vardict.vcf"), (r"lofreq", "lofreq.vcf"),
    (r"freebayes", "freebayes.vcf"), (r"coverage\.regions|regions\.bed|mosdepth.*region", "coverage.regions.bed"),
    (r"fusions\.discarded|discarded", "fusions.discarded.tsv"), (r"fusions", "fusions.tsv"),
    (r"sj\.out|splice.*junction", "SJ.out.tab"), (r"gene_counts|reads_per_gene", "rna_gene_counts.tsv"),
    (r"dna_umi_extract", "dna_umi_extract.log"), (r"dna_dedup", "dna_dedup.log"),
    (r"rna_umi_extract", "rna_umi_extract.log"), (r"rna_star\.log|star.*log", "rna_star.log"),
    (r"coverage\.summary", "coverage.summary.txt"),
    (r"f1r2", "f1r2.tar.gz"),
    (r"mutect.*stats|stats.*mutect", "mutect2.vcf.stats"),
]
SKIP_EXT = {".bam", ".cram", ".sam", ".fastqsanger", ".fastqsanger.gz", ".fastq.gz"}


def api(server, key, path, params=None):
    q = dict(params or {}); q["key"] = key
    url = f"{server}/api/{path}?" + urllib.parse.urlencode(q)
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "tumor-profiler/1.0"}),
                                timeout=120) as r:
        return json.load(r)


def download(server, key, dataset_id, dest):
    url = f"{server}/api/datasets/{dataset_id}/display?to_ext=data&key={key}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "tumor-profiler/1.0"}),
                                timeout=900) as r, open(dest, "wb") as fh:
        while chunk := r.read(1 << 20):
            fh.write(chunk)
    return dest.stat().st_size


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--server", default=os.environ.get("GALAXY_SERVER", DEFAULT_SERVER))
    ap.add_argument("--history", help="history name (substring match) — omit with --list")
    ap.add_argument("--out", default="results/galaxy_import")
    ap.add_argument("--list", action="store_true", help="list histories and exit")
    ap.add_argument("--all", action="store_true", help="download every dataset, not just the wanted ones")
    a = ap.parse_args()
    key = os.environ.get("GALAXY_API_KEY", "")
    if not key:
        print("Set GALAXY_API_KEY first (Galaxy → User → Preferences → Manage API Key).", file=sys.stderr)
        return 2
    hists = api(a.server, key, "histories")
    if a.list or not a.history:
        for h in hists:
            print(f"{h['id']}  {h['name']}")
        return 0
    match = [h for h in hists if a.history.lower() in h["name"].lower() or h["id"] == a.history]
    if not match:
        print(f"no history matching {a.history!r}", file=sys.stderr); return 1
    hist = match[0]
    print(f"history: {hist['name']} ({hist['id']})")
    contents = api(a.server, key, f"histories/{hist['id']}/contents", {"v": "dev", "keys": "id,name,extension,deleted,visible,state"})
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for d in contents:
        if d.get("deleted") or not d.get("visible") or d.get("state") != "ok":
            continue
        name = d.get("name", "")
        ext = "." + (d.get("extension") or "")
        if ext in SKIP_EXT and not a.all:
            print(f"  skip (large): {name}")
            continue
        target = None
        for pattern, filename in WANTED:
            if re.search(pattern, name, re.I):
                target = filename
                break
        if target is None:
            if not a.all:
                continue
            target = re.sub(r"[^A-Za-z0-9._-]+", "_", name)
        size = download(a.server, key, d["id"], out / target)
        print(f"  {name}  ->  {target}  ({size/1e6:.1f} MB)")
        manifest.append({"galaxy_dataset": name, "galaxy_id": d["id"], "local": target, "bytes": size})
    (out / "manifest.json").write_text(json.dumps(
        {"server": a.server, "history": hist["name"], "history_id": hist["id"], "files": manifest}, indent=1))
    print(f"\n{len(manifest)} file(s) -> {out}/ (manifest.json records the Galaxy dataset ids for provenance)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
