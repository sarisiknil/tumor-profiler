#!/usr/bin/env python3
"""Tumour-only somatic filtering: separate likely-somatic from likely-germline and from artefacts WITHOUT a
matched normal, and label (never silently delete) every call.

Rules implemented (each one is reported per variant so the reasoning is inspectable/teachable):
  QUALITY     depth >= min_depth, alt reads >= min_alt_reads, VAF >= min_vaf, >= min_callers callers agree
  GERMLINE    gnomAD popmax AF > germline_popmax_af  (PCGR's default tumour-only rule, 0.1 %)
              OR VAF inside the heterozygous window AND present in dbSNP/gnomAD AND benign in ClinVar
  ARTEFACT    FFPE-type C>T / G>A at low VAF (flag only), low-complexity/homopolymer context if provided,
              single-caller calls at VAF < 5 %
Caveat that MUST be reported (Sun et al. 2018, GATK FAQ): population-frequency filtering cannot remove RARE,
PRIVATE germline variants; without a matched normal a residual germline fraction always remains, and
incidental germline findings raise consent duties.

Input : merged variants TSV (merge_callers.py) joined with the VEP annotation TSV (annotate_vep_rest.py)
Output: <out>.tsv with columns class (SOMATIC_LIKELY / GERMLINE_LIKELY / ARTEFACT_LIKELY / LOW_QUALITY),
        reasons, plus <out>_summary.json counts.
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json, load_config

def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

BENIGN = ("benign", "likely_benign")
PATHOGENIC = ("pathogenic", "likely_pathogenic")
# Genes where a heterozygous pathogenic variant is plausibly an inherited (germline) finding. A tumour-only
# assay cannot tell germline from somatic here; flagging them is an ethical requirement, because an incidental
# germline finding has consequences for the patient's relatives and needs genetic counselling, not a pipeline.
PREDISPOSITION = {"BRCA1","BRCA2","PALB2","ATM","CHEK2","TP53","PTEN","STK11","CDH1","MLH1","MSH2","MSH6",
                  "PMS2","EPCAM","APC","MUTYH","RB1","NF1","NF2","VHL","RET","MEN1","SDHA","SDHB","SDHC",
                  "SDHD","BAP1","CDKN2A","BRIP1","RAD51C","RAD51D","BARD1","FH","TSC1","TSC2","SMARCA4"}

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True, help="merged callers TSV")
    ap.add_argument("--vep", required=True, help="VEP annotation TSV")
    ap.add_argument("--config", default=None)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    cfg = load_config(a.config)["thresholds"]
    ann = {r["key"]: r for r in read_tsv(a.vep)}
    rows = []
    for v in read_tsv(a.variants):
        an = ann.get(v["key"], {})
        def num(x, d=0.0):
            try: return float(x)
            except (TypeError, ValueError): return d
        vaf, dp, alt_reads = num(v.get("vaf")), num(v.get("depth")), num(v.get("alt_reads"))
        ncall = int(num(v.get("n_callers"), 0))
        af_pop = num(an.get("gnomad_af_max"))
        clin = (an.get("clinvar_sig") or "").lower()
        ref, alt = v["ref"], v["alt"]
        reasons, cls = [], "SOMATIC_LIKELY"
        # quality gate
        if dp and dp < cfg["min_alt_reads"] * 5:
            reasons.append(f"low_depth({int(dp)})")
        if alt_reads and alt_reads < cfg["min_alt_reads"]:
            reasons.append(f"few_alt_reads({int(alt_reads)})")
        if vaf and vaf < cfg["min_vaf"]:
            reasons.append(f"low_vaf({vaf:.3f})")
        if ncall < cfg["min_callers_agree"]:
            reasons.append(f"single_caller({v.get('callers','')})")
        if any(r.startswith(("low_depth", "few_alt_reads", "low_vaf")) for r in reasons) or \
           (ncall < cfg["min_callers_agree"] and vaf and vaf < 0.05):
            cls = "LOW_QUALITY"
        # germline
        lo, hi = cfg["het_germline_vaf_window"]
        if af_pop > cfg["germline_popmax_af"]:
            cls = "GERMLINE_LIKELY"; reasons.append(f"gnomad_af={af_pop:g}")
        elif (vaf and lo <= vaf <= hi and an.get("dbsnp")
              and any(b in clin for b in BENIGN) and not any(p_ in clin for p_ in PATHOGENIC)):
            # only when ClinVar is *unambiguously* benign: VEP aggregates the significance of every ClinVar
            # record at the position, so a mixed string such as "benign,pathogenic" must not trigger this rule
            cls = "GERMLINE_LIKELY"; reasons.append("het_vaf+dbsnp+clinvar_benign_only")
        elif vaf and 0.90 <= vaf <= 1.0 and af_pop > 0:
            cls = "GERMLINE_LIKELY"; reasons.append("homozygous_common_variant")
        # FFPE / orientation artefact flag (never a hard filter for actionable variants: Diossy 2021)
        ffpe = (ref, alt) in (("C", "T"), ("G", "A"))
        if ffpe and vaf and vaf < 0.10:
            reasons.append("possible_FFPE_deamination(C>T/G>A, VAF<10%)")
            if cls == "SOMATIC_LIKELY" and ncall < cfg["min_callers_agree"]:
                cls = "ARTEFACT_LIKELY"
        # ambiguity flag: heterozygous-looking, rare, and either ClinVar-pathogenic or truncating in a known
        # cancer-predisposition gene -> could be an inherited variant; a tumour-only assay cannot decide.
        ambiguous = bool(vaf and lo <= vaf <= hi and af_pop <= cfg["germline_popmax_af"]
                         and (any(p_ in clin for p_ in PATHOGENIC) or an.get("gene") in PREDISPOSITION)
                         and cls == "SOMATIC_LIKELY")
        if ambiguous:
            reasons.append("germline_ambiguous: heterozygous VAF in a cancer-predisposition gene - "
                           "a matched normal would be needed; possible incidental germline finding")
        rows.append({**v, "germline_ambiguous": "yes" if ambiguous else "no",
                     "gene": an.get("gene", ""), "consequence": an.get("consequence", ""),
                     "hgvsp": an.get("hgvsp", ""), "hgvsc": an.get("hgvsc", ""),
                     "impact": an.get("impact", ""), "gnomad_af_max": an.get("gnomad_af_max", ""),
                     "clinvar_sig": an.get("clinvar_sig", ""),
                     "class": cls, "reasons": ";".join(reasons) or "passes_all_filters",
                     "ffpe_context": "yes" if ffpe else "no"})
    order = {"SOMATIC_LIKELY": 0, "ARTEFACT_LIKELY": 1, "GERMLINE_LIKELY": 2, "LOW_QUALITY": 3}
    rows.sort(key=lambda r: (order[r["class"]], -float(r.get("vaf") or 0)))
    cols = ["key","gene","hgvsc","hgvsp","consequence","impact","class","germline_ambiguous","reasons","vaf",
            "depth","alt_reads","n_callers","callers","gnomad_af_max","clinvar_sig","ffpe_context",
            "chrom","pos","ref","alt"]
    write_tsv(rows, a.out + ".tsv", cols)
    counts = {}
    for r in rows: counts[r["class"]] = counts.get(r["class"], 0) + 1
    write_json({"counts": counts, "thresholds": cfg,
                "germline_ambiguous": sum(1 for r in rows if r["germline_ambiguous"] == "yes"),
                "caveat": "Population-frequency filtering cannot exclude rare private germline variants "
                          "(Sun et al. 2018; GATK Mutect2 FAQ). Residual germline calls are possible; "
                          "germline-looking variants are labelled, not deleted."},
               a.out + "_summary.json")
    print(f"{len(rows)} variants -> {counts}", file=sys.stderr)

if __name__ == "__main__":
    main()
