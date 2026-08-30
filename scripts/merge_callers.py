#!/usr/bin/env python3
"""Merge VCFs from several somatic callers into one table with a per-variant caller-agreement count.

Why: on tumour-only amplicon data no single caller is reliable; the 2026 tumour-only amplicon benchmark
(bioRxiv 2026.04.08.717310) and Lai et al. 2016 recommend multi-caller concordance with flexible depth/VAF filters.

Usage: merge_callers.py --vcf mutect2=m2.vcf.gz --vcf vardict=vd.vcf --vcf lofreq=lf.vcf -o results/dna/variants_merged.tsv
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_vcf, vaf_and_depth, variant_key, write_tsv, write_json

def norm(rec):
    """Left-trim shared prefix so callers that pad indels differently still match."""
    chrom, pos, ref, alt = rec["chrom"], rec["pos"], rec["ref"], rec["alt"].split(",")[0]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if not chrom.startswith("chr"):
        chrom = "chr" + chrom
    return f"{chrom}:{pos}:{ref}>{alt}", chrom, pos, ref, alt

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", action="append", required=True, metavar="NAME=PATH")
    ap.add_argument("--pass-only", action="store_true", help="keep only records with FILTER PASS/.")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    merged = {}
    per_caller = {}
    for spec in a.vcf:
        name, _, path = spec.partition("=")
        n = 0
        for rec in read_vcf(path):
            if a.pass_only and rec["filter"] not in ("PASS", ".", ""):
                continue
            key, chrom, pos, ref, alt = norm(rec)
            vaf, dp, alt_reads = vaf_and_depth(rec)
            e = merged.setdefault(key, {"key": key, "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
                                        "callers": [], "vaf": {}, "depth": {}, "alt_reads": {}, "filters": {}})
            e["callers"].append(name)
            if vaf is not None: e["vaf"][name] = round(vaf, 4)
            if dp is not None: e["depth"][name] = dp
            if alt_reads is not None: e["alt_reads"][name] = alt_reads
            e["filters"][name] = rec["filter"]
            n += 1
        per_caller[name] = n
        print(f"{name}: {n} records", file=sys.stderr)
    rows = []
    for e in merged.values():
        vafs = [v for v in e["vaf"].values()]
        dps = [v for v in e["depth"].values()]
        rows.append({"key": e["key"], "chrom": e["chrom"], "pos": e["pos"], "ref": e["ref"], "alt": e["alt"],
                     "n_callers": len(set(e["callers"])), "callers": ",".join(sorted(set(e["callers"]))),
                     "vaf": round(sum(vafs)/len(vafs), 4) if vafs else "",
                     "vaf_per_caller": ";".join(f"{k}={v}" for k, v in sorted(e["vaf"].items())),
                     "depth": max(dps) if dps else "",
                     "alt_reads": max(e["alt_reads"].values()) if e["alt_reads"] else "",
                     "filters": ";".join(f"{k}={v}" for k, v in sorted(e["filters"].items()))})
    rows.sort(key=lambda r: (-r["n_callers"], r["chrom"], r["pos"]))
    write_tsv(rows, a.out, ["key","chrom","pos","ref","alt","n_callers","callers","vaf","vaf_per_caller","depth","alt_reads","filters"])
    write_json({"per_caller_records": per_caller, "merged_variants": len(rows),
                "by_agreement": {str(k): sum(1 for r in rows if r["n_callers"] == k) for k in sorted({r["n_callers"] for r in rows})}},
               str(a.out).replace(".tsv", "_summary.json"))
    print(f"merged {len(rows)} variants -> {a.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
