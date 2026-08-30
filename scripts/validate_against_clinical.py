#!/usr/bin/env python3
"""Compare this pipeline's calls against the accredited laboratory's report — the project's key result.

The laboratory used a different aligner, a different reference build (hg19), the vendor's closed analysis
software and its own reporting thresholds. Agreement is therefore evidence that an open, reproducible pipeline
reaches the same clinical conclusion; disagreement is equally informative, and is classified rather than hidden:

  CONCORDANT             reported by the laboratory and found here
  MISSED                 reported by the laboratory, not found here            <- the serious failure mode
  ADDITIONAL_SOMATIC     found here, not in the report (often below the laboratory's 5 % reporting threshold,
                         or a Tier III/IV variant the laboratory does not report at all)
  ADDITIONAL_FILTERED    found here but classified germline/artefact/low-quality

The laboratory's file is private and lives outside the repository; only the summary counts and anonymised
comparison are meant for the report.

Usage:
  python3 scripts/validate_against_clinical.py --reported ../clinical/reported_findings.tsv \
      --variants results/dna/variants_tiered.tsv --all-variants results/dna/variants_filtered.tsv \
      --fusions results/rna/fusions_summary.tsv -o results/validation/concordance
"""
import argparse, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json


def read_tsv(p):
    if not p or not Path(p).exists():
        return []
    rows, hdr = [], None
    for line in open(p):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = f
            continue
        rows.append(dict(zip(hdr, f)))
    return rows


def key_of(chrom, pos, ref, alt):
    c = chrom if str(chrom).startswith("chr") else f"chr{chrom}"
    return f"{c}:{pos}:{ref}>{alt}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reported", required=True, help="the laboratory's findings (private TSV)")
    ap.add_argument("--variants", required=True, help="our tiered variants")
    ap.add_argument("--all-variants", help="our full classified table, to explain misses")
    ap.add_argument("--fusions", help="our fusion summary")
    ap.add_argument("--vaf-tolerance", type=float, default=0.10)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    reported = read_tsv(a.reported)
    ours = read_tsv(a.variants)
    every = read_tsv(a.all_variants) if a.all_variants else ours
    fusions = read_tsv(a.fusions) if a.fusions else []

    by_key = {r.get("key"): r for r in ours}
    all_by_key = {r.get("key"): r for r in every}
    rows, missed, concordant = [], 0, 0

    for r in reported:
        if (r.get("type") or "SNV").upper() in ("SNV", "INDEL"):
            k = key_of(r["chrom_grch38"], r["pos_grch38"], r["ref"], r["alt"])
            mine = by_key.get(k) or all_by_key.get(k)
            if mine:
                concordant += 1
                try:
                    dv = abs(float(mine.get("vaf") or 0) - float(r.get("vaf") or 0))
                except ValueError:
                    dv = None
                rows.append({"status": "CONCORDANT", "gene": r["gene"], "variant": r.get("hgvs_p") or k,
                             "key": k, "lab_vaf": r.get("vaf"), "our_vaf": mine.get("vaf"),
                             "vaf_difference": round(dv, 3) if dv is not None else "",
                             "lab_tier": r.get("tier"), "our_tier": mine.get("amp_tier", ""),
                             "our_class": mine.get("class", ""), "our_callers": mine.get("callers", ""),
                             "note": ("VAF agrees within tolerance"
                                      if dv is not None and dv <= a.vaf_tolerance else
                                      "VAF differs by more than the tolerance — check tumour content and "
                                      "deduplication")})
            else:
                missed += 1
                rows.append({"status": "MISSED", "gene": r["gene"], "variant": r.get("hgvs_p") or k, "key": k,
                             "lab_vaf": r.get("vaf"), "our_vaf": "", "vaf_difference": "",
                             "lab_tier": r.get("tier"), "our_tier": "", "our_class": "", "our_callers": "",
                             "note": "not present in our call set — check coverage at this position, the "
                                     "calling intervals, and whether primer clipping removed the site"})
        elif (r.get("type") or "").upper() == "FUSION":
            want = f"{r.get('gene')}"
            hit = [f for f in fusions if want in (f.get("fusion") or "")]
            rows.append({"status": "CONCORDANT" if hit else "MISSED", "gene": want,
                         "variant": r.get("note", ""), "key": "", "lab_vaf": r.get("vaf"),
                         "our_vaf": hit[0].get("support_total") if hit else "", "vaf_difference": "",
                         "lab_tier": r.get("tier"), "our_tier": "", "our_class": "fusion",
                         "our_callers": "arriba" if hit else "",
                         "note": "fusion detected" if hit else "fusion not detected by Arriba"})
            concordant += bool(hit); missed += (not hit)

    reported_keys = {key_of(r["chrom_grch38"], r["pos_grch38"], r["ref"], r["alt"])
                     for r in reported if r.get("chrom_grch38")}
    extra_somatic = [v for v in ours if v.get("key") not in reported_keys]
    extra_filtered = [v for v in every
                      if v.get("key") not in reported_keys and v.get("class") != "SOMATIC_LIKELY"]
    for v in extra_somatic:
        rows.append({"status": "ADDITIONAL_SOMATIC", "gene": v.get("gene", ""),
                     "variant": (v.get("hgvsp") or "").split(":")[-1], "key": v.get("key", ""),
                     "lab_vaf": "", "our_vaf": v.get("vaf", ""), "vaf_difference": "",
                     "lab_tier": "", "our_tier": v.get("amp_tier", ""), "our_class": "SOMATIC_LIKELY",
                     "our_callers": v.get("callers", ""),
                     "note": "not in the laboratory's report — below its 5 % reporting threshold, a tier it "
                             "does not report, or a false positive of ours"})

    write_tsv(rows, a.out + ".tsv",
              ["status", "gene", "variant", "key", "lab_vaf", "our_vaf", "vaf_difference", "lab_tier",
               "our_tier", "our_class", "our_callers", "note"])
    n_rep = len(reported)
    write_json({"reported_by_laboratory": n_rep, "concordant": concordant, "missed": missed,
                "sensitivity_vs_report": round(concordant / n_rep, 3) if n_rep else None,
                "additional_somatic_calls": len(extra_somatic),
                "additional_calls_filtered_out": len(extra_filtered),
                "vaf_tolerance": a.vaf_tolerance,
                "interpretation": [
                    "Sensitivity is measured against the laboratory's *reported* variants, not against ground "
                    "truth: the laboratory reports only what passes its clinical thresholds, so 'additional' "
                    "calls here are not necessarily false positives.",
                    "The laboratory aligned to hg19 with the vendor's closed software; coordinates were "
                    "converted to GRCh38 before comparison.",
                    "A missed variant is the failure that matters — check coverage at the position, the "
                    "calling intervals and the primer-clipping step before concluding anything else."]},
               a.out + "_summary.json")
    print(f"reported {n_rep} | concordant {concordant} | missed {missed} | "
          f"additional somatic {len(extra_somatic)}", file=sys.stderr)
    for r in rows[:20]:
        print(f"  {r['status']:20s} {r['gene']:8s} {r['variant']:16s} lab={r['lab_vaf']} ours={r['our_vaf']}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
