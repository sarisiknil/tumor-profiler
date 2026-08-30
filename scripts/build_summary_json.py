#!/usr/bin/env python3
"""Collect every stage's output into one machine-readable summary (results/summary.json) plus a human-readable
Markdown report (results/report.md). The Streamlit dashboard and the internship report both read this file, so
there is a single source of truth for every number that gets quoted.
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_json, load_config

def jload(p, default=None):
    p = Path(p)
    if not p.exists():
        return default
    try:
        return json.loads(p.read_text())
    except Exception:
        return default

def tload(p, limit=None):
    p = Path(p)
    if not p.exists():
        return []
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]
    return rows[:limit] if limit else rows

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default="results", help="results directory")
    ap.add_argument("--config", default=None)
    ap.add_argument("--sample", default=None, help="public alias for the sample")
    a = ap.parse_args()
    R = Path(a.results)
    cfg = load_config(a.config)
    alias = a.sample or cfg.get("sample_id", "TUMOR01")

    tiers = tload(R / "dna" / "variants_tiered.tsv")
    filtered = tload(R / "dna" / "variants_filtered.tsv")
    fusions = tload(R / "rna" / "fusions_summary.tsv")
    pathways = tload(R / "pathways" / "pathway_hits.tsv")
    evidence = tload(R / "evidence" / "annot_evidence.tsv")
    trials = tload(R / "evidence" / "annot_trials.tsv")
    summary = {
        "sample": alias,
        "library": cfg.get("library", {}),
        "panel": cfg.get("panel", {}),
        "qc": {"dna": jload(R / "qc" / "read_structure_DNA.json"),
               "rna": jload(R / "qc" / "read_structure_RNA.json")},
        "variants": {
            "counts_by_class": (jload(R / "dna" / "variants_filtered_summary.json") or {}).get("counts", {}),
            "counts_by_tier": (jload(R / "dna" / "variants_tiered_summary.json") or {}).get("tier_counts", {}),
            "caller_concordance": jload(R / "dna" / "variants_merged_summary.json"),
            "tier_I_II": [v for v in tiers if v.get("amp_tier") in ("I", "II")],
            "top_somatic": [v for v in filtered if v.get("class") == "SOMATIC_LIKELY"][:50],
        },
        "biomarkers": jload(R / "dna" / "biomarkers.json"),
        "copy_number": jload(R / "dna" / "cnv_summary.json"),
        "fusions": {"summary": jload(R / "rna" / "fusions_summary_summary.json"), "table": fusions},
        "exon_skipping": jload(R / "rna" / "exon_skipping.json"),
        "pathways": {"table": pathways, "summary": jload(R / "pathways" / "pathway_hits_summary.json")},
        "evidence": {"n_items": len(evidence), "items": evidence[:200], "trials": trials[:50]},
        "validation": jload(R / "validation" / "concordance.json"),
        "disclaimer": ("Educational output of a student pipeline. Not a clinical report and not a basis for "
                       "medical decisions. Somatic findings only; no identifying information."),
    }
    write_json(summary, R / "summary.json")

    # ---- Markdown report ----------------------------------------------------
    L = [f"# Tumour profiling summary — sample `{alias}`", "",
         "> Educational output. Not a clinical report.", "",
         "## 1. Library and QC", ""]
    for arm in ("dna", "rna"):
        q = (summary["qc"].get(arm) or {}).get("R1")
        if q:
            L.append(f"- **{arm.upper()}**: {q['reads_analyzed']:,} reads inspected on `{q['instrument']}`; "
                     f"UMI {q['umi_len_inferred']} nt + common region `{q['anchor_inferred']}`; "
                     f"adapter read-through {q['adapter']['frac_reads_with_adapter']*100:.1f} %")
    bm = summary.get("biomarkers") or {}
    if bm.get("tmb"):
        t = bm["tmb"]
        L += ["", "## 2. Biomarkers", "",
              f"- **TMB**: {t['tmb_per_mb']} mutations/Mb "
              f"({t['nonsynonymous_somatic_variants']} nonsynonymous over {t['assessable_mb']} Mb assessable; "
              f"95 % CI {t['ci95_per_mb'][0]}–{t['ci95_per_mb'][1]})"
              + ("" if t.get("meaningful", True) else
                 "  \n  **This TMB is not interpretable**: the assessable territory is far below the ~1 Mb "
                 "that panel TMB estimation assumes. The value is shown only to demonstrate the calculation."),
              f"- FFPE indicator (C>T/G>A share of SNVs): {(bm.get('ffpe_indicator') or {}).get('C>T_or_G>A_fraction_of_snvs')}"]
    L += ["", "## 3. Variants", ""]
    cbc = summary["variants"]["counts_by_class"]; cbt = summary["variants"]["counts_by_tier"]
    if cbc: L.append("- Classification: " + ", ".join(f"{k} {v}" for k, v in cbc.items()))
    if cbt: L.append("- AMP/ASCO/CAP tiers: " + ", ".join(f"Tier {k}: {v}" for k, v in sorted(cbt.items())))
    if summary["variants"]["tier_I_II"]:
        L += ["", "| Gene | Variant | Tier | ESCAT | VAF | Therapies |", "|---|---|---|---|---|---|"]
        for v in summary["variants"]["tier_I_II"][:20]:
            L.append(f"| {v.get('gene','')} | {v.get('hgvsp') or v.get('hgvsc','')} | {v.get('amp_tier','')} | "
                     f"{v.get('escat','')} | {v.get('vaf','')} | {(v.get('therapies') or '')[:60]} |")
    if fusions:
        L += ["", "## 4. Fusions and splicing", "", "| Fusion | Confidence | Frame | Support | Targetable |",
              "|---|---|---|---|---|"]
        for f in fusions[:15]:
            L.append(f"| {f.get('fusion','')} | {f.get('confidence','')} | {f.get('reading_frame','')} | "
                     f"{f.get('support_total','')} | {f.get('targetable','')} |")
    es = summary.get("exon_skipping") or {}
    for e in (es.get("events") or []):
        if e.get("called") == "yes":
            L.append(f"- **{e['event']}** detected: {e['skipping_reads']} skipping reads "
                     f"(ratio {e['skipping_ratio']}) — {e['significance']}")
    if pathways:
        L += ["", "## 5. Pathways", ""]
        for p in pathways:
            if int(p.get("n_altered", 0) or 0):
                L.append(f"- **{p['pathway']}**: {p['n_altered']} altered gene(s) — {p['altered_genes']}")
    if trials:
        L += ["", "## 6. Clinical-trial context (recruiting)", ""]
        for t in trials[:10]:
            L.append(f"- [{t.get('nct_id','')}]({t.get('url','')}) — {t.get('title','')} ({t.get('phase','')})")
    L += ["", "## 7. Caveats", ""]
    for c in ((bm.get("tmb") or {}).get("caveats") or []):
        L.append(f"- {c}")
    for c in ((summary["fusions"]["summary"] or {}).get("caveats") or []):
        L.append(f"- {c}")
    L += ["", f"_Generated by tumor-profiler; see results/summary.json for the machine-readable version._"]
    (R / "report.md").write_text("\n".join(L) + "\n")
    print(f"wrote {R}/summary.json and {R}/report.md")

if __name__ == "__main__":
    main()
