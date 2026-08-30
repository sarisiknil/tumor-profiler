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

def parse_umi_logs(results_dir):
    """Reads in, reads surviving the UMI pattern, and unique molecules out — the library-complexity numbers.

    The fraction of reads matching the expected `UMI + common region` pattern is also an independent check on
    the read-structure inference: if the structure were wrong, most reads would fail the pattern."""
    import re
    out = {}
    ex = Path(results_dir) / "galaxy_import" / "dna_umi_extract.log"
    if ex.exists():
        txt = ex.read_text(errors="replace")
        m_in = re.search(r"Input Reads:\s*(\d+)", txt)
        m_out = re.search(r"Reads output:\s*(\d+)", txt)
        if m_in and m_out:
            a, b = int(m_in.group(1)), int(m_out.group(1))
            out["umi_pattern"] = {"reads_in": a, "reads_matching_pattern": b,
                                  "fraction_matching": round(b / a, 4) if a else None}
    rex = Path(results_dir) / "galaxy_import" / "rna_umi_extract.log"
    if rex.exists():
        txt = rex.read_text(errors="replace")
        m_in = re.search(r"Input Reads:\s*(\d+)", txt)
        m_out = re.search(r"Reads output:\s*(\d+)", txt)
        if m_in and m_out:
            a, b = int(m_in.group(1)), int(m_out.group(1))
            out["rna_umi_pattern"] = {"reads_in": a, "reads_matching_pattern": b,
                                      "fraction_matching": round(b / a, 4) if a else None}
    star = Path(results_dir) / "galaxy_import" / "rna_star.log"
    if star.exists():
        txt = star.read_text(errors="replace")

        def g(pattern, cast=float):
            m = re.search(pattern + r"\s*\|\s*([\d.]+)%?", txt)
            return cast(m.group(1)) if m else None

        out["rna_alignment"] = {
            "input_reads": g(r"Number of input reads", int),
            "uniquely_mapped_pct": g(r"Uniquely mapped reads %"),
            "multi_mapped_pct": g(r"% of reads mapped to multiple loci"),
            "unmapped_too_short_pct": g(r"% of reads unmapped: too short"),
            "chimeric_pct": g(r"% of chimeric reads"),
            "note": "A very high chimeric rate is characteristic of anchored multiplex PCR, not of the tumour: "
                    "Read 1 starts at an arbitrary ligation point and Read 2 at a fixed gene-specific primer, "
                    "and with short inserts STAR's permissive chimeric threshold (chimSegmentMin 10, set by "
                    "the Arriba preset) classifies a large share of reads as chimeric. It is the reason "
                    "whole-transcriptome fusion callers over-call on this chemistry, and the reason the fusion "
                    "summary keeps the panel-membership check.",
        }
    dd = Path(results_dir) / "galaxy_import" / "dna_dedup.log"
    if dd.exists():
        txt = dd.read_text(errors="replace")
        m_out = re.search(r"Number of reads out:\s*(\d+)", txt)
        m_pos = re.search(r"Total number of positions deduplicated:\s*(\d+)", txt)
        m_umi = re.search(r"Mean number of unique UMIs per position:\s*([\d.]+)", txt)
        m_in = re.search(r"Input Reads:\s*(\d+)", txt)
        d = {}
        if m_in: d["reads_in"] = int(m_in.group(1))
        if m_out: d["unique_molecules_out"] = int(m_out.group(1))
        if m_pos: d["positions_deduplicated"] = int(m_pos.group(1))
        if m_umi: d["mean_unique_umis_per_position"] = float(m_umi.group(1))
        if d.get("reads_in") and d.get("unique_molecules_out"):
            d["reads_per_molecule"] = round(d["reads_in"] / d["unique_molecules_out"], 2)
        if d:
            out["deduplication"] = d
    return out


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
               "rna": jload(R / "qc" / "read_structure_RNA.json"),
               "library_complexity": parse_umi_logs(R)},
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
        "validation": {"summary": jload(R / "validation" / "concordance_summary.json"),
                       "table": tload(R / "validation" / "concordance.tsv")},
        "disclaimer": ("Educational output of a student pipeline. Not a clinical report and not a basis for "
                       "medical decisions. Somatic findings only; no identifying information."),
    }
    write_json(summary, R / "summary.json")

    # ---- Markdown report ----------------------------------------------------
    section = {"n": 0}

    def head(title):
        section["n"] += 1
        return ["", f"## {section['n']}. {title}", ""]

    L = [f"# Tumour profiling summary — sample `{alias}`", "",
         "> Educational output. Not a clinical report.", ""] + head("Library and QC")
    for arm in ("dna", "rna"):
        q = (summary["qc"].get(arm) or {}).get("R1")
        if q:
            L.append(f"- **{arm.upper()}**: {q['reads_analyzed']:,} reads inspected on `{q['instrument']}`; "
                     f"UMI {q['umi_len_inferred']} nt + common region `{q['anchor_inferred']}`; "
                     f"adapter read-through {q['adapter']['frac_reads_with_adapter']*100:.1f} %")
    lc = (summary["qc"].get("library_complexity") or {})
    if lc.get("umi_pattern"):
        u = lc["umi_pattern"]
        L.append(f"- **{u['fraction_matching']*100:.1f} %** of DNA reads carry the expected "
                 f"`12-nt UMI + common region` structure ({u['reads_matching_pattern']:,} of "
                 f"{u['reads_in']:,}) — an independent confirmation of the library identification")
    if lc.get("deduplication"):
        d = lc["deduplication"]
        L.append(f"- Deduplication: {d.get('reads_in', 0):,} reads collapse to "
                 f"**{d.get('unique_molecules_out', 0):,} unique molecules** "
                 f"({d.get('reads_per_molecule', '?')} reads per molecule, "
                 f"{d.get('mean_unique_umis_per_position', '?')} molecules per start position). "
                 f"Molecules, not reads, are what the variant calls actually rest on.")
    if lc.get("rna_umi_pattern"):
        u = lc["rna_umi_pattern"]
        L.append(f"- **{u['fraction_matching']*100:.1f} %** of RNA reads carry the expected structure "
                 f"({u['reads_matching_pattern']:,} of {u['reads_in']:,})")
    if lc.get("rna_alignment"):
        ra = lc["rna_alignment"]
        L.append(f"- RNA alignment: {ra.get('uniquely_mapped_pct')} % uniquely mapped, "
                 f"**{ra.get('chimeric_pct')} % chimeric** — expected for anchored multiplex PCR, "
                 f"not a property of the tumour (see the note in summary.json)")
    bm = summary.get("biomarkers") or {}
    if bm.get("tmb"):
        t = bm["tmb"]
        L += head("Biomarkers") + [
              f"- **TMB**: {t['tmb_per_mb']} mutations/Mb "
              f"({t['nonsynonymous_somatic_variants']} nonsynonymous over {t['assessable_mb']} Mb assessable; "
              f"95 % CI {t['ci95_per_mb'][0]}–{t['ci95_per_mb'][1]})"
              + ("" if t.get("meaningful", True) else
                 "  \n  **This TMB is not interpretable**: the assessable territory is far below the ~1 Mb "
                 "that panel TMB estimation assumes. The value is shown only to demonstrate the calculation."),
              ] + ([f"- FFPE indicator (C>T/G>A share of somatic SNVs): "
                    f"{(bm.get('ffpe_indicator') or {}).get('C>T_or_G>A_fraction_of_snvs')}"]
                   if (bm.get('ffpe_indicator') or {}).get('interpretable') else
                   [f"- FFPE indicator not reported: only "
                    f"{(bm.get('ffpe_indicator') or {}).get('n_snvs', 0)} somatic SNV(s), too few to interpret"])
    L += head("Variants")
    cbc = summary["variants"]["counts_by_class"]; cbt = summary["variants"]["counts_by_tier"]
    if cbc: L.append("- Classification: " + ", ".join(f"{k} {v}" for k, v in cbc.items()))
    if cbt: L.append("- AMP/ASCO/CAP tiers: " + ", ".join(f"Tier {k}: {v}" for k, v in sorted(cbt.items())))
    if summary["variants"]["tier_I_II"]:
        L += ["", "| Gene | Variant | Tier | ESCAT | VAF | Therapies |", "|---|---|---|---|---|---|"]
        for v in summary["variants"]["tier_I_II"][:20]:
            L.append(f"| {v.get('gene','')} | {v.get('hgvsp') or v.get('hgvsc','')} | {v.get('amp_tier','')} | "
                     f"{v.get('escat','')} | {v.get('vaf','')} | {(v.get('therapies') or '')[:60]} |")
    if fusions:
        L += head("Fusions and splicing") + ["| Fusion | Confidence | Frame | Support | Targetable |",
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
        L += head("Pathways")
        for p in pathways:
            if int(p.get("n_altered", 0) or 0):
                L.append(f"- **{p['pathway']}**: {p['n_altered']} altered gene(s) — {p['altered_genes']}")
    if trials:
        L += head("Clinical-trial context (recruiting)")
        for t in trials[:10]:
            L.append(f"- [{t.get('nct_id','')}]({t.get('url','')}) — {t.get('title','')} ({t.get('phase','')})")
    val = (summary.get("validation") or {}).get("summary") or {}
    vtab = (summary.get("validation") or {}).get("table") or []
    if val:
        L += head("Validation against the accredited laboratory's report") + [
              f"- Reported by the laboratory: **{val.get('reported_by_laboratory')}** · "
              f"recovered here: **{val.get('concordant')}** · missed: **{val.get('missed')}** "
              f"(sensitivity {val.get('sensitivity_vs_report')})",
              f"- Additional somatic calls not in the report: {val.get('additional_somatic_calls')}"]
        if vtab:
            L += ["", "| Status | Gene | Variant | Lab VAF | Our VAF | Lab tier | Our tier |", "|---|---|---|---|---|---|---|"]
            for r in vtab[:15]:
                L.append(f"| {r.get('status','')} | {r.get('gene','')} | {r.get('variant','')} | "
                         f"{r.get('lab_vaf','')} | {r.get('our_vaf','')} | {r.get('lab_tier','')} | "
                         f"{r.get('our_tier','')} |")
        for c in (val.get("interpretation") or []):
            L.append(f"- {c}")
    L += head("Caveats")
    for c in ((bm.get("tmb") or {}).get("caveats") or []):
        L.append(f"- {c}")
    for c in ((summary["fusions"]["summary"] or {}).get("caveats") or []):
        L.append(f"- {c}")
    L += ["", f"_Generated by tumor-profiler; see results/summary.json for the machine-readable version._"]
    (R / "report.md").write_text("\n".join(L) + "\n")
    print(f"wrote {R}/summary.json and {R}/report.md")

if __name__ == "__main__":
    main()
