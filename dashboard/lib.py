"""Shared helpers for the dashboard: results loading, layout blocks, and the consistent visual language."""
import argparse, json
from pathlib import Path

import pandas as pd
import streamlit as st

REPO = Path(__file__).resolve().parents[1]
DISCLAIMER = ("Educational output of a student pipeline. Not a clinical report; not a basis for medical "
              "decisions. Somatic findings only, no identifying information.")

TIER_COLOR = {"I": "#2f6f4f", "II": "#7a5c1e", "III": "#4a5568", "IV": "#6b7280"}
CLASS_HELP = {
    "SOMATIC_LIKELY": "Passed every filter: rare in the population, adequate depth and allele fraction, "
                      "supported by more than one caller.",
    "GERMLINE_LIKELY": "Looks inherited: common in gnomAD, or a heterozygous benign dbSNP/ClinVar variant.",
    "ARTEFACT_LIKELY": "Pattern of a technical artefact — typically a low-VAF C>T/G>A change from formalin "
                       "fixation, seen by only one caller.",
    "LOW_QUALITY": "Too little evidence: low depth, few alternate reads, or a single caller at low VAF.",
}

def _args():
    p = argparse.ArgumentParser()
    p.add_argument("--results", default=None)
    a, _ = p.parse_known_args()
    return a

@st.cache_data(show_spinner=False)
def load_results(results_dir: str | None = None) -> dict:
    d = results_dir or _args().results
    if not d:
        for cand in ("results", "results_example"):
            if (REPO / cand / "summary.json").exists():
                d = cand
                break
    d = d or "results_example"
    root = REPO / d
    out = {"dir": str(root), "summary": {}}
    sj = root / "summary.json"
    if sj.exists():
        out["summary"] = json.loads(sj.read_text())
    for key, rel in {
        "variants_tiered": "dna/variants_tiered.tsv",
        "variants_filtered": "dna/variants_filtered.tsv",
        "vep": "dna/vep.tsv",
        "cnv": "dna/cnv.tsv",
        "fusions": "rna/fusions_summary.tsv",
        "exon_skipping": "rna/exon_skipping.tsv",
        "pathways": "pathways/pathway_hits.tsv",
        "evidence": "evidence/annot_evidence.tsv",
        "trials": "evidence/annot_trials.tsv",
        "druggability": "evidence/annot_druggability.tsv",
        "primers": "primers/gsp2_DNA.tsv",
    }.items():
        p = root / rel
        if p.exists() and p.stat().st_size > 0:
            try:
                out[key] = pd.read_csv(p, sep="\t")
            except Exception:
                out[key] = pd.DataFrame()
        else:
            out[key] = pd.DataFrame()
    for key, rel in {"biomarkers": "dna/biomarkers.json", "qc_dna": "qc/read_structure_DNA.json",
                     "qc_rna": "qc/read_structure_RNA.json", "cnv_summary": "dna/cnv_summary.json",
                     "fusion_summary": "rna/fusions_summary_summary.json",
                     "filter_summary": "dna/variants_filtered_summary.json",
                     "merge_summary": "dna/variants_merged_summary.json",
                     "pathway_summary": "pathways/pathway_hits_summary.json",
                     "tier_summary": "dna/variants_tiered_summary.json"}.items():
        p = root / rel
        out[key] = json.loads(p.read_text()) if p.exists() else {}
    return out

def header(title: str, R: dict):
    s = R.get("summary", {})
    st.title(f"🧬 {title}")
    left, right = st.columns([3, 1])
    with left:
        st.caption(f"Sample **{s.get('sample', '—')}** · library: "
                   f"{(s.get('library') or {}).get('chemistry', 'unknown chemistry')} · "
                   f"results from `{R['dir']}`")
    with right:
        st.caption(f"generated {(s.get('_provenance') or {}).get('generated_utc', '—')}")

def kpi_row(items):
    cols = st.columns(len(items))
    for col, (label, value, help_text) in zip(cols, items):
        col.metric(label, value, help=help_text)

def teaching_note(title: str, body: str):
    with st.expander(f"💡 {title}", expanded=False):
        st.markdown(body)

def caveats(items, title="Caveats"):
    if not items:
        return
    st.markdown(f"**{title}**")
    for c in items:
        st.markdown(f"- {c}")
