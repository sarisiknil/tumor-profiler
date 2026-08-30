"""Educational tumour-profiling dashboard.

Reads ONLY the JSON/TSV files a pipeline run produced (results/ or results_example/). No patient identifiers,
no raw reads, no BAMs — so the dashboard can be demonstrated without exposing anything sensitive.

Run:  streamlit run dashboard/app.py -- --results results_example
"""
import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "dashboard"))
from lib import load_results, header, kpi_row, teaching_note, DISCLAIMER  # noqa: E402

st.set_page_config(page_title="Tumour profiling — educational dashboard", page_icon="🧬", layout="wide")

R = load_results()
header("Tumour DNA–RNA profiling", R)

st.markdown(
    """
This dashboard walks through what a tumour biopsy's **DNA** and **RNA** can and cannot tell us about a patient's
disease. Each page mirrors one stage of the pipeline and states, in plain language, *what the step measures*,
*what the result means*, and *how it can mislead you*.
"""
)

summary = R.get("summary", {})
var = summary.get("variants", {})
bm = summary.get("biomarkers") or {}
fus = (summary.get("fusions") or {}).get("table") or []
paths = (summary.get("pathways") or {}).get("table") or []

kpi_row([
    ("Somatic candidates", (var.get("counts_by_class") or {}).get("SOMATIC_LIKELY", 0),
     "Variants that survived tumour-only filtering"),
    ("Tier I / II", sum((var.get("counts_by_tier") or {}).get(t, 0) for t in ("I", "II")),
     "Clinically significant under AMP/ASCO/CAP"),
    ("Fusions", len(fus), "Gene fusions reported by the RNA arm"),
    ("Pathways hit", sum(1 for p in paths if int(p.get("n_altered", 0) or 0)),
     "Canonical oncogenic pathways with an altered gene"),
    ("TMB (mut/Mb)", (bm.get("tmb") or {}).get("tmb_per_mb", "n/a"),
     "Tumour mutational burden — read the caveats"),
])

st.divider()
c1, c2 = st.columns([3, 2])
with c1:
    st.subheader("The question this project asks")
    st.markdown(
        """
> *What can we understand about a person's disease by looking at their tumour's genetic material?*

**From DNA** we learn which genes carry mutations, how large a fraction of the tumour carries each one
(the variant allele frequency, a clue to clonal structure), whether a mutation is a known driver, and whether a
drug targets it. **From RNA** we learn which fusion genes exist — events that DNA panels frequently miss because
the DNA breakpoints sit in large introns — and whether exons are being skipped.

What we **cannot** learn from this sample: anything requiring a matched normal (a definitive somatic/germline
split), genome-wide mutational processes (the panel is too small), or the tumour's tissue of origin.
        """
    )
    teaching_note(
        "Why tumour-only analysis is hard",
        "Without a blood or normal-tissue sample from the same patient, every variant could in principle be "
        "inherited. We filter on population frequency (gnomAD), ClinVar and allele fraction, but rare private "
        "germline variants cannot be removed this way (Sun et al. 2018). Variants that stay ambiguous are "
        "flagged rather than deleted — and a heterozygous pathogenic variant in a cancer-predisposition gene is "
        "an ethical matter, not just a technical one.",
    )
with c2:
    st.subheader("Pipeline stages")
    st.markdown(
        """
| Stage | Page |
|---|---|
| Library structure, UMIs, primers | **1 QC & Library** |
| Alignment, deduplication, coverage | **1 QC & Library** |
| Somatic calling and filtering | **2 Variants** |
| Fusions and splicing | **3 RNA & Fusions** |
| Pathway interpretation | **4 Pathways** |
| Therapies and trials | **5 Therapies** |
| Method, limits, provenance | **6 Method & Limits** |
        """
    )
    st.caption(DISCLAIMER)

if not summary:
    st.warning("No results found. Run `snakemake -c4` first, then start the dashboard with "
               "`streamlit run dashboard/app.py -- --results results_example`.")
