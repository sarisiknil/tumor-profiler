"""Page 5 — from a variant to the treatment literature (and its limits)."""
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import load_results, header, teaching_note  # noqa: E402

st.set_page_config(page_title="Therapies", page_icon="💊", layout="wide")
R = load_results()
header("5 · Therapeutic and trial context", R)

st.error(
    "This page summarises published evidence linked to the alterations found. It is **not** a treatment "
    "recommendation. Turning a variant into a therapy decision requires the tumour type, the patient's history "
    "and prior treatments, confirmation of the finding in an accredited laboratory, and a molecular tumour "
    "board.", icon="⚕️")

ev = R.get("evidence")
if isinstance(ev, pd.DataFrame) and not ev.empty:
    st.subheader("Evidence items")
    lvl = st.multiselect("Evidence level", sorted(x for x in ev["level"].dropna().unique() if x),
                         default=sorted(x for x in ev["level"].dropna().unique() if x))
    view = ev[ev["level"].isin(lvl)] if lvl else ev
    st.dataframe(view[["gene", "alteration", "source", "type", "level", "significance", "direction",
                       "therapies", "disease", "url"]],
                 use_container_width=True, hide_index=True,
                 column_config={"url": st.column_config.LinkColumn("source", display_text="open")})
    st.caption("CIViC evidence levels run from A (validated in clinical practice) through B (clinical trial), "
               "C (case study), D (preclinical) to E (inferential). OncoKB levels 1–2 mark standard-of-care "
               "biomarkers; 3–4 are investigational. The same alteration can be level A in one tumour type and "
               "have no evidence at all in another — which is why the disease column matters as much as the "
               "level.")

    st.subheader("Same variant, different disease")
    if "disease" in ev.columns:
        cross = ev.groupby(["gene", "alteration", "disease"]).size().reset_index(name="evidence items")
        st.dataframe(cross.sort_values("evidence items", ascending=False).head(30),
                     use_container_width=True, hide_index=True)

dr = R.get("druggability")
if isinstance(dr, pd.DataFrame) and not dr.empty:
    st.subheader("Drug–gene interactions (DGIdb)")
    st.dataframe(dr.head(60), use_container_width=True, hide_index=True)
    st.caption("A drug–gene interaction is a much weaker statement than clinical evidence: it means somebody "
               "has reported that a compound acts on this gene product, not that it helps this patient.")

tr = R.get("trials")
if isinstance(tr, pd.DataFrame) and not tr.empty:
    st.subheader("Recruiting clinical trials mentioning these genes")
    st.dataframe(tr[["matched_gene", "nct_id", "title", "phase", "status", "url"]],
                 use_container_width=True, hide_index=True,
                 column_config={"url": st.column_config.LinkColumn("registry", display_text="open")})
    st.caption("Retrieved live from ClinicalTrials.gov. Eligibility depends on far more than a gene name.")

teaching_note(
    "Where the knowledge comes from, and what it costs",
    "**CIViC** is community-curated and released under CC0, so its evidence can be redistributed freely — it is "
    "the backbone of this page. **OncoKB** is free for academic research with a personal token, but its "
    "annotations may not be redistributed, so they are written only to the local results folder. **COSMIC** is "
    "free for academics yet its licence excludes for-profit settings and forbids re-publication, so this "
    "pipeline does not query it at all. Choosing the free, redistributable source wherever possible is what "
    "makes the repository shareable.",
)
