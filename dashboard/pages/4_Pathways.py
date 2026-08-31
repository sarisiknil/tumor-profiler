"""Page 4 — from a list of genes to a picture of what is broken in the cell."""
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import load_results, header, teaching_note  # noqa: E402

st.set_page_config(page_title="Pathways", page_icon="🕸️", layout="wide")
R = load_results()
header("4 · Oncogenic pathways", R)

st.markdown(
    "A mutation list is not yet an explanation. Grouping altered genes into the signalling pathways they belong "
    "to shows *which control systems of the cell have been broken* — and pathways, not individual genes, are "
    "what most targeted drugs act on. The gene sets used here are the ten canonical oncogenic pathways curated "
    "across 9,125 tumours by Sanchez-Vega et al. (2018), plus two supplementary sets those ten deliberately "
    "exclude (DNA repair and chromatin regulators)."
)

pw = R.get("pathways")
if isinstance(pw, pd.DataFrame) and not pw.empty:
    hit = pw[pw["n_altered"] > 0].copy()
    if not hit.empty:
        fig = px.bar(hit.sort_values("n_altered"), x="n_altered", y="pathway", orientation="h",
                     text="altered_genes", labels={"n_altered": "altered genes", "pathway": ""})
        fig.update_traces(marker_color="#3f6fa8", textposition="outside", cliponaxis=False)
        fig.update_layout(height=110 + 40 * len(hit), margin=dict(l=10, r=120, t=30, b=10))
        st.plotly_chart(fig, width='stretch')
        for _, r in hit.iterrows():
            with st.expander(f"{r['pathway']} — {r['n_altered']} of {r['n_genes_in_set']} genes altered"):
                st.markdown(r["alterations"] or "—")
                st.caption(r["source"])
    st.dataframe(pw[["pathway", "n_altered", "n_genes_in_set", "altered_genes"]],
                 width='stretch', hide_index=True)

teaching_note(
    "Read this as a map, not as a statistic",
    "This panel is **not** an enrichment analysis, and it would be wrong to present it as one. Enrichment asks "
    "whether a pathway contains more altered genes than chance would predict — but the gene universe here was "
    "chosen by whoever designed the panel, and a single sample provides no null distribution. What the figure "
    "does show is honest and useful: *which* pathways carry an alteration. Just as important is what it cannot "
    "show — a pathway with no hit may simply have no genes on the panel.",
)

ps = R.get("pathway_summary") or {}
if ps.get("unassigned_genes"):
    st.info("Altered genes not in any curated set: " + ", ".join(ps["unassigned_genes"]))
