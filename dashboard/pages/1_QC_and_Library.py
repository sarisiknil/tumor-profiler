"""Page 1 — what kind of experiment produced these reads, and can we trust them?"""
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import load_results, header, teaching_note  # noqa: E402

st.set_page_config(page_title="QC & Library", page_icon="🔬", layout="wide")
R = load_results()
header("1 · Quality control and library structure", R)

st.markdown(
    "Before any biology, the data must be identified: **which assay produced it**, whether the reads carry "
    "molecular barcodes, where the primers sit, and how much of the target region was actually sequenced deeply "
    "enough to call variants."
)

for arm, key in (("DNA", "qc_dna"), ("RNA", "qc_rna")):
    q = (R.get(key) or {}).get("R1")
    q2 = (R.get(key) or {}).get("R2")
    if not q:
        continue
    st.subheader(f"{arm} library")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Instrument", q.get("instrument", "—"))
    c2.metric("Read length (mode)", q["read_length"]["mode"])
    c3.metric("UMI length", q["umi_len_inferred"] or "none")
    c4.metric("Adapter read-through", f"{q['adapter']['frac_reads_with_adapter']*100:.1f} %",
              help="Reads that ran off the end of a short fragment into the adapter. High values indicate "
                   "short inserts, typical of degraded (e.g. formalin-fixed) material.")
    st.markdown(
        f"**Read 1 structure** &nbsp; `[{q['umi_len_inferred']}-nt UMI]"
        f"[{q['anchor_inferred']}{q['anchor_tail_inferred']}]"
        f"[insert from base {q['insert_starts_at_base']}]`"
    )
    if q2:
        cov = q2.get("start_kmers_covering", {})
        st.markdown(
            f"**Read 2 structure** &nbsp; starts at a gene-specific primer — "
            f"{cov.get('80pct')} distinct 20-mers cover 80 % of reads out of "
            f"{q2['distinct_start_kmers']:,} observed."
        )
        # positional base composition: shows the random UMI, then the fixed common region
        pc = pd.DataFrame(q["positional_consensus"])
        fig = go.Figure()
        for base, colour in (("A", "#4f8a5b"), ("C", "#3f6fa8"), ("G", "#b5843a"), ("T", "#a34f4f")):
            fig.add_bar(x=pc["pos"], y=pc[base], name=base, marker_color=colour)
        fig.update_layout(barmode="stack", height=260, margin=dict(l=10, r=10, t=30, b=10),
                          title=f"{arm} Read 1 — base composition by position",
                          xaxis_title="position in read", yaxis_title="fraction",
                          legend=dict(orientation="h", y=1.15))
        st.plotly_chart(fig, width='stretch')
        st.caption("Positions 1–12 are random (all four bases ~25 %): that is the unique molecular identifier. "
                   "The block that follows is a single fixed base per position: the adapter's common region. "
                   "Everything after it is genomic or transcript sequence.")
    st.divider()

teaching_note(
    "Why the read structure matters so much",
    "This library was built by **anchored multiplex PCR**: one gene-specific primer defines one end of every "
    "fragment, and a ligated adapter carrying a random 12-nt barcode defines the other. Two consequences drive "
    "the whole pipeline:\n\n"
    "1. **Deduplication must use the barcode.** Every molecule from one amplicon starts at the same coordinate, "
    "so ordinary positional duplicate marking would throw away real molecules — measured on commercial panels, "
    "the unique-molecule rate rose from 10 % to 52 % when UMIs were used instead (Kim et al. 2019).\n"
    "2. **Primer bases are not evidence.** The first ~25 bases of Read 2 come from a synthetic primer, not from "
    "the patient. If they were left in, a mismatch designed into the primer would look like a mutation. They are "
    "clipped before calling.",
)

prim = R.get("primers")
if isinstance(prim, pd.DataFrame) and not prim.empty:
    st.subheader("Recovering the panel design from the data")
    st.markdown(
        "The panel's target list is proprietary, but every Read 2 begins with a primer — so counting distinct "
        "read-2 prefixes reconstructs the primer set, and mapping those primers to the genome reconstructs "
        "the panel."
    )
    prim = prim.copy()
    fig = go.Figure(go.Scatter(x=prim["rank"], y=prim["cum_frac"], mode="lines", line=dict(width=2)))
    fig.update_layout(height=260, margin=dict(l=10, r=10, t=30, b=10),
                      title="Cumulative share of reads explained by the top-N primer candidates",
                      xaxis_title="primer candidates (ranked)", yaxis_title="cumulative fraction of reads")
    st.plotly_chart(fig, width='stretch')
    st.dataframe(prim.head(20), width='stretch', hide_index=True)

merge = R.get("merge_summary") or {}
if merge:
    st.subheader("Caller agreement")
    ba = merge.get("by_agreement", {})
    st.bar_chart(pd.DataFrame({"variants": ba.values()}, index=[f"{k} caller(s)" for k in ba]))
    st.caption("Variants found by two or three independent callers are far more likely to be real. "
               "No single caller is reliable on tumour-only amplicon data, so the pipeline requires "
               "agreement between callers rather than trusting one of them.")
