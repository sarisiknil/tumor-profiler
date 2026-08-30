"""Page 2 — the DNA findings: which mutations, how confident, and what they mean."""
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import load_results, header, teaching_note, caveats, CLASS_HELP, TIER_COLOR  # noqa: E402

st.set_page_config(page_title="Variants", page_icon="🧪", layout="wide")
R = load_results()
header("2 · DNA variants", R)

tier = R.get("variants_tiered")
filt = R.get("variants_filtered")
if not isinstance(filt, pd.DataFrame) or filt.empty:
    st.warning("No variant table found — run the pipeline first.")
    st.stop()

st.subheader("From raw calls to interpretable variants")
counts = (R.get("filter_summary") or {}).get("counts", {})
cols = st.columns(max(1, len(counts)))
for col, (k, v) in zip(cols, counts.items()):
    col.metric(k.replace("_", " ").title(), v, help=CLASS_HELP.get(k, ""))

amb = (R.get("filter_summary") or {}).get("germline_ambiguous", 0)
if amb:
    st.warning(
        f"**{amb} variant(s) cannot be assigned to tumour or germline.** They sit at a heterozygous allele "
        "fraction in a cancer-predisposition gene and are rare in the population. Distinguishing them requires "
        "sequencing normal tissue from the same patient. If one of these turned out to be inherited it would "
        "matter for the patient's relatives, which is why such findings are governed by consent and genetic "
        "counselling rather than by a pipeline.",
        icon="⚠️",
    )

if isinstance(tier, pd.DataFrame) and not tier.empty:
    st.subheader("Clinically significant variants")
    show = tier.copy()
    show["variant"] = show["gene"].astype(str) + " " + show["hgvsp"].astype(str).str.split(":").str[-1]
    fig = px.scatter(
        show, x="vaf", y="variant", color="amp_tier", size="depth",
        color_discrete_map=TIER_COLOR, size_max=22,
        labels={"vaf": "variant allele fraction", "variant": "", "amp_tier": "AMP/ASCO/CAP tier"},
        hover_data=["consequence", "escat", "n_callers", "gnomad_af_max"],
    )
    fig.update_layout(height=90 + 42 * len(show), margin=dict(l=10, r=10, t=30, b=10),
                      xaxis_range=[0, max(0.6, float(show["vaf"].max() or 0.6) * 1.15)])
    st.plotly_chart(fig, use_container_width=True)
    st.caption("The horizontal position is the fraction of DNA molecules carrying the mutation. Variants near "
               "the highest fraction were probably present in the founding tumour clone; low-fraction variants "
               "are either subclonal (present in part of the tumour) or diluted by normal cells in the biopsy.")

    st.dataframe(
        show[["gene", "hgvsp", "consequence", "amp_tier", "escat", "vaf", "depth", "n_callers",
              "therapies", "tier_rationale"]],
        use_container_width=True, hide_index=True,
        column_config={"hgvsp": st.column_config.TextColumn("protein change"),
                       "amp_tier": st.column_config.TextColumn("tier", width="small"),
                       "vaf": st.column_config.NumberColumn("VAF", format="%.3f"),
                       "tier_rationale": st.column_config.TextColumn("why this tier", width="large")},
    )

teaching_note(
    "How a tier is decided",
    "**AMP/ASCO/CAP (Li et al. 2017)** sorts variants by the strength of clinical evidence, not by how "
    "biologically interesting they are:\n\n"
    "- **Tier I** — a drug approved for this alteration in this tumour type, or a professional guideline.\n"
    "- **Tier II** — a drug approved in a different tumour type, or clinical-trial evidence.\n"
    "- **Tier III** — variant of unknown clinical significance (often a real driver with no drug).\n"
    "- **Tier IV** — benign or likely benign.\n\n"
    "**ESCAT (Mateo et al. 2018)** answers a different question — *how ready is this target for routine use* — "
    "which is why a variant can be Tier II but ESCAT IV. Both are shown so the difference is visible.",
)

st.subheader("All calls, including the ones that were rejected")
st.markdown("Nothing is silently deleted: every variant keeps the reason it was classified the way it was. "
            "This is what makes the pipeline auditable.")
klass = st.multiselect("Show classes", sorted(filt["class"].unique()), default=sorted(filt["class"].unique()))
st.dataframe(filt[filt["class"].isin(klass)][
    ["gene", "hgvsp", "consequence", "class", "germline_ambiguous", "vaf", "depth", "alt_reads",
     "n_callers", "callers", "gnomad_af_max", "clinvar_sig", "reasons"]],
    use_container_width=True, hide_index=True)

bm = R.get("biomarkers") or {}
if bm:
    st.divider()
    st.subheader("Biomarkers")
    tmb = bm.get("tmb", {})
    c1, c2, c3 = st.columns(3)
    c1.metric("TMB (mutations/Mb)", tmb.get("tmb_per_mb", "—"),
              help="Nonsynonymous somatic variants divided by the assessable target size.")
    ci = tmb.get("ci95_per_mb") or [None, None]
    c2.metric("95 % confidence interval", f"{ci[0]} – {ci[1]}")
    c3.metric("Assessable territory", f"{tmb.get('assessable_mb', 0)} Mb")
    if tmb.get("meaningful") is False:
        st.error("The assessable territory is far too small for a meaningful TMB. The number is displayed only "
                 "to show how the calculation works.", icon="🚫")
    spec = bm.get("mutation_spectrum") or {}
    if spec:
        fig = go.Figure(go.Bar(x=list(spec.keys()), y=list(spec.values()), marker_color="#3f6fa8"))
        fig.update_layout(height=240, margin=dict(l=10, r=10, t=30, b=10),
                          title="Substitution spectrum of somatic SNVs", yaxis_title="count")
        st.plotly_chart(fig, use_container_width=True)
        ff = bm.get("ffpe_indicator", {})
        st.caption(f"C>T / G>A share: **{ff.get('C>T_or_G>A_fraction_of_snvs')}** "
                   f"({ff.get('low_vaf_C>T_count')} of them below 10 % VAF). {ff.get('note','')}")
    caveats(tmb.get("caveats"), "What the TMB number does and does not mean")
    nc = bm.get("not_computed") or {}
    if nc:
        st.markdown("**Deliberately not computed**")
        for k, v in nc.items():
            st.markdown(f"- **{k}** — {v}")
