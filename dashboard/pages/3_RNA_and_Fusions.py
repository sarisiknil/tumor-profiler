"""Page 3 — what the RNA adds that DNA cannot see."""
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import load_results, header, teaching_note, caveats  # noqa: E402

st.set_page_config(page_title="RNA & Fusions", page_icon="🔗", layout="wide")
R = load_results()
header("3 · RNA: fusions and splicing", R)

st.markdown(
    "A DNA panel reads the genome; an RNA panel reads what the tumour is actually transcribing. That difference "
    "matters most for **gene fusions**: the DNA breakpoint usually lies deep inside a large intron that no panel "
    "covers, while the fused transcript is short, abundant and unambiguous. In lung adenocarcinomas where DNA "
    "sequencing found no driver, RNA sequencing recovered a targetable kinase fusion in a substantial fraction "
    "of cases (Benayed et al. 2019) — which is exactly why this assay pairs a DNA panel with an RNA panel."
)

fus = R.get("fusions")
fsum = R.get("fusion_summary") or {}
if isinstance(fus, pd.DataFrame) and not fus.empty:
    st.subheader("Candidate fusions, after artefact screening")
    flagged = (fus["artefact_risk"] != "not flagged").sum() if "artefact_risk" in fus.columns else 0
    if flagged == len(fus) and len(fus):
        st.error(
            f"**All {len(fus)} candidate fusions carry an artefact flag.** They are the paralogous FGFR "
            "family joined in every possible combination, in both orientations, on a handful of reads each — "
            "the signature of cross-priming between genes the assay amplifies in the same tube. The "
            "accredited laboratory reported no fusion, and after screening neither does this pipeline. "
            "An unscreened reading would have reported an FGFR2 fusion, the single most consequential false "
            "positive available in this tumour type.", icon="🚫")
    elif flagged:
        st.warning(f"{flagged} of {len(fus)} candidates carry an artefact flag.", icon="⚠️")
    cols = [c for c in ["fusion", "artefact_risk", "confidence", "reading_frame", "targetable",
                        "support_total", "artefact_reasons", "breakpoint1", "breakpoint2", "notes"]
            if c in fus.columns]
    st.dataframe(
        fus[cols],
        width='stretch', hide_index=True,
        column_config={"support_total": st.column_config.NumberColumn("supporting reads"),
                       "notes": st.column_config.TextColumn("interpretation", width="large")},
    )
    for _, f in fus.iterrows():
        if f.get("targetable") and f.get("artefact_risk", "not flagged") == "not flagged":
            st.success(f"**{f['fusion']}** involves {f['targetable']}, a kinase with approved inhibitors. "
                       f"Confidence: {f['confidence']}; {f['support_total']} supporting reads. "
                       "A finding like this would need orthogonal confirmation (FISH, IHC or RT-PCR) before "
                       "it could influence treatment.", icon="🎯")
else:
    st.info("No fusion table available for this run. In patient mode the RNA reads are aligned with STAR and "
            "screened with Arriba on usegalaxy.eu; the resulting `fusions.tsv` is then imported here.")

st.divider()
sk = R.get("exon_skipping")
if isinstance(sk, pd.DataFrame) and not sk.empty:
    st.subheader("Exon-skipping events")
    st.dataframe(sk, width='stretch', hide_index=True)
    st.caption("Fusion callers do not look for exon skipping — it is a splicing change inside a single gene, "
               "not a join between two genes. MET exon 14 skipping is the clinically important example: it is "
               "targetable, and RNA detects it far more reliably than DNA, because the causal DNA variants are "
               "scattered across the flanking introns and splice sites.")

teaching_note(
    "Why the free fusion caller is not the vendor's caller",
    "This library was designed for the manufacturer's analysis software, which knows the primer set. Open-source "
    "callers do not, and their filters were tuned for whole-transcriptome data. Benchmarked on real Archer "
    "FusionPlex samples, Arriba recovered 86 % (lung panel) and 57 % (sarcoma panel) of the vendor's fusion "
    "calls, while STAR-Fusion recovered 33 % and 7 % (Capone et al. 2022). The pipeline therefore uses Arriba, "
    "keeps low-confidence calls instead of discarding them, and states this sensitivity gap in the report rather "
    "than pretending the open pipeline is equivalent.",
)
caveats((fsum or {}).get("caveats"))
