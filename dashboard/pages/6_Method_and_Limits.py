"""Page 6 — how the numbers were produced, and what they cannot support."""
import json
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import load_results, header  # noqa: E402

st.set_page_config(page_title="Method & Limits", page_icon="📋", layout="wide")
R = load_results()
header("6 · Method, provenance and limits", R)
s = R.get("summary", {})

st.subheader("Where each step ran")
st.markdown(
    """
| Stage | Tool | Where it runs | Why |
|---|---|---|---|
| Read structure, primer inference | project scripts (pure Python) | laptop | seconds; no reference needed |
| UMI extraction | `umi_tools extract` (regex) | laptop | the fixed common region must be discarded, not re-attached |
| Alignment to GRCh38 | `bwa mem` | **usegalaxy.eu** | the human genome index needs ~6 GB RAM |
| UMI deduplication | `umi_tools dedup --method directional` | Galaxy / laptop | positional dedup would collapse real molecules |
| Primer clipping | `samtools ampliconclip` | Galaxy / laptop | primer bases are synthetic, not patient sequence |
| Somatic calling | Mutect2 (tumour-only) · VarDict · LoFreq · FreeBayes | Galaxy / laptop | no single caller is reliable; agreement is the filter |
| RNA alignment + fusions | `STAR` + `Arriba` | **usegalaxy.eu** | STAR's human index needs ~31 GB RAM |
| Annotation | Ensembl VEP REST API | laptop | avoids a 20 GB offline cache; a panel is only a few hundred variants |
| Evidence | CIViC · OncoKB · DGIdb · ClinicalTrials.gov | laptop | live APIs, free for academic use |
| Interpretation, report, dashboard | project scripts + Streamlit | laptop | reads only the pipeline's own JSON/TSV |
"""
)

st.subheader("What this analysis cannot tell us")
st.markdown(
    """
- **Whether a variant is inherited.** No normal tissue was sequenced. Population-frequency filtering removes
  common polymorphisms but cannot remove rare, family-private germline variants (Sun et al. 2018).
- **A trustworthy mutational burden or mutational signature.** Both need far more sequenced territory than a
  targeted panel provides; signature fitting needs on the order of 100–200 mutations.
- **Microsatellite instability**, unless the panel covers enough microsatellite loci and a panel-specific model
  exists — it does not here.
- **Copy number, quantitatively.** Without a reference cohort processed identically, apparent gains and losses
  are confounded by primer efficiency.
- **The tumour's tissue of origin or its whole transcriptome.** The RNA assay is a small fusion panel, not
  RNA-seq: expression is measured only for primed targets.
- **Anything about prognosis for this individual.** Evidence linked here describes cohorts, not a person.
"""
)

st.subheader("Data protection")
st.markdown(
    """
The patient's raw reads, alignments and per-sample files never leave the local machine and are excluded from
version control. Only **somatic** findings are shown, under an alias; germline-looking variants are flagged but
not published. Under GDPR — and the Turkish KVKK, which classes genetic data as special-category data —
pseudonymised genomic data remain personal data, so 'we removed the name' is not anonymisation. The reference
practice followed here is the one used by large cancer consortia: share somatic variants, withhold raw reads
and germline genotypes.
"""
)

prov = (s.get("_provenance") or {})
if prov:
    st.subheader("Provenance of this run")
    st.json(prov)

st.subheader("Machine-readable summary")
st.caption("Everything on these pages is derived from a single file, so the report and the dashboard can never "
           "disagree.")
st.download_button("Download summary.json",
                   data=json.dumps(s, indent=1),
                   file_name="summary.json", mime="application/json")
md = Path(R["dir"]) / "report.md"
if md.exists():
    st.download_button("Download report.md", data=md.read_text(), file_name="report.md", mime="text/markdown")
    with st.expander("Preview the generated report"):
        st.markdown(md.read_text())
