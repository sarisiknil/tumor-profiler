# The dashboard

Six pages that read **only** the files a pipeline run produced — no database, no raw reads, no patient
identifiers. Everything shown comes from `results/summary.json` and the tables beside it, so the dashboard and
the written report can never disagree.

## Running it

```bash
conda activate tp-py
streamlit run dashboard/app.py                                  # picks up results/ if present
streamlit run dashboard/app.py -- --results results_example     # the synthetic example instead
```

It opens at <http://localhost:8501>. Use `--results results_example` when demonstrating the project to anyone
who should not see patient-derived output.

## The pages

| Page | What it answers |
|---|---|
| **Home** | What the project asks, and what a tumour-only panel can and cannot establish |
| **1 · QC & Library** | What assay produced these reads? Base composition showing the barcode and adapter, primer-concentration curve, caller agreement |
| **2 · Variants** | Agreement with the accredited laboratory, artefact screening, tiered variants, every rejected call with its reason, TMB and mutation spectrum |
| **3 · RNA & Fusions** | Candidate fusions with artefact flags, exon-skipping screen, why the open caller differs from the vendor's |
| **4 · Pathways** | Which oncogenic pathways carry an alteration — membership, explicitly not enrichment |
| **5 · Therapies** | CIViC/OncoKB evidence, drug–gene interactions, recruiting trials, and the licensing behind each source |
| **6 · Method & Limits** | Where each step ran, what the analysis cannot support, data protection, and the provenance block |

## Design rule

Every page states what its numbers mean *and* how they mislead. The two findings that matter most in this
sample are surfaced rather than buried: the exact allele-fraction agreement with the laboratory on page 2, and
the artefact screen that turns twenty candidate fusions into a negative result on page 3.

## Tests

`tests/test_dashboard.py` renders every page headlessly with Streamlit's `AppTest` and fails if any raises, so
a change to the results schema breaks the build rather than the demo.
