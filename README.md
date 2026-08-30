# Educational DNA–RNA Tumour Profiling Pipeline

A reproducible, documented pipeline that turns raw tumour sequencing reads into an interpretable answer to one
question:

> **What can we understand about a person's disease — and their possible treatments — from their tumour's DNA and RNA?**

It is built around a real clinical assay (an **anchored multiplex PCR** targeted DNA panel plus a targeted RNA
fusion panel, tumour-only, no matched normal) but ships with a **synthetic dataset on a real 230 kB GRCh38
mini-reference**, so the entire workflow — alignment, three somatic callers, live annotation, tiering, pathway
mapping, report and dashboard — runs on a laptop in under a minute, with no patient data involved.

```bash
git clone <this repo> && cd tumor-profiler
bash environment/create_envs.sh        # conda envs from Bioconda
bash resources/fetch_resources.sh      # cancer-gene list (not redistributed here)
python3 scripts/check_tools.py         # verifies every binary is present and runnable
snakemake -c4 --use-conda              # runs the example end to end -> results_example/
streamlit run dashboard/app.py -- --results results_example
pytest tests/ -q                       # 21 tests, including truth-set recovery
```

## What it produces

- `results_*/report.md` — a readable summary: library QC, variants by tier, biomarkers, fusions, pathways, trials
- `results_*/summary.json` — the same content, machine-readable, with a provenance block
- a **six-page Streamlit dashboard**: QC & Library · Variants · RNA & Fusions · Pathways · Therapies · Method & Limits
- per-stage tables (`variants_tiered.tsv`, `fusions_summary.tsv`, `pathway_hits.tsv`, …)

## How it works

| Stage | Tool | Runs on |
|---|---|---|
| Read-structure and primer inference | project scripts, pure Python | laptop |
| UMI extraction and deduplication | `umi_tools` (regex extract, directional dedup) | laptop / Galaxy |
| Alignment to GRCh38 | `bwa mem` | **usegalaxy.eu** (patient) / laptop (example) |
| Primer clipping, coverage | `samtools ampliconclip`, `mosdepth` | laptop / Galaxy |
| Somatic calling | Mutect2 tumour-only · VarDict · LoFreq · FreeBayes | Galaxy / laptop |
| Tumour-only filtering | project script (gnomAD, ClinVar, VAF, FFPE) | laptop |
| Annotation | Ensembl VEP REST | laptop |
| Evidence and tiering | CIViC · OncoKB · DGIdb · ClinicalTrials.gov → AMP/ASCO/CAP + ESCAT | laptop |
| Biomarkers | TMB with confidence interval, mutation spectrum, FFPE indicator | laptop |
| RNA fusions and exon skipping | `STAR` + `Arriba`, splice-junction counting | **usegalaxy.eu** |
| Pathways | Sanchez-Vega ten-pathway membership map | laptop |
| Report and dashboard | project scripts + Streamlit | laptop |

The split exists for a concrete reason: a human `bwa` index needs ~6 GB of RAM and a human STAR index ~31 GB, so
those steps run on Galaxy's free public service (250 GB per account, GRCh38 indexes pre-built), while everything
that fits in a laptop's memory stays local and fast. See [docs/11_reproducibility.md](docs/11_reproducibility.md).

## Two execution modes

```bash
snakemake -c4                          # mode=example  (default): synthetic data, fully local
snakemake -c4 --config mode=patient    # patient data: starts from files imported from Galaxy
```

## Documentation

Each pipeline stage has a page explaining *what it measures*, *what the result means* and *how it can mislead*:

[Overview](docs/00_overview.md) · [Read structure](docs/01_read_structure.md) · [UMIs](docs/02_umi.md) ·
[Primers](docs/03_primers.md) · [Alignment](docs/04_alignment.md) · [Variant calling](docs/05_variant_calling.md) ·
[Tumour-only filtering](docs/06_tumour_only.md) · [Annotation and tiering](docs/07_annotation_tiering.md) ·
[Biomarkers](docs/08_biomarkers.md) · [RNA and fusions](docs/09_rna_fusions.md) · [Pathways](docs/10_pathways.md) ·
[Reproducibility](docs/11_reproducibility.md) · [Ethics and data protection](docs/12_ethics.md) ·
[Data licences](docs/licensing.md)

## The example dataset

`scripts/fetch_example_reference.py` downloads ~230 kB of genuine GRCh38 sequence around six driver hotspots;
`scripts/simulate_amp.py` builds an anchored-multiplex library over it with PCR duplicates, sequencing errors and
six spiked mutations at known allele fractions:

| Gene | Variant | Spiked VAF | What it tests |
|---|---|---|---|
| KRAS | G12C | 0.32 | tier I actionability (sotorasib, adagrasib) |
| BRAF | V600E | 0.18 | tier I in melanoma, tier II elsewhere |
| TP53 | R175H | 0.55 | prognostic evidence, no approved drug |
| EGFR | L858R | 0.07 | low-allele-fraction sensitivity |
| PIK3CA | H1047R | 0.24 | pathway mapping (PI3K) |
| PTEN | R130* | 0.40 | truncating variant; germline-ambiguity flag |

Because the coordinates are real, VEP, ClinVar, gnomAD and CIViC all work on the example exactly as they do on a
patient sample. `tests/` asserts that each one is recovered, classified, annotated and tiered correctly.

## Data protection

No patient data is in this repository and none can be: raw reads, alignments, germline variant lists and the
laboratory's identifiers are excluded by `.gitignore`, the clinical report is kept outside the repository
entirely, and a test fails the build if a laboratory identifier ever appears in a shareable output. Pseudonymised
genomic data remain personal data under GDPR and KVKK — see [docs/12_ethics.md](docs/12_ethics.md).

**This is an educational project. It is not accredited, not clinically validated, and not a medical device.
Nothing it produces should inform a treatment decision.**

## Licence

MIT for the code. Third-party data sources keep their own licences — see [docs/licensing.md](docs/licensing.md).
