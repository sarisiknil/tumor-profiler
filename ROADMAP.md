# Project roadmap — what is done and what is left

Internship project XX395, Summer 2026. Report deadline: mid/late September 2026.

## Phase 1 · Understand the data — **done**

- [x] Identified the assay from the raw reads alone (12-nt UMI + fixed common region, one-sided primers,
      8+10 nt indices) and confirmed it against the laboratory report: **Archer VariantPlex Expanded Solid Tumor**
      (DNA) and **Archer FusionPlex Lung v2** (RNA), AMP chemistry, NextSeq 550 + NovaSeq 6000
- [x] Quantified the libraries: 5.54 M DNA read pairs, 13.87 M RNA read pairs, 29–33 % adapter read-through on
      the DNA arm (consistent with FFPE), RNA already adapter-trimmed
- [x] Reconstructed the panel design from Read 2: 4,201 DNA and 477 RNA primer candidates
- [x] Built the real target BEDs from the assay's gene lists: 74 DNA genes / 1,366 intervals / 0.49 Mb,
      16 RNA genes / 338 intervals / 0.13 Mb
- [x] Extracted the laboratory's findings: **IDH1 p.R132C, 28.1 % VAF, 501×, Tier I-A**; no fusions;
      MSI not calculable; diagnosis cholangiocarcinoma, 70 % tumour content
- [x] Converted the reported hg19 coordinate to GRCh38 (chr2:208248389 G>A), verified both directions

## Phase 2 · Build the pipeline — **done**

- [x] Snakemake workflow, two modes (`example` fully local, `patient` via Galaxy), Bioconda environments
- [x] UMI extraction and UMI-aware deduplication; primer clipping; coverage
- [x] Three somatic callers with concordance filtering; tumour-only classification with an explicit
      germline-ambiguity flag; FFPE artefact indicator
- [x] Annotation via Ensembl VEP REST; evidence from CIViC, OncoKB, DGIdb, ClinicalTrials.gov
- [x] AMP/ASCO/CAP + ESCAT tiering; TMB with a confidence interval; copy-number hints; pathway mapping
- [x] RNA fusion summary and exon-skipping detection (MET exon 14 is in this assay's scope)
- [x] Report generator, `summary.json`, six-page Streamlit dashboard
- [x] Validation module comparing our calls with the laboratory's report
- [x] Synthetic example dataset on a real 230 kB GRCh38 mini-reference with six spiked hotspots
- [x] 21 tests, GitHub Actions CI, clean-clone verification
- [x] 14 documentation pages, README, licensing page, 29 papers downloaded

## Phase 3 · Run the patient sample — **next, ~1 day**

- [ ] **Step 1** Galaxy account + API key
- [ ] **Step 2** upload (`scripts/galaxy_upload.py`, ~1 h unattended)
- [ ] **Step 3** DNA arm on Galaxy: UMI extract → BWA-MEM → dedup → ampliconclip → mosdepth → Mutect2 + VarDict + LoFreq
- [ ] **Step 4** RNA arm on Galaxy: UMI extract → STAR → dedup → Arriba (+ discarded fusions, + SJ.out.tab)
- [ ] **Step 5** `scripts/import_galaxy.py`
- [ ] **Step 6** `snakemake -c4 --config mode=patient`, then the dashboard
- [ ] **Step 7** export both Galaxy workflows as `.ga` and commit them

Everything in Phase 3 is written out click-by-click in `galaxy/UPLOAD_CHECKLIST.md`.

## Phase 4 · Interpret the result — ~1 day after Phase 3

- [ ] Concordance table: did we recover IDH1 R132C, and at what VAF? Explain any difference (hg19 vs GRCh38,
      vendor consensus reads vs our deduplication, 70 % tumour content)
- [ ] Explain every additional call we make that the laboratory did not report — most will be below its 5 %
      reporting threshold or Tier III/IV, which the laboratory does not report at all
- [ ] Confirm the negative findings: no fusion, no MET exon 14 skipping, MSI not assessable
- [ ] Per-gene coverage table: which of the 74 genes were actually assessable at ≥100×
- [ ] Runtime and cost table: Galaxy CPU-hours versus what a commercial re-analysis would cost

## Phase 5 · Publish the repository — half a day

- [ ] `git remote add origin` and push to GitHub (public)
- [ ] Confirm CI passes on GitHub
- [x] Final audit: history rewritten with `git filter-repo` to purge the laboratory sample code,
      the patient's first name (which was the workspace folder name), and local paths. Verified
      across all 263 blobs in every commit. Re-run with:
      `git rev-list --objects --all | ... | git cat-file blob | grep -i <pattern>`
- [ ] Add the exported `.ga` workflows and a short screencast or screenshots of the dashboard

## Phase 6 · Write the deliverables — ~4 days, start by 13 September at the latest

Three files, named as the faculty guidelines require
(`BIO395_FinalReport_FIRSTNAME_LASTNAME_DDMonthYYYY.doc` etc.):

- [ ] **Report** — section mapping:
  - 3.3 Motivation: vendor software is closed and per-sample priced; the analysis behind a clinical report
    should be inspectable and re-runnable
  - 3.4 Related literature (≤3 pages): AMP chemistry and UMIs · FFPE artefacts · tumour-only calling ·
    AMP/ASCO/CAP and ESCAT tiering · DNA+RNA complementarity for fusions · reproducibility
    (Galaxy, Snakemake, Bioconda, FAIR). All 29 PDFs are in `workspace/literature/` with a README mapping each
    paper to where it is used
  - 4.3 Methodology: the stage table from `docs/00_overview.md`
  - 4.5 Details (≤10 pages): one subsection and one figure per stage
  - 4.6 Results (1 page): reads → molecules, coverage, variants per tier, concordance with the laboratory,
    runtimes
  - 5.3 Difficulties: 8 GB RAM versus a 31 GB STAR index (solved with Galaxy); tumour-only germline ambiguity;
    open-source fusion callers being less sensitive than the vendor's on amplicon data
  - Appendix: full variant table (anonymised), Galaxy workflow screenshots, JSON schema
- [ ] **Presentation** (10–15 min, Times New Roman) — motivation, objective, methods, deliverables, results
- [ ] **Digest** — one slide: title, company, dates, objectives, six outcomes

## Points worth making in the discussion

1. The kit was identified from the raw reads before the report was read — the read structure alone was enough.
2. **A lung fusion panel was applied to a bile-duct cancer.** It does cover FGFR2, the main actionable fusion in
   cholangiocarcinoma, but that is worth stating rather than assuming.
3. **The DNA panel omits several of the genes most often mutated in cholangiocarcinoma** — ARID1A, BAP1, PBRM1,
   ELF3, KMT2D. A negative result there is an untested region, not an absence of mutation. This is the clearest
   possible illustration of the project's central point: *what you can learn is bounded by what the assay looked at.*
4. The laboratory's MSI call was "not calculable" — matching this pipeline's a-priori refusal to report MSI from
   a panel without a matched baseline. Agreeing about what *cannot* be measured is also a result.
5. Tumour content is 70 %, so a clonal heterozygous somatic variant should sit near 35 % VAF. The reported
   28.1 % is consistent with a clonal IDH1 mutation — a small piece of quantitative reasoning worth showing.
