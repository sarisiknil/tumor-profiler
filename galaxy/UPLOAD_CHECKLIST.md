# Galaxy run — do this one step at a time

Tool versions below were read from the usegalaxy.eu API on 2026-08-30. Search the tool name in Galaxy's tool
panel (top left) and check the version shown under the tool title matches.

---

## Step 1 · Account and API key (5 min)

1. Go to <https://usegalaxy.eu> → **Login or Register** → register with your Sabancı address.
2. Confirm the e-mail Galaxy sends you.
3. Go to <https://usegalaxy.eu/user/api_key> → **Create a new key** → copy it.
4. In a terminal:

```bash
cd ~/projects/workspace/tumor-profiler
export GALAXY_API_KEY=<paste the key>
conda activate tp-py
python3 scripts/galaxy_upload.py --check
```

You should see your username and `disk usage 0 bytes of quota 250.0 GB`. Keep this terminal open — the key
lives only in it. (If you close it, `export` it again.)

---

## Step 2 · Upload everything with one command (~40–60 min, unattended)

```bash
python3 scripts/galaxy_upload.py --dry-run          # shows the 6 files and their total size
python3 scripts/galaxy_upload.py --history "TUMOR01 tumour profiling"
```

This creates the history, uploads the four FASTQ files as `fastqsanger.gz` with genome `hg38`, uploads
`resources/panel_dna.bed` and `resources/panel_rna.bed`, and asks Galaxy to download the gnomAD germline
resource directly from Google's server (it never touches your laptop). It prints a link to watch progress and
waits until every dataset is green.

**If the script fails for any reason**, do it in the browser instead: `Upload Data → Choose local files`, add the
four FASTQ files, set **Type `fastqsanger.gz`** and **Genome `Human Dec. 2013 (GRCh38/hg38)`** for each, then
`Start`. Add the two BED files as type `bed`. Then `Paste/Fetch data` with this URL, type `vcf_bgzip`:

```
https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
```

While it uploads, read `docs/02_umi.md` and `docs/05_variant_calling.md` — those two pages are what Step 4 does.

---

## Step 3 · DNA arm

### The short way

```bash
export GALAXY_API_KEY=<your key>
python3 scripts/galaxy_run.py --history "TUMOR01 tumour profiling" --arm dna --dry-run   # validates only
python3 scripts/galaxy_run.py --history "TUMOR01 tumour profiling" --arm dna
```

`galaxy_run.py` submits the six tools in order, waits for each, renames the outputs to the names
`import_galaxy.py` expects, and skips any step whose output already exists — so an interrupted run just
restarts. Every parameter path is checked against the tool definition Galaxy currently serves *before*
anything is submitted, so a renamed parameter fails loudly instead of silently falling back to a default.
`tests/test_galaxy_params.py` runs the same check in CI.

### The long way — the same thing in the web interface

Run each tool, wait for it to turn green, then start the next. Total ~1–2 h including queue time.

### 3.1 UMI-tools extract — `UMI-tools extract` v1.1.6+galaxy0

| Field | Value |
|---|---|
| Library type | **Paired-end** |
| Forward reads | `TUMOR01_DNA_R1_001.fastq.gz` |
| Reverse reads | `TUMOR01_DNA_R2_001.fastq.gz` |
| Barcode pattern for first read | `^(?P<umi_1>.{12})(?P<discard_1>AGTCGTCTCGAAGT?){s<=2}` |
| Barcode pattern for second read | *(leave empty)* |
| Extract method | **regex** |
| Output log | Yes |

This moves the 12-nt molecular barcode into the read name and deletes the fixed adapter region. Do **not** use
fastp here: Galaxy's fastp wrapper exposes `umi_len` but not `umi_skip`, so it cannot remove the common region.

### 3.2 Map with BWA-MEM — `Map with BWA-MEM` v0.7.19+galaxy1

| Field | Value |
|---|---|
| Reference genome source | Use a built-in genome index |
| Using reference genome | **Human: hg38** |
| Single or Paired-end reads | Paired |
| Select first / second set | the two outputs of step 3.1 |
| Set read groups information | **Set read groups (SAM/BAM specification)** → ID `TUMOR01`, SM `TUMOR01`, PL `ILLUMINA`, LB `VariantPlex` |
| Select analysis mode | 1. Simple Illumina mode |

### 3.3 UMI-tools deduplicate — `UMI-tools deduplicate` v1.1.6+galaxy0

| Field | Value |
|---|---|
| Reads to deduplicate | the BAM from 3.2 |
| Library type | **Paired-end** |
| Method to identify group of reads | **directional** |
| Output log | Yes |

Check the log afterwards: reads-in versus unique-molecules-out is your library-complexity number, and it belongs
in the report.

### 3.4 Samtools ampliconclip — `Samtools ampliconclip` v1.22+galaxy2

| Field | Value |
|---|---|
| BAM file | the deduplicated BAM |
| BED file of amplicon primers | `panel_dna.bed` |
| Clipping | **Hard clip** |
| Clip both ends | Yes |

### 3.5 mosdepth — `mosdepth` v0.3.8+galaxy0

| Field | Value |
|---|---|
| BAM/CRAM | the clipped BAM |
| Compute coverage over | **regions in a BED file** → `panel_dna.bed` |
| Per-base output | **No** |

### 3.6 GATK4 Mutect2 — `GATK4 Mutect2` v4.6.2.0+galaxy0

| Field | Value |
|---|---|
| Analysis mode | **tumor_only** |
| Input Tumor BAM | the clipped BAM |
| Reference source | built-in → **hg38** |
| Germline Resource | `af-only-gnomad.hg38.vcf.gz` |
| Intervals File | `panel_dna.bed` |
| GZIP Output | No |

### 3.7 VarDict — `VarDict` v1.8.4+galaxy0

| Field | Value |
|---|---|
| Reference | built-in **hg38** |
| BAM | the clipped BAM |
| Regions | `panel_dna.bed` |
| Sample name | `TUMOR01` |
| Allele frequency threshold | **0.01** |
| Output type | VCF |

### 3.8 LoFreq Call variants — `Call variants` v2.1.5+galaxy3

| Field | Value |
|---|---|
| Reference | built-in **hg38** |
| BAM | the clipped BAM |
| Call indels | Yes |
| Variant calling parameters | defaults are fine |

---

## Step 4 · RNA arm (~1 h)

```bash
python3 scripts/galaxy_run.py --history "TUMOR01 tumour profiling" --arm rna
```

It fetches the GENCODE annotation into the history itself, then runs STAR with the wrapper's built-in
**`arriba` preset** (which sets every chimeric parameter for you), Arriba Get Filters for hg38, and Arriba
with the discarded-fusion output kept. The manual equivalent is below.

### 4.1 UMI-tools extract — same tool as 3.1

Same settings, but the RNA files and this pattern:

```
^(?P<umi_1>.{12})(?P<discard_1>CTGGATAGTACGCT){s<=2}
```

### 4.2 Get a gene annotation

`Upload Data → Paste/Fetch data`, type `gtf`, genome `hg38`:

```
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
```

### 4.3 RNA STAR — `RNA STAR` v2.7.8a+galaxy1

| Field | Value |
|---|---|
| Single-end or paired-end | Paired-end (as individual datasets) → the two outputs of 4.1 |
| Reference genome | **Use a built-in index** → *genome with no built-in gene-model* → **hg38**, then supply the GENCODE GTF as the gene model, splice-junction overhang **150** |
| Per gene read counts | **GeneCounts** |
| Chimeric alignments | output type **WithinBAM SoftClip** |
| Computational settings | **Use parameters suggested for Arriba** — this one preset sets every chimeric parameter correctly |
| Output splice junctions (SJ.out.tab) | Yes — **you need this file** |

### 4.4 UMI-tools deduplicate — on the STAR BAM, paired, directional.

### 4.5 Arriba Get Filters — `Arriba Get Filters` v2.5.1+galaxy1 → genome **hg38**.

### 4.6 Arriba — `Arriba` v2.5.1+galaxy1

| Field | Value |
|---|---|
| STAR BAM | the deduplicated STAR BAM |
| Genome | built-in **hg38** FASTA |
| Gene annotation | the GENCODE GTF |
| Blacklist / known fusions / protein domains | the three files from 4.5 |
| Output discarded fusions | **Yes** — Arriba's filters are tuned for whole-transcriptome data and will over-filter amplicon panels |

---

## Step 5 · Bring the results back (10 min)

```bash
export GALAXY_API_KEY=<your key>
python3 scripts/import_galaxy.py --list
python3 scripts/import_galaxy.py --history "TUMOR01" --out results/galaxy_import
ls results/galaxy_import        # mutect2.vcf, vardict.vcf, lofreq.vcf, coverage.regions.bed.gz,
                                # fusions.tsv, fusions.discarded.tsv, SJ.out.tab
```

Only small text files come down; the BAMs stay on Galaxy. If a file is missing, rename the Galaxy dataset so its
name contains `mutect`, `vardict`, `lofreq`, `regions`, `fusions` or `SJ` and re-run the import.

**`FilterMutectCalls` is skipped on purpose.** It needs a local 3 GB copy of hg38 and your disk is nearly full;
the pipeline instead keeps Mutect2's own FILTER column and requires two of three callers to agree. Say so in the
methods section — it is a defensible, documented choice, not an omission.

---

## Step 6 · Finish locally (5 min)

```bash
snakemake -c4 --config mode=patient
streamlit run dashboard/app.py -- --results results
```

Then open, in this order:

1. `results/validation/concordance.tsv` — **the headline: is IDH1 p.R132C there at ~28 % VAF?**
2. `results/dna/variants_tiered.tsv` — everything else, and at which tier.
3. `results/rna/fusions_summary.tsv`, `results/rna/exon_skipping.tsv` — the laboratory reported no fusion.
4. `results/report.md` and the dashboard — your report figures.

---

## Step 7 · Export the workflows (5 min)

In Galaxy: `Workflow → Extract workflow from history`, once for the DNA history and once for the RNA one.
Then `Workflow → Download` and save them here as `galaxy/dna_panel.ga` and `galaxy/rna_fusion.ga`, and commit.
Those two files let anyone re-run the heavy half exactly as you did — they are the reproducibility deliverable.
