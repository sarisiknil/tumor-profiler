# Do this next — step-by-step Galaxy run

Everything local is finished. These are the only remaining actions, in order.

## 0 · Before you start (2 min)

```bash
cd ~/projects/workspace/tumor-profiler
python3 scripts/check_tools.py          # should print "All required tools are present"
ls resources/panel_dna.bed              # 1,366 intervals, ~0.49 Mb — the calling target
```
Free some disk if you can: the archive `workspace/archive.zip` (2.6 GB) is byte-for-byte
the same four FASTQ files that are already unpacked. Move it to an external drive **after** the upload succeeds.

## 1 · Create the account (5 min)

<https://usegalaxy.eu> → Register with your Sabancı address → confirm by e-mail.
Then `User → Preferences → Manage API Key → Create a new key` and copy it.

## 2 · Upload (30–60 min, unattended)

`Upload Data` → `Choose local files`, add all four:

```
workspace/TUMOR01/DNA/TUMOR01_DNA_R1_001.fastq.gz
workspace/TUMOR01/DNA/TUMOR01_DNA_R2_001.fastq.gz
workspace/TUMOR01/RNA/TUMOR01_RNA_R1_001.fastq.gz
workspace/TUMOR01/RNA/TUMOR01_RNA_R2_001.fastq.gz
```

Set **Type = `fastqsanger.gz`** and **Genome = `Human Dec. 2013 (GRCh38/hg38)`** for all four before pressing
Start. Also upload `resources/panel_dna.bed` and `resources/panel_rna.bed` (type `bed`).

Name the history `TUMOR01 DNA panel`. While it uploads, continue to step 3.

## 3 · Fetch the germline resource inside Galaxy (no local download)

`Upload Data → Paste/Fetch data`, paste this URL, type `vcf_bgzip`:

```
https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
```

## 4 · DNA workflow (~1–2 h including queue)

Run these tools in order; the exact parameters are in `galaxy/README.md`.

1. **fastp** — UMI: on, location `read1`, length `12`, skip `14`
2. **Map with BWA-MEM** — reference `hg38`, paired, read group `SM=TUMOR01`
3. **UMI-tools dedup** — paired, method `directional`
4. **samtools ampliconclip** — the primer BED, hard clip
5. **mosdepth** — `--by resources/panel_dna.bed`
6. **GATK4 Mutect2** — tumour-only, germline resource = the gnomAD VCF, intervals = `panel_dna.bed`, emit `f1r2`
7. **VarDict** — `-f 0.01`, intervals = `panel_dna.bed`
8. **LoFreq call** — call indels

## 5 · RNA workflow (~1 h)

1. **fastp** — UMI: on, `read1`, length `12`, skip `14`
2. **RNA STAR** — `hg38` + GTF, with the chimeric parameters listed in `galaxy/README.md`
3. **UMI-tools dedup**
4. **Arriba Get Filters**, then **Arriba** — keep `fusions.tsv` **and** `fusions.discarded.tsv`
5. Keep STAR's `SJ.out.tab` — MET exon 14 skipping is in this assay's scope and Arriba cannot see it

## 6 · Bring it back and finish locally (10 min)

```bash
export GALAXY_API_KEY=<your key>
python3 scripts/import_galaxy.py --list                       # find the history name
python3 scripts/import_galaxy.py --history "TUMOR01" --out results/galaxy_import
gatk LearnReadOrientationModel -I results/galaxy_import/f1r2.tar.gz -O results/galaxy_import/orientation.tar.gz
gatk FilterMutectCalls -R <hg38.fa> -V results/galaxy_import/mutect2.vcf \
     --ob-priors results/galaxy_import/orientation.tar.gz -O results/galaxy_import/mutect2.filtered.vcf
snakemake -c4 --config mode=patient
streamlit run dashboard/app.py -- --results results
```

(If you would rather skip `FilterMutectCalls`, the pipeline still works: the merge step keeps Mutect2's own
FILTER column and the caller-concordance rule does the heavy lifting.)

## 7 · What to check first in the output

- `results/validation/concordance.tsv` — **does the pipeline find IDH1 p.R132C at ~28 % VAF?**
  That single line is the project's headline result.
- `results/dna/variants_tiered.tsv` — anything else worth reporting, and at which tier.
- `results/rna/fusions_summary.tsv` and `results/rna/exon_skipping.tsv` — the laboratory reported no fusion;
  agreeing with a negative is also a result.
- `results/report.md` and the dashboard — the figures for the internship report.

## 8 · Export the workflows

`Workflow → Extract workflow from history → Download` → save as `galaxy/dna_panel.ga` and
`galaxy/rna_fusion.ga`. Those two files are the reproducibility deliverable.
