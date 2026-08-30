# Running the heavy steps on usegalaxy.eu

The patient sample is aligned against the whole human genome, which needs more memory than a laptop has
(`bwa` ~6 GB, `bwa-mem2` ~10 GB, `STAR` ~31 GB). Galaxy's free public service provides those resources, keeps a
full record of every tool version and parameter, and exports the whole thing as a re-runnable `.ga` workflow —
which is exactly what a reproducibility deliverable needs.

- Server: <https://usegalaxy.eu> (250 GB per registered account)
- Register with an institutional address; no payment, no quota request needed for this project's ~3 GB.

## A · Upload

Upload the four FASTQ files (`Datasets → Upload`, or FTP for large files). Set the datatype to `fastqsanger.gz`
so Galaxy does not try to re-compress them. Also upload, **by URL** (Galaxy fetches server-side, nothing is
downloaded to the laptop):

```
https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
```

## B · DNA workflow

| # | Tool (Galaxy tool id) | Key parameters |
|---|---|---|
| 1 | **fastp** `iuc/fastp` | UMI processing: enabled · location `read1` · length `12` · skip `14`; report HTML+JSON |
| 2 | **Map with BWA-MEM** `devteam/bwa_mem` | reference: built-in `hg38` · paired · read group `SM=TUMOR01` · algorithm `mem` |
| 3 | **UMI-tools dedup** `iuc/umi_tools_dedup` | paired `--paired` · method `directional` · output stats |
| 4 | **samtools ampliconclip** | primer BED from `results/primers/gsp2_DNA.bed` · hard clip · both ends |
| 5 | **mosdepth** | `--by` the panel BED · no per-base |
| 6a | **GATK4 Mutect2** `iuc/gatk4_mutect2` | tumour-only (single input) · germline resource: the gnomAD VCF · emit `--f1r2-tar-gz` |
| 6b | **VarDict** `iuc/vardict_java` | `-f 0.01 -c 1 -S 2 -E 3 -g 4`, panel BED |
| 6c | **LoFreq call** `iuc/lofreq_call` | call indels · no default filter |
| 7 | **Ensembl VEP** *(optional — the local pipeline uses the REST API instead)* | cache `homo_sapiens_vep_106_GRCh38` |

`FilterMutectCalls` and `LearnReadOrientationModel` are **not** available as Galaxy tools; download Mutect2's
VCF, stats and `f1r2` tarball and run them locally:

```bash
gatk LearnReadOrientationModel -I f1r2.tar.gz -O orientation.tar.gz
gatk FilterMutectCalls -R GRCh38.fa -V mutect2.vcf.gz --ob-priors orientation.tar.gz -O filtered.vcf.gz
```

### Naming the panel's genes

Upload `results/primers/gsp2_DNA.fa` and `results/primers/gsp2_RNA.fa` (produced locally by
`scripts/infer_primers.py`) and run **Map with BWA-MEM** against built-in `hg38` on each, then
**BAM-to-BED** (`bedtools bamtobed`). Download the two BEDs as `primers_mapped.bed` and run:

```bash
python3 scripts/map_primers.py --bed results/galaxy_import/primers_mapped.bed -o results/primers/panel
```

This produces the panel's gene list and the number of primers per gene — the answer to "which genes could this
assay have seen at all?", which every negative result depends on.

## C · RNA workflow

| # | Tool | Key parameters |
|---|---|---|
| 1 | **fastp** | UMI: `read1`, length `12`, skip `14` (the RNA common region is `CTGGATAGTACGCT`) |
| 2 | **RNA STAR** `iuc/rgrnastar` | built-in `hg38` index + GTF; chimeric settings below |
| 3 | **UMI-tools dedup** | paired, directional |
| 4 | **Arriba** `iuc/arriba` (2.5.1) | genome FASTA + GTF + the STAR BAM; blacklist / known-fusions / protein-domains from **Arriba Get Filters** |
| 5 | **Arriba Draw Fusions** | produces the publication-quality fusion diagrams |
| 6 | **featureCounts** | per-target expression (relative only — this is a fusion panel, not RNA-seq) |

STAR chimeric parameters Arriba expects:

```
--chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10
--chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1
--alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 --peOverlapNbasesMin 10
--outSAMunmapped Within
```

Expect Arriba's coverage-based filters to be conservative on amplicon data; keep `fusions.discarded.tsv` too —
`scripts/fusion_summary.py` scans it for discarded calls that involve a targetable kinase.

## D · Bring the results back

```bash
python3 scripts/import_galaxy.py --history "Tumour profiling — DNA" --out results/galaxy_import
snakemake -c4 --config mode=patient
```

`import_galaxy.py` uses the Galaxy API (`GALAXY_API_KEY` from `User → Preferences → Manage API Key`) and
downloads only the small files the local half needs: the VCFs, `mosdepth` regions, Arriba's `fusions.tsv` and
`fusions.discarded.tsv`, and STAR's `SJ.out.tab` — not the BAMs, which stay on the server.

## E · Export the workflow (the reproducibility deliverable)

`Workflow → Extract workflow from history`, then `Download` → save the `.ga` file into this folder as
`dna_panel.ga` / `rna_fusion.ga`. A `.ga` file records every tool version and parameter, can be re-imported by
anyone with a Galaxy account, and can be published to WorkflowHub or Dockstore.
