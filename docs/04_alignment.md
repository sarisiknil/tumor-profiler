# 4 · Alignment and coverage: which parts of the genome can we actually judge?

**What the step does.** Aligns the UMI-stripped reads to GRCh38 with `bwa mem`, deduplicates by barcode, clips
primer bases, and measures per-target depth with `mosdepth`.

**Where it runs.** The human genome index needs about 6 GB of RAM for `bwa` (and ~10 GB for `bwa-mem2`, whose
index alone is 10 GB on disk), so for the patient sample this step runs on **usegalaxy.eu**, which gives every
registered account 250 GB of storage and pre-built GRCh38 indexes. The example dataset uses a 230 kB
mini-reference and therefore runs on any laptop — that is what continuous integration executes.

```bash
bwa mem -Y -K 150000000 -t 4 -R '@RG\tID:S\tSM:S\tPL:ILLUMINA' GRCh38.fa R1.fq.gz R2.fq.gz | samtools sort -o aligned.bam
```
`-Y` keeps supplementary alignments soft-clipped (needed for fusion and structural evidence); `-K` fixes the
chunk size so the result is deterministic; `-M` must **not** be used, as it breaks downstream tools.

**Coverage defines the report's scope.** A variant can only be called where there are reads. The pipeline
computes the territory covered at ≥100× and calls it *assessable*; everything else is reported as **not
assessable**, which is a different statement from "no mutation found". This distinction is the single most
common misreading of a panel report.

**How it can mislead.** Amplicon coverage is extremely uneven — primer efficiency varies by an order of
magnitude — so a mean depth of 800× can hide targets at 20×. Always look at the per-target distribution, not
the mean.
