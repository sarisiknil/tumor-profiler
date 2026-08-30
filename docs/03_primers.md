# 3 · Recovering the panel design from the reads

**What the step does.** Counts the most frequent 25-base prefixes of Read 2 (`scripts/infer_primers.py`), writes
them as a FASTA, and — after mapping them to GRCh38 — turns them into a primer BED file.

**Why.** The panel's target list is a trade secret, but it is *implicit in the data*: in anchored multiplex PCR
every Read 2 starts at a gene-specific primer. Counting read-2 prefixes therefore reconstructs the primer set,
and mapping those primers reconstructs the panel's genes. In this sample about 2,200 distinct prefixes account
for 80 % of the DNA reads and only 58 for 80 % of the RNA reads — the RNA assay is a small fusion panel, and
that single number changes what the RNA arm is allowed to claim.

Low-complexity prefixes are discarded first: on two-colour instruments (NextSeq, NovaSeq) a dark cluster is
called as `G`, so poly-G runs are an artefact of the sequencer, not a primer.

**Why primer bases must be clipped.**

```
 reference   ...ACGTTGCAGGATCCATTAGGC...
 primer      ...ACGTTGCAGGATCC          <- synthesised oligo: always matches the design, never the patient
 read 2      ...ACGTTGCAGGATCCATTAGGC...
                └── these bases carry no patient information ──┘
```

If a patient's true variant lies under the primer, the primer overwrites it and the variant disappears; if the
primer was designed against a different allele, it can create a false variant at high allele fraction.
`samtools ampliconclip -b primers.bed --hard-clip` removes them before calling.

**How it can mislead.** Inferred primers are *candidates*: sequencing errors create near-duplicates, and a primer
that overlaps a repetitive region can map to several places. The reconstructed panel is good enough to define
which genes are assessable; it is not the manufacturer's design file, and the report says so.
