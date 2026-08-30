# 2 · Unique molecular identifiers: counting molecules, not reads

**What the step does.** Moves the 12-nt barcode from the start of Read 1 into the read name, discards the fixed
common region, and later collapses reads that share a barcode and a mapping position into one molecule.

```bash
umi_tools extract --extract-method=regex \
  --bc-pattern='^(?P<umi_1>.{12})(?P<discard_1>AGTCGTCTCGAAGT?){s<=2}' \
  -I R1.fastq.gz --read2-in R2.fastq.gz -S out_R1.fastq.gz --read2-out out_R2.fastq.gz
umi_tools dedup -I aligned.bam -S dedup.bam --paired --method=directional
```

**Why a regex and not the simple string pattern.** In UMI-tools' string syntax, `N` marks barcode bases and `X`
marks bases that are *kept and re-attached to the read*. `NNNNNNNNNNNNXXXXXXXXXXXXX` would therefore leave the
adapter's common region inside the read, where it cannot align and drags the alignment score down. The regex
form has an explicit `discard_` group; `{s<=2}` tolerates up to two sequencing errors in the fixed region.

**Why deduplication must use the barcode.** In an amplicon assay every molecule from one primer starts at the
same coordinate, so position-based duplicate marking treats independent molecules as copies of each other. On
commercial UMI panels this pushed the apparent unique-molecule rate from 10 % up to 52 % once barcodes were used
(Kim et al. 2019, *BMC Genomics*). The ratio of reads to molecules is also the honest measure of library
complexity: 5 million reads over 300,000 molecules is a 300,000-molecule experiment.

**Going further: consensus reads.** Deduplication keeps one representative read per molecule. `fgbio` instead
*builds* a consensus base by base across the family, which suppresses PCR and sequencing errors and is what
makes variant detection below ~1 % allele fraction credible (Salk et al. 2018, *Nature Reviews Genetics*):

```bash
fgbio FastqToBam --read-structures 12M14S+T +T -i R1.fq.gz R2.fq.gz -o unmapped.bam -s SAMPLE
# align, then:
fgbio GroupReadsByUmi --strategy Identity   # Identity, not Adjacency: amplicon start sites are fixed
fgbio CallMolecularConsensusReads --min-reads 1
fgbio FilterConsensusReads --min-reads 3 --min-base-quality 45 --max-base-error-rate 0.2
```

**How it can mislead.** A 12-nt barcode has 16.7 million combinations, but collisions still happen at high depth;
`--method=directional` accounts for barcode errors, `--method=unique` does not. And deduplication cannot rescue a
library that was low-complexity to begin with: if 5 million reads collapse to 20,000 molecules, no algorithm can
recover the information that was never there.
