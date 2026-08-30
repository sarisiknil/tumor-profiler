# 8 · Biomarkers: numbers that are easy to compute and easy to misuse

## Tumour mutational burden

TMB is the count of nonsynonymous somatic mutations divided by the assessable territory in megabases. It is used
because tumours with many mutations present more neoantigens and respond better to immune-checkpoint inhibitors;
≥10 mutations/Mb is the threshold in the tissue-agnostic pembrolizumab indication, measured with a large
commercial panel.

Three things make a panel TMB fragile, and the pipeline reports all three:

1. **Assay dependence.** The Friends of Cancer Research TMB Harmonization Project showed that the same tumour
   yields different TMB values on different panels; values are not portable without harmonisation
   (Merino et al. 2020).
2. **Germline leakage.** Every inherited variant that survives tumour-only filtering inflates the count.
3. **Counting statistics.** Over 1 Mb, five mutations gives a 95 % confidence interval of roughly 1.6–11.7 —
   spanning the clinical threshold. The pipeline prints the interval next to the point estimate, and refuses to
   treat the number as meaningful when the assessable territory is below 0.1 Mb.

## What the pipeline deliberately does not compute

- **Microsatellite instability.** MSIsensor2 needs a model trained for the specific panel; MSIsensor-pro's
  normal-free mode needs a baseline built from ≥20 normal samples on the same assay. Neither exists here, and a
  number produced without them would be unfalsifiable.
- **Mutational signatures.** Refitting the COSMIC SBS catalogue needs on the order of 100–200 SNVs (PCGR's
  default minimum is 200). A panel yielding a few dozen cannot support it. Reporting "signature 3, homologous
  recombination deficiency" from ten mutations would be the single most dangerous thing this pipeline could do.

The FFPE indicator *is* computed: the share of C>T/G>A among somatic SNVs, and how many of those sit below 10 %
allele fraction. It is a quality-control signal about the specimen, not a biological finding.
