# What can we learn about a person's disease from their tumour's DNA and RNA?

This pipeline answers that question for one biopsy, and — just as importantly — shows where the answer stops.

## The chain of reasoning

```
 tumour biopsy
      │  DNA and RNA extracted, amplified with gene-specific primers (anchored multiplex PCR)
      ▼
 raw reads  ──► what assay is this?          (read structure: UMI, common region, primers)
      │        is the material intact?       (adapter read-through, insert size, depth)
      ▼
 alignment  ──► which molecules were unique? (UMI deduplication - not positional dedup)
      │        which regions can we judge?   (coverage: only >100x territory is "assessable")
      ▼
 variants   ──► what changed in the DNA?     (three callers; agreement is the filter)
      │        is it tumour or inherited?    (population frequency, ClinVar, allele fraction)
      │        is it an artefact?            (FFPE C>T deamination, strand-orientation bias)
      ▼
 meaning    ──► is this a known driver?      (VEP, ClinVar, cancer-gene list)
      │        does a drug target it?        (CIViC, OncoKB -> AMP/ASCO/CAP tier, ESCAT)
      │        which pathway is broken?      (Sanchez-Vega ten-pathway map)
      ▼
 RNA        ──► is there a fusion?           (STAR + Arriba; DNA panels miss most fusions)
      │        is an exon being skipped?     (splice-junction counting, e.g. MET exon 14)
      ▼
 report + dashboard
```

## What a single tumour-only panel *can* establish

| Question | Evidence in this data | Confidence |
|---|---|---|
| Which cancer genes are mutated? | variant calls in panel genes | high, within covered regions |
| How large is the mutant clone? | variant allele fraction | moderate; confounded by purity and copy number |
| Is a mutation a known driver? | ClinVar, cancer-gene list, hotspot position | high for recurrent hotspots |
| Is there an approved drug? | CIViC/OncoKB evidence levels, AMP/ASCO/CAP tier | high for tier I; the tumour type matters |
| Is there a targetable fusion? | RNA fusion calls | moderate; open callers are less sensitive than the vendor's |
| Which pathways are affected? | pathway membership of altered genes | qualitative only |

## What it cannot establish

- **Germline versus somatic.** No normal tissue was sequenced. See [06_tumour_only.md](06_tumour_only.md).
- **Mutational burden or signatures.** Too little territory; see [08_biomarkers.md](08_biomarkers.md).
- **Microsatellite instability**, without a panel-specific model.
- **Copy number, quantitatively**, without a cohort reference.
- **Tissue of origin or global expression** — the RNA assay is a small fusion panel, not RNA-seq.
- **Prognosis for this individual.** Every evidence item describes a cohort.

Each stage has its own page: [01](01_read_structure.md) read structure · [02](02_umi.md) UMIs ·
[03](03_primers.md) primers · [04](04_alignment.md) alignment · [05](05_variant_calling.md) calling ·
[06](06_tumour_only.md) tumour-only filtering · [07](07_annotation_tiering.md) annotation and tiers ·
[08](08_biomarkers.md) biomarkers · [09](09_rna_fusions.md) RNA and fusions · [10](10_pathways.md) pathways ·
[11](11_reproducibility.md) reproducibility · [12](12_ethics.md) ethics and data protection ·
[licensing.md](licensing.md) data licences.
