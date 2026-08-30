# 10 · From a gene list to a picture of the cell

**What the step does.** Maps every altered gene onto the ten canonical oncogenic signalling pathways curated by
Sanchez-Vega et al. (2018) across 9,125 TCGA tumours — cell cycle, Hippo, MYC, Notch, NRF2, PI3K, RTK-RAS,
TGF-β, TP53 and WNT — plus two supplementary sets those ten deliberately exclude: DNA-damage repair/mismatch
repair, and chromatin regulators. The gene sets live in `resources/oncogenic_pathways.gmt`.

**Why it helps.** Drugs act on pathways. Two patients with different mutations in *EGFR* and *KRAS* share a
broken RTK-RAS pathway, and that shared fact is often what determines the therapeutic strategy. In the TCGA
analysis, 89 % of tumours carried at least one driver alteration in these ten pathways and 57 % carried at least
one potentially targetable alteration.

**What this is not.** It is a *membership map*, not an enrichment analysis, and the distinction is not
pedantic:

- Enrichment asks whether a pathway contains more altered genes than chance predicts. That requires a gene
  universe — and here the universe was chosen by whoever designed the panel, which is precisely the set of genes
  most likely to be mutated in cancer. Any enrichment test would be circular.
- A single sample provides no null distribution.
- **A pathway with no hit may simply have no genes on the panel.** Absence of evidence, not evidence of absence.

The dashboard states this on the page itself, because a bar chart of pathways looks like a statistical result
whether or not it is one.
