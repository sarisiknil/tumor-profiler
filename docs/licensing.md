# Data sources, licences and what may be redistributed

| Source | Licence / access | Used how | May we redistribute? |
|---|---|---|---|
| **Ensembl VEP REST** | free, no key; ~55,000 requests/hour | annotation of every variant | results are derived data; the API is public |
| **ClinVar** (via VEP) | public domain (NIH) | clinical significance | yes |
| **gnomAD frequencies** (via VEP) | CC0 / open | germline filtering | yes |
| **CIViC** | CC0 | evidence items, tiering | **yes** — this is why CIViC is the backbone |
| **OncoKB** | free for academic research with a personal token; three gene lists downloadable without a licence | oncogenicity, therapeutic levels | **no** — annotations must not be redistributed; written only to `results/`, which is git-ignored. The token lives in `ONCOKB_TOKEN`, never in the repository |
| **COSMIC** | free for academics, but excludes for-profit settings and forbids exposing the data on a free-to-access platform | **not used** | n/a |
| **DGIdb** | open | drug–gene interactions | yes |
| **ClinicalTrials.gov API v2** | public domain | recruiting trials | yes |
| **Sanchez-Vega et al. 2018 pathway gene lists** | from the paper's supplement | pathway mapping | gene symbols only, with citation |
| **GRCh38 sequence** (Ensembl REST) | open | ~230 kB example reference | yes |
| **Galaxy (usegalaxy.eu)** | free public service, 250 GB per account | alignment, calling, STAR/Arriba | workflows exported as `.ga` |

**Practical consequences for this repository**

1. `resources/cancer_genes.tsv` is *fetched* by `resources/fetch_resources.sh` at setup time and git-ignored,
   rather than committed.
2. OncoKB annotations appear only under `results/`, never in `examples/` or in any committed output.
3. COSMIC is not queried, so no COSMIC-derived table can leak into a public dashboard.
4. The internship is hosted by a company. OncoKB and COSMIC both draw a line between academic and commercial
   use; work done under a company rather than the university may fall on the commercial side. The pipeline is
   built so that removing the OncoKB token leaves it fully functional on CIViC alone.
