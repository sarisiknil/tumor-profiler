# 6 · Tumour or inherited? The central problem of a sample without a matched normal

**The problem.** Every cell in the body carries ~4–5 million germline variants. A tumour adds a comparatively
small number of somatic mutations. With only tumour DNA sequenced, both appear in the same VCF, and the
distinction is not visible in the reads themselves.

**What the pipeline does** (`scripts/filter_tumor_only.py`), in order, labelling rather than deleting:

1. **Population frequency.** A variant seen in gnomAD above 0.1 % is almost certainly inherited (this is PCGR's
   default tumour-only threshold). → `GERMLINE_LIKELY`.
2. **Allele fraction plus database evidence.** ~50 % allele fraction, a dbSNP identifier, and an *unambiguously*
   benign ClinVar record. The word "unambiguously" matters: VEP aggregates the significance of every ClinVar
   record at a position, so the string `benign,likely_benign,likely_pathogenic,pathogenic` must not be read as
   "benign" — an earlier version of this filter made exactly that mistake and wrongly discarded a pathogenic
   *PTEN* nonsense variant.
3. **Quality.** Depth, alternate-read count, allele fraction and caller agreement. → `LOW_QUALITY`.
4. **Artefact pattern.** Low-fraction C>T/G>A in a single caller, the signature of formalin damage.
   → `ARTEFACT_LIKELY` (flagged, never hard-filtered when the variant is actionable — the SOBDetector authors
   are explicit that a probability should be reported alongside a targetable mutation, not used to remove it).
5. **Ambiguity flag.** Heterozygous allele fraction, rare in the population, in a cancer-predisposition gene or
   ClinVar-pathogenic → `germline_ambiguous = yes`. This is the honest answer, and it has consequences beyond
   the pipeline (see [12_ethics.md](12_ethics.md)).

**What cannot be fixed by filtering.** Population databases cannot contain a variant that is private to one
family. Sun et al. (2018) built the reference method for this problem, SGZ, which models each variant's expected
allele fraction from tumour purity, ploidy and local copy number; it reaches 95–99 % accuracy but only assigns a
call for 83–85 % of variants, needs >500× depth over a genome-wide SNP backbone, and fails above ~90 % tumour
purity. A small amplicon panel does not provide that backbone. The residual germline fraction is therefore a
property of the assay, not a bug in the code, and the report states it.

**The correct fix**, when it is available: sequence the patient's blood alongside the tumour.
