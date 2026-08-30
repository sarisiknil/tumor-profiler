# 7 · From a position to a meaning, and from a meaning to a decision

**Annotation** (`scripts/annotate_vep_rest.py`). Ensembl VEP is queried through its REST API rather than the
20 GB offline cache: a targeted panel produces a few hundred variants, the API accepts 180 at a time and allows
~55,000 requests an hour, so the whole panel annotates in seconds on a laptop. For each variant the pipeline
keeps the consequence on the MANE Select transcript, the HGVS notations, gnomAD population frequency, ClinVar
significance, dbSNP identifier and the SIFT/PolyPhen predictions.

**Evidence** (`scripts/annotate_evidence.py`). The protein change is converted to one-letter notation
(`p.Val600Glu` → `V600E`) because that is how knowledge bases index variants, then:

- **CIViC** (GraphQL, CC0): every accepted evidence item for that variant — type, level A–E, significance,
  therapies and disease.
- **OncoKB** (free academic token in `ONCOKB_TOKEN`): oncogenicity and therapeutic levels 1–4/R1–R2.
- **DGIdb**: drug–gene interactions, a much weaker claim than clinical evidence.
- **ClinicalTrials.gov v2**: recruiting trials that genuinely mention the gene.
- **COSMIC is not queried at all** — free for academics, but the licence forbids redistribution and excludes
  for-profit settings, which would make the repository unshareable. See [licensing.md](licensing.md).

**Tiering** (`scripts/tiering.py`). Two systems, both rule-based and both printing their reasoning:

| | AMP/ASCO/CAP (Li et al. 2017) | ESCAT (Mateo et al. 2018) |
|---|---|---|
| asks | how strong is the clinical evidence? | how ready is this target for routine use? |
| I | approved drug / guideline, this tumour type | I-A/I-B ready for routine use |
| II | approved in another tumour type, or trial evidence | II investigational, clinical data exist |
| III | unknown clinical significance | III benefit shown in a different tumour type |
| IV | benign | IV preclinical only |

The two disagree on purpose, and seeing where they disagree is the point: *BRAF* V600E is tier I in melanoma and
tier II in a tumour type where the drug is not approved — the same molecule, a different decision.

**The disease context is part of the tier, so the pipeline enforces it.** Evidence is only allowed to support
Tier I if it comes from the patient's own tumour type; identical evidence from a different tumour type is capped
at Tier II, and the rationale column says which disease each evidence item came from. This has a practical
consequence worth understanding: a variant can be genuinely Tier I under a clinical guideline while the *freely
available* evidence only supports Tier II, simply because the open knowledge bases cover tumour types unevenly.
When that happens the honest output is Tier II with the reason attached — not a Tier I the evidence does not
carry.

**How it can mislead.** Automated tiering reproduces the *evidence*, not the judgement of a molecular tumour
board. It does not know the patient's prior treatments, cannot weigh a resistance mutation against a
sensitivity one, and inherits every gap and lag in the underlying databases.
