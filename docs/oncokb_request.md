# Requesting an OncoKB academic licence

Register at <https://www.oncokb.org/account/register> with your **institutional** address (`...@sabanciuniv.edu`). The form asks for a short description of intended use. Text you can adapt:

## Short version (for a small form field)

> Undergraduate research project at Sabancı University (Faculty of Engineering and Natural Sciences). I am
> building an open-source, educational pipeline that interprets somatic variants from a targeted tumour
> sequencing panel. I would use the OncoKB Web API to annotate a few hundred variants with oncogenicity and
> therapeutic evidence levels, feeding an AMP/ASCO/CAP and ESCAT tiering step. Non-commercial and educational
> only; no OncoKB data will be redistributed and none will be used for patient care.

## Longer version (if there is room)

> I am an undergraduate student at Sabancı University, Faculty of Engineering and Natural Sciences, carrying out
> a required internship project (course XX395) whose deliverables are a written report submitted to the
> university and an open-source software repository.
>
> The project is a reproducible bioinformatics pipeline that takes tumour DNA and RNA sequencing data from a
> targeted cancer panel and turns it into an interpretable summary: read quality control, UMI-aware
> deduplication, tumour-only somatic variant calling, annotation, clinical tiering and a small visualisation
> dashboard. Its purpose is educational — to make each step of a clinical-grade workflow explicit and
> inspectable — rather than diagnostic.
>
> I would use the OncoKB Web API programmatically, through the `/annotate/mutations/byProteinChange` endpoint,
> to obtain the oncogenicity classification and therapeutic evidence level (1–4, R1/R2) for individual protein
> changes. Those annotations feed a rule-based tiering module that implements AMP/ASCO/CAP (Li et al. 2017) and
> ESCAT (Mateo et al. 2018). OncoKB matters here specifically because it records tumour-type-specific
> indications that the other free knowledge bases I use (CIViC, ClinVar, DGIdb) do not: my current tiering can
> only reach Tier II for a variant whose approved indication exists in the patient's tumour type, because the
> open sources carry that evidence only for other tumour types.
>
> Scale is very small: one anonymised clinical sample and one synthetic example dataset, a few hundred variants
> in total, a few hundred API calls in the lifetime of the project.
>
> Regarding the terms: no OncoKB data will be redistributed. The public GitHub repository contains code only;
> OncoKB responses are written to a results directory that is excluded from version control, and the pipeline
> runs fully without an OncoKB token (falling back to CIViC alone). I will not create a local copy of the
> database, will not use OncoKB data to train machine-learning models, and will not use it for clinical care or
> patient reporting. Both OncoKB papers (Chakravarty et al. 2017; Suehnholz et al. 2024) will be cited in the
> report and in the repository.

## Be accurate about the academic/commercial line

The internship is hosted by a company, and OncoKB's terms treat *research in a commercial setting* as
commercial use requiring a paid licence. Do not describe the project as purely academic if the company intends
to use the pipeline commercially. What is true, and what the text above says, is that:

- the project is coursework for a Sabancı University degree, with the report submitted to the university;
- the repository is open-source (MIT) and contains no OncoKB-derived data;
- OncoKB annotations are used only to produce the academic analysis, not a product.

If the company later wants to use the pipeline in a product or service, that needs its own OncoKB licence. If
you are unsure how your supervisor sees this, ask before submitting the form, and say so in the description —
OncoKB's team answer these questions directly.

## After approval

```bash
export ONCOKB_TOKEN=<your token>          # never commit it; scripts/annotate_evidence.py reads it from the env
snakemake -c4 --config mode=patient       # re-runs annotation and tiering with OncoKB evidence included
```

The pipeline is designed to work without the token, so nothing is blocked while you wait: `annotate_evidence.py`
simply skips OncoKB and records `oncokb_token_present: false` in its summary.

Three OncoKB files are downloadable with **no** licence at all, and the pipeline already uses one of them —
`resources/fetch_resources.sh` pulls the Cancer Gene List. Those are the Cancer Gene List, the All Curated Genes
List and the Biomarker-Drug Association List.
