# 12 · Ethics, consent and data protection

**The data.** One real tumour biopsy from a clinical genetics laboratory, analysed with the patient's consent for
analysis and for the publication of anonymised findings in a student report.

**Pseudonymised is not anonymous.** Under the GDPR, and under the Turkish KVKK — which classes genetic and health
data as special-category data requiring explicit consent — removing a name does not take genomic data out of
scope. Genomes are identifying by nature: published work has shown re-identification from as few as ~25
independent SNP genotypes. The operative question, per GDPR Recital 26, is whether identification is possible by
"means reasonably likely to be used", and for raw sequence data the honest answer is yes.

**What follows in practice.** The project therefore separates two categories:

| Never leaves the local machine | May appear in the report / dashboard |
|---|---|
| FASTQ, BAM, CRAM, gVCF | somatic variants, with the sample under an alias |
| the laboratory's sample identifier and dates | fusion calls, pathway summaries, tiering |
| germline variant lists, common-SNP genotypes | aggregate QC statistics |
| the clinical report from the laboratory | method descriptions and runtimes |

This mirrors the practice of large cancer consortia, which treat tumour-specific somatic variants as effectively
anonymous and distribute them openly, while germline data stay under controlled access. `.gitignore` enforces the
left-hand column, the private clinical report lives outside the repository, and
`tests/test_pipeline_outputs.py` asserts that the laboratory identifier appears in no shareable output.

**Incidental germline findings.** The pipeline flags variants it cannot assign to tumour or germline —
heterozygous, rare, pathogenic, in a cancer-predisposition gene such as *BRCA1/2*, *TP53* or *PTEN*. If such a
variant were inherited it would matter for the patient's relatives. That is a matter for consent and genetic
counselling, not for a filtering rule, and the correct action for a student project is to flag it, report it as
ambiguous, and not to pursue it further.

**What this pipeline is not.** It is an educational re-analysis. It is not accredited, not validated against a
reference standard, and not a diagnostic device. Every output carries that statement, and the clinically
reported results from the accredited laboratory remain the authoritative ones — the comparison against them is
part of the project's validation, not a challenge to them.
