# 5 · Finding variants: why three callers instead of one

**What the step does.** Runs several independent somatic callers over the same BAM and records, for every
position, which of them saw a variant and at what allele fraction (`scripts/merge_callers.py`).

| Caller | Model | Strength |
|---|---|---|
| GATK **Mutect2** (tumour-only) | haplotype reassembly + somatic likelihood | indels, read-orientation (FFPE) model |
| **VarDict** | local realignment, amplicon-aware | designed for deep targeted panels, low allele fractions |
| **LoFreq** | quality-aware Poisson-binomial | very low-frequency substitutions |
| **FreeBayes** | Bayesian haplotype | sensitive; needs filtering |

**Why concordance.** A benchmark on tumour-only amplicon panel data against a commercial reference standard
found that no single caller dominates: FreeBayes, VarScan, Mutect2 and Pisces each performed best on some
variant classes, and the authors' recommendation was explicit — *combine several callers with a concordance
metric and flexible depth/VAF thresholds* rather than trusting one. FreeBayes alone produced 2,433 calls on a
sample where Platypus produced 55; agreement between callers is what separates signal from each caller's own
failure modes.

The pipeline therefore records `n_callers` per variant, and the tumour-only filter treats a single-caller call
at low allele fraction as low quality rather than as a finding.

**Mutect2 in tumour-only mode** (used on Galaxy for the patient sample):

```bash
gatk Mutect2 -R GRCh38.fa -I tumour.bam -L panel.bed \
     --germline-resource af-only-gnomad.hg38.vcf.gz \
     --f1r2-tar-gz f1r2.tar.gz -O unfiltered.vcf
gatk LearnReadOrientationModel -I f1r2.tar.gz -O orientation.tar.gz
gatk FilterMutectCalls -R GRCh38.fa -V unfiltered.vcf --ob-priors orientation.tar.gz -O filtered.vcf
```
The read-orientation model is the FFPE defence: formalin deaminates cytosine on one strand, so the resulting
C>T changes appear preferentially on reads of one orientation. GATK recommends this filter for *all* FFPE and
NovaSeq samples. Note that `FilterMutectCalls` is not available as a Galaxy tool, so it is run locally on the
VCF that Galaxy produces.

**How it can mislead.** Requiring agreement costs sensitivity: a true subclonal variant seen by one caller only
will be demoted. That is a deliberate trade — on a tumour-only sample with no matched normal, precision is worth
more than the last few percent of recall, and every demoted call is still listed with its reason.
