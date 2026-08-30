# 9 · What RNA adds that DNA cannot see

**Why sequence RNA at all.** A gene fusion is created by a DNA rearrangement whose breakpoint usually falls in a
large intron — tens of kilobases of sequence that no targeted DNA panel covers. The *transcript*, in contrast,
splices the two partners directly together, so a single read can span the junction. In lung adenocarcinomas
where DNA sequencing found no driver, adding RNA sequencing recovered targetable kinase fusions in a substantial
proportion of cases (Benayed et al. 2019); sequential DNA-then-RNA testing is now standard practice
(Cohen et al. 2020).

**How the fusion arm works.** UMI extraction (common region `CTGGATAGTACGCT`), then STAR with the chimeric
alignment parameters Arriba expects, then Arriba with its blacklist, known-fusion list and protein-domain
annotation. STAR's human genome index needs ~31 GB of RAM — the alignment step therefore runs on usegalaxy.eu,
where both STAR and Arriba 2.5.1 are installed.

**Why the RNA arm does *not* UMI-deduplicate before Arriba**, although the DNA arm does. With
`--chimOutType WithinBAM`, the chimeric evidence for a fusion lives in *supplementary* alignments inside the
main BAM. `umi_tools dedup` works on primary alignments and does not preserve those supplementary records
reliably, so deduplicating here can quietly delete the very reads the fusion call rests on. Arriba is therefore
run on the STAR output directly, as its own workflow intends, and its supporting-read counts are read counts
rather than molecule counts. The UMI is still in each read name, so a molecule-level count can be recovered
afterwards from the read names Arriba reports — the ordering is a deliberate trade, not an oversight.

`scripts/fusion_summary.py` then re-ranks Arriba's output by clinical relevance: whether either partner is a
kinase with approved inhibitors (ALK, ROS1, RET, NTRK1/2/3, MET, FGFR1-3, BRAF, NRG1 …), whether the fusion is
in frame, which protein domains are retained, and whether only one partner is on the panel — which is *expected*
for anchored multiplex PCR, since the assay primes one known gene and discovers the partner without prior
knowledge. That is the property that lets a fixed panel find fusion partners nobody has described before.

**What the alignment statistics look like, and why they are not alarming.** On this sample STAR reported
~21 % uniquely mapped reads and **~76 % chimeric reads**. On whole-transcriptome data that would signal a
serious problem. Here it is a property of the chemistry: Read 1 begins at an arbitrary ligation point and
Read 2 at a fixed gene-specific primer, inserts are short, and the Arriba preset sets `--chimSegmentMin 10`,
so a 10-base segment mapping elsewhere is enough to classify a read as chimeric. The consequence is practical
rather than cosmetic — the fusion caller receives an enormous pile of candidate chimeras, most of them
artefacts of library construction, which is exactly why the panel-membership check ("is either partner
actually a target of this assay?") does more filtering work here than Arriba's own confidence score.

**Screening the output for artefacts — the step that matters most on this chemistry.** On this sample Arriba
reported 20 fusions and discarded 166,694 more. Every one of the 20 involved a gene the assay primes, and the
pattern was unmistakable once the calls were laid side by side: *FGFR1–FGFR3*, *FGFR3–FGFR1*, *FGFR2–FGFR3*,
*FGFR3–FGFR2* — the same paralogous family joined in every possible combination, in **both orientations**, on a
handful of reads each. `scripts/fusion_summary.py` therefore screens for the patterns this chemistry produces:

| Flag | Why it indicates an artefact |
|---|---|
| reciprocal call (A–B *and* B–A reported) | a real fusion has one orientation |
| paralogue pair (FGFR1/2/3, NTRK1/2/3, …) | reads and PCR products cross between near-identical genes |
| both partners are panel targets | AMP primes one gene and *discovers* the partner; a join between two primed genes is a classic PCR chimera, since all targets are amplified in one tube |
| promiscuous partner (≥3 different partners) | an artefact hub, not a driver |
| partner with no gene symbol (ENSG…) | uncharacterised locus, usually a mapping artefact |
| 5'-5' orientation | the two 5' ends are joined; no functional protein is possible |
| ≤2 supporting reads | below any credible threshold |

This matters more than it may appear. The most clinically exciting fusion in cholangiocarcinoma is an **FGFR2
fusion** — a target with approved inhibitors. An unscreened reading of this output would have reported exactly
that, and it would have been wrong. The accredited laboratory reported no fusion in this sample, and after
screening the pipeline agrees. **Agreeing about a negative is a real result**, and the reasoning is recorded per
call rather than asserted.

**Honest sensitivity.** Benchmarked on real Archer FusionPlex FFPE samples against the vendor's software, Arriba
recovered 86 % of fusions on the lung panel and 57 % on the sarcoma panel; STAR-Fusion recovered 33 % and 7 %
(Capone et al. 2022). The pipeline therefore uses Arriba, keeps low-confidence calls rather than discarding
them, and states the gap instead of implying that the open-source route is equivalent. Any fusion that would
influence treatment needs orthogonal confirmation by FISH, IHC or RT-PCR.

**Exon skipping is a separate problem.** MET exon 14 skipping is targetable, common enough to matter in lung
cancer, and *invisible to fusion callers* — it is a splicing change within one gene, not a join between two.
`scripts/exon_skipping.py` counts reads across the skipping junction versus the canonical junctions in STAR's
`SJ.out.tab` and reports the ratio. RNA detects these events far more reliably than DNA, where the causal
variants are scattered across the flanking introns and splice sites.

**Expression from a fusion panel.** Only primed targets are measured, so per-gene counts are informative only in
relative terms, and single-sample pathway or tissue-of-origin inference is not defensible: published work on
targeted RNA panels shows that even lineage classification needs a supervised reference cohort and still reached
only ~81 % accuracy.
