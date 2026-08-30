# 1 · Reading the library before reading the biology

**What the step does.** Counts, for the first few hundred thousand reads, the base composition at each position,
the read-length distribution, the adapter content, and how concentrated the read start positions are.
`scripts/read_structure.py`, no reference genome needed, runs in seconds.

**Why it comes first.** A FASTQ file does not say which kit produced it, and the correct processing depends
entirely on that. The composition profile answers it directly:

```
position:  1 ......... 12 | 13 .......... 26 | 27 ...............
Read 1  :  random bases   | one fixed base    | patient sequence
           = molecular      per position      = the insert
             barcode        = adapter's
             (UMI)            common region
Read 2  :  gene-specific primer | patient sequence
```

A position where all four bases occur at ~25 % is random; a position where one base occurs at >98 % is a fixed
synthetic sequence. The transition between the two marks the end of the barcode.

**What we found in this sample.** DNA: 12-nt UMI, then `AGTCGTCTCGAAG` + `T`. RNA: 12-nt UMI, then
`CTGGATAGTACGCT`. Read 2 begins at a gene-specific primer in both. Index reads are 8 nt (i7) + 10 nt (i5).
Those four facts together identify the chemistry as **IDT/Archer anchored multiplex PCR** with the
"Liquid P5 MBC + P7" adapters, whose product insert specifies a 12-nt molecular barcode with a 10-nt i5 index
and an 8-nt i7 index. QIAseq is excluded (its UMI is on Read 2), CleanPlex and Pillar are excluded (primers at
*both* read starts), and hybrid-capture assays such as TSO500 are excluded (no primer signature at all).

**How it can mislead.** The heuristic needs a fixed region after the barcode; on a library without one it must
report "no inline UMI" rather than inventing a length — a mistake that would silently corrupt every downstream
step. The code does exactly that, and the unit test in `tests/test_scripts.py` pins the behaviour.

**Read further.** Zheng et al. 2014 (anchored multiplex PCR, *Nature Medicine*); the ArcherDX molecular-barcode
technical note; `literature/` in this project.
