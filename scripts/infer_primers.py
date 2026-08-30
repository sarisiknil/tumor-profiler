#!/usr/bin/env python3
"""Reverse-engineer the gene-specific primer (GSP2) set of a one-sided (AMP-style) panel from Read 2.

In anchored multiplex PCR every Read 2 begins with a gene-specific primer, so the most frequent Read-2 prefixes
ARE the primers. We count k-mer prefixes (default k=25, > typical primer length so that mapping is unambiguous)
and write (a) a TSV with counts, (b) a FASTA of candidates (count >= --min-count) that can be mapped to GRCh38
(e.g. BWA-MEM on Galaxy) to obtain primer coordinates/strand and, via gene annotation, the panel gene list.

Usage: infer_primers.py --r2 R2.fastq.gz [-n 2000000] [-k 25] [--min-count 30] -o out_prefix
"""
import argparse, gzip, collections

def open_fq(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r2", required=True); ap.add_argument("-n", type=int, default=2000000)
    ap.add_argument("-k", type=int, default=25); ap.add_argument("--min-count", type=int, default=30)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    def low_complexity(seq):
        """Drop NextSeq/NovaSeq two-colour poly-G artefacts and other homopolymer/dinucleotide junk."""
        if len(set(seq)) <= 2:
            return True
        for b in "ACGT":
            if b * 10 in seq:
                return True
        return False

    cnt = collections.Counter(); n = 0; n_lowcomplex = 0
    with open_fq(a.r2) as fh:
        while n < a.n:
            h = fh.readline()
            if not h: break
            s = fh.readline().rstrip("\n"); fh.readline(); fh.readline()
            n += 1
            if len(s) >= a.k and "N" not in s[:a.k]:
                if low_complexity(s[:a.k]):
                    n_lowcomplex += 1
                    continue
                cnt[s[:a.k]] += 1
    tot = sum(cnt.values())
    cands = [(km, c) for km, c in cnt.most_common() if c >= a.min_count]
    with open(a.out + ".tsv", "w") as t, open(a.out + ".fa", "w") as f:
        t.write("rank\tkmer\tcount\tfrac_of_reads\tcum_frac\n")
        cum = 0
        for i, (km, c) in enumerate(cands, start=1):
            cum += c
            t.write(f"{i}\t{km}\t{c}\t{c/tot:.5f}\t{cum/tot:.4f}\n")
            f.write(f">gsp2_{i:05d} count={c}\n{km}\n")
    covered = sum(c for _, c in cands) / tot if tot else 0
    print(f"reads scanned: {n:,} ({n_lowcomplex:,} low-complexity/poly-G discarded); distinct {a.k}-mers: {len(cnt):,}; "
          f"candidates (>= {a.min_count} reads): {len(cands):,} covering {covered*100:.1f} % of usable reads "
          f"-> {a.out}.tsv / .fa")

if __name__ == "__main__":
    main()
