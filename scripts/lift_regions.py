#!/usr/bin/env python3
"""Translate coordinates on the example's region contigs (`chr12_25204800_25250900`) back to real GRCh38
coordinates, so that annotation services see genuine chromosome positions.

This mirrors a real-world concern: a pipeline's coordinate system must be stated explicitly, and every
annotation step must agree on it. Works on VCF and on TSV files with chrom/pos columns.

Usage: lift_regions.py --regions examples/reference/regions.tsv -i in.vcf -o out.vcf
"""
import argparse, sys
from pathlib import Path

def load(regions):
    m = {}
    with open(regions) as fh:
        hdr = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4:
                m[f[0]] = (f[1], int(f[2]))
    return m

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True)
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    m = load(a.regions)
    is_vcf = a.input.endswith(".vcf")
    n = 0
    with open(a.input) as fin, open(a.out, "w") as fout:
        if is_vcf:
            for line in fin:
                if line.startswith("##contig"):
                    continue                      # contig lines refer to the mini reference
                if line.startswith("#"):
                    fout.write(line); continue
                f = line.rstrip("\n").split("\t")
                if f[0] in m:
                    chrom, off = m[f[0]]
                    f[0], f[1] = chrom, str(int(f[1]) + off - 1)
                    n += 1
                fout.write("\t".join(f) + "\n")
        else:
            hdr = fin.readline().rstrip("\n").split("\t")
            fout.write("\t".join(hdr) + "\n")
            ci, pi = hdr.index("chrom"), hdr.index("pos")
            for line in fin:
                f = line.rstrip("\n").split("\t")
                if f[ci] in m:
                    chrom, off = m[f[ci]]
                    f[ci], f[pi] = chrom, str(int(f[pi]) + off - 1)
                    if "key" in hdr:
                        ki = hdr.index("key")
                        parts = f[ki].split(":")
                        if len(parts) >= 3:
                            f[ki] = f"{chrom}:{f[pi]}:{parts[2]}"
                    n += 1
                fout.write("\t".join(f) + "\n")
    print(f"lifted {n} records -> {a.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
