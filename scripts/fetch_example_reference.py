#!/usr/bin/env python3
"""Build a *tiny* but REAL GRCh38 reference for the public example dataset.

Instead of shipping a synthetic genome (which no annotation service understands), we fetch a few dozen kilobases
of genuine GRCh38 sequence around well-known driver-gene hotspots through the Ensembl REST API. The result is
~250 kB, indexes in a second, and — crucially — every coordinate is a real genomic coordinate, so VEP, ClinVar,
gnomAD, OncoKB and CIViC all work on the example exactly as they do on a real sample.

Contigs are named `<chrom>_<start>_<end>` so that a lift-over back to genome coordinates is explicit and
teachable (see scripts/lift_regions.py); `regions.tsv` records the mapping.

Usage: fetch_example_reference.py -o examples/reference
"""
import argparse, json, sys, time, urllib.request
from pathlib import Path

# GRCh38 windows around canonical hotspots (chrom, start, end, gene, note)
REGIONS = [
    ("12", 25204800, 25250900, "KRAS",  "exon 2-3; G12/G13/Q61 hotspots (chr12:25245350 = G12C)"),
    ("7", 140730000, 140784000, "BRAF",  "exon 11-15; V600E at chr7:140753336"),
    ("17", 7668000, 7690000,   "TP53",  "whole coding region; R175H at chr17:7675088"),
    ("7", 55174000, 55200000,  "EGFR",  "exon 18-21; L858R at chr7:55191822, exon 19 deletions"),
    ("10", 87925000, 87975000, "PTEN",  "exon 5-9; R130* at chr10:87933147"),
    ("3", 179198000, 179235000,"PIK3CA","exon 9-20; E545K at chr3:179218303, H1047R at chr3:179234297"),
]
SERVER = "https://rest.ensembl.org"

def fetch(chrom, start, end, tries=4):
    url = f"{SERVER}/sequence/region/human/{chrom}:{start}..{end}?coord_system_version=GRCh38"
    req = urllib.request.Request(url, headers={"Content-Type": "text/x-fasta",
                                               "User-Agent": "tumor-profiler/1.0"})
    for i in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                txt = r.read().decode()
            return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">")).upper()
        except Exception as e:
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="examples/reference")
    a = ap.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    fa, regions = out / "mini_GRCh38.fa", out / "regions.tsv"
    if fa.exists() and regions.exists():
        print(f"[skip] {fa} already exists"); return
    with open(fa, "w") as f, open(regions, "w") as t:
        t.write("contig\tchrom\tstart\tend\tgene\tnote\n")
        for chrom, start, end, gene, note in REGIONS:
            seq = fetch(chrom, start, end)
            name = f"chr{chrom}_{start}_{end}"
            f.write(f">{name} {gene} GRCh38 chr{chrom}:{start}-{end}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
            t.write(f"{name}\tchr{chrom}\t{start}\t{end}\t{gene}\t{note}\n")
            print(f"[ok] {gene:7s} {name} ({len(seq):,} bp)")
            time.sleep(0.4)
    print(f"wrote {fa} and {regions}")

if __name__ == "__main__":
    main()
