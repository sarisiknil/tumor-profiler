#!/usr/bin/env python3
"""Detect clinically important exon-skipping events from STAR's splice-junction table (SJ.out.tab).

Motivation: MET exon 14 skipping is an approved drug target (capmatinib, tepotinib) and is a *splicing* event,
not a fusion, so fusion callers such as Arriba and STAR-Fusion do not report it (Capone et al. 2022). RNA-based
detection is more reliable than DNA-based detection of the underlying intronic variants (Davies et al. 2019).

Method: for each configured skipping event we count reads supporting the canonical junctions flanking the exon
versus reads supporting the junction that skips it, then report a skipping ratio. Genomic coordinates are
GRCh38. The approach generalises: add an entry to EVENTS to screen another exon.

Input : STAR SJ.out.tab (col1 chrom, col2 intron start, col3 intron end, col7 uniquely-mapped reads)
Output: <out>.tsv and <out>.json
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

# GRCh38 intron coordinates, derived from the Ensembl canonical transcript and pinned here so the tool needs
# no network access. `tests/test_exon_skipping_coords.py` re-derives them from Ensembl and fails if they drift -
# a wrong coordinate would silently report "not detected", which is the worst possible failure mode for a
# targetable event.
#   MET  ENST00000397752: exon 13 ends 116771654, exon 14 = 116771849-116771989 (141 bp), exon 15 starts 116774881
#   EGFR ENST00000275493: exon 1 ends 55019365, exon 2 starts 55142286, exon 8 starts 55155830
EVENTS = {
    "MET_exon14_skipping": {
        "chrom": "chr7", "gene": "MET", "transcript": "ENST00000397752",
        "canonical_junctions": [(116771655, 116771848), (116771990, 116774880)],
        "skipping_junction": (116771655, 116774880),
        "significance": "Targetable by MET inhibitors (capmatinib, tepotinib); ~3-4 % of lung adenocarcinoma. "
                        "RNA detects it more reliably than DNA, because the causal variants are scattered "
                        "across the flanking introns and splice sites (Davies et al. 2019). This event is "
                        "explicitly within the FusionPlex Lung v2 assay's scope.",
        "tolerance": 12,
    },
    "EGFRvIII_exon2-7_deletion": {
        "chrom": "chr7", "gene": "EGFR", "transcript": "ENST00000275493",
        "canonical_junctions": [(55019366, 55142285)],
        "skipping_junction": (55019366, 55155829),
        "significance": "EGFRvIII: in-frame loss of exons 2-7, characteristic of glioblastoma.",
        "tolerance": 20,
    },
}

def load_sj(path):
    js = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        chrom = f[0] if f[0].startswith("chr") else "chr" + f[0]
        try:
            js.append((chrom, int(f[1]), int(f[2]), int(f[6]), int(f[7]) if len(f) > 7 else 0))
        except ValueError:
            continue
    return js

def count_near(js, chrom, start, end, tol):
    return sum(u + m for c, s, e, u, m in js if c == chrom and abs(s - start) <= tol and abs(e - end) <= tol)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sj", required=True, help="STAR SJ.out.tab")
    ap.add_argument("--min-reads", type=int, default=3)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    js = load_sj(a.sj)
    rows = []
    for name, ev in EVENTS.items():
        tol = ev["tolerance"]
        skip = count_near(js, ev["chrom"], *ev["skipping_junction"], tol)
        canon = [count_near(js, ev["chrom"], s, e, tol) for s, e in ev["canonical_junctions"]]
        canon_max = max(canon) if canon else 0
        ratio = skip / (skip + canon_max) if (skip + canon_max) else 0.0
        rows.append({"event": name, "gene": ev["gene"], "skipping_reads": skip,
                     "canonical_reads": ";".join(map(str, canon)), "canonical_max": canon_max,
                     "skipping_ratio": round(ratio, 3),
                     "called": "yes" if skip >= a.min_reads and ratio >= 0.10 else "no",
                     "significance": ev["significance"]})
    write_tsv(rows, a.out + ".tsv",
              ["event","gene","called","skipping_reads","canonical_max","skipping_ratio","canonical_reads","significance"])
    write_json({"events": rows, "min_reads": a.min_reads, "junctions_in_file": len(js),
                "method": "Counts of uniquely+multi-mapped reads spanning the skipping junction versus the "
                          "canonical flanking junctions in STAR's SJ.out.tab (GRCh38 coordinates).",
                "caveats": ["Coordinates are transcript-model dependent; verify against the panel's target "
                            "transcript before reporting.",
                            "A targeted RNA panel only covers primed regions: absence of evidence here is not "
                            "evidence of absence unless the exon is covered.",
                            "Positive calls require orthogonal confirmation (RT-PCR) before clinical use."]},
               a.out + ".json")
    for r in rows:
        print(f"{r['event']}: skip={r['skipping_reads']} canonical={r['canonical_max']} "
              f"ratio={r['skipping_ratio']} called={r['called']}", file=sys.stderr)

if __name__ == "__main__":
    main()
