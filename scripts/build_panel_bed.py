#!/usr/bin/env python3
"""Turn a panel's gene list into a BED of exonic target regions, using the Ensembl REST API.

Why: knowing exactly which genomic territory an assay interrogates is what makes a negative result meaningful
and what gives the TMB calculation an honest denominator. The manufacturer's design file is not public, but the
gene list is printed in the assay documentation, and the MANE Select transcript's exons are the right
approximation of what a panel targets.

Usage:
  python3 scripts/build_panel_bed.py --genes resources/panel_variantplex_expanded_solid_tumor.txt \
                                     -o resources/panel_variantplex_expanded_solid_tumor
Produces <out>.bed (exons, padded), <out>_genes.tsv (per-gene span and exon count) and <out>_summary.json.
"""
import argparse, json, sys, time, urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

ENSEMBL = "https://rest.ensembl.org"
UA = {"User-Agent": "tumor-profiler/1.0 (educational pipeline)", "Accept": "application/json"}


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90) as r:
                return json.load(r)
        except Exception as e:
            if i == tries - 1:
                print(f"    request failed: {e}", file=sys.stderr)
                return None
            time.sleep(2 * (i + 1))


# HGNC has renamed some genes since the assay documentation was written; keep the assay's label but resolve
# the current symbol.
ALIASES = {"H3F3A": "H3-3A", "HIST1H3B": "H3C2", "MLL2": "KMT2D", "MLL3": "KMT2C", "FAM123B": "AMER1",
           "C11orf30": "EMSY", "MRE11A": "MRE11", "PARK2": "PRKN", "RAD21B": "RAD21L1", "WHSC1": "NSD2"}


def mane_exons(symbol):
    symbol = ALIASES.get(symbol, symbol)
    """Exons of the MANE Select transcript (falling back to the canonical one) for a gene symbol."""
    g = get(f"{ENSEMBL}/lookup/symbol/homo_sapiens/{symbol}?expand=1;content-type=application/json")
    if not g or "Transcript" not in g:
        return None, []
    txs = g["Transcript"]
    pick = None
    for t in txs:
        if t.get("is_canonical") and t.get("biotype") == "protein_coding":
            pick = t
    for t in txs:                                   # MANE Select wins over canonical when both are present
        if "MANE_Select" in json.dumps(t.get("Parent", "")) or t.get("mane") == "MANE_Select":
            pick = t
    pick = pick or (txs[0] if txs else None)
    if not pick:
        return None, []
    exons = [(e["start"], e["end"]) for e in pick.get("Exon", [])]
    return {"gene": symbol, "chrom": "chr" + str(g["seq_region_name"]), "transcript": pick["id"],
            "strand": g.get("strand"), "start": g["start"], "end": g["end"]}, sorted(exons)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", required=True, help="one gene symbol per line (# comments allowed)")
    ap.add_argument("--pad", type=int, default=20, help="bases of intronic padding on each side of an exon")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    symbols = [l.strip() for l in open(a.genes) if l.strip() and not l.startswith("#")]
    print(f"resolving {len(symbols)} genes via Ensembl", file=sys.stderr)
    bed, rows, missing = [], [], []
    for i, s in enumerate(symbols, start=1):
        info, exons = mane_exons(s)
        if not info or not exons:
            missing.append(s)
            print(f"  ! {s}: not resolved", file=sys.stderr)
            continue
        span = 0
        for st, en in exons:
            bed.append((info["chrom"], max(0, st - 1 - a.pad), en + a.pad, s))
            span += (en - st + 1) + 2 * a.pad
        rows.append({"gene": s, "chrom": info["chrom"], "transcript": info["transcript"],
                     "n_exons": len(exons), "gene_start": info["start"], "gene_end": info["end"],
                     "exonic_bp_padded": span})
        if i % 20 == 0:
            print(f"  {i}/{len(symbols)}", file=sys.stderr)
        time.sleep(0.08)

    def sort_key(b):
        c = b[0][3:]
        return (int(c) if c.isdigit() else 100 + ord(c[0]), b[1])
    bed.sort(key=sort_key)
    with open(a.out + ".bed", "w") as fh:
        for c, s, e, n in bed:
            fh.write(f"{c}\t{s}\t{e}\t{n}\n")
    write_tsv(rows, a.out + "_genes.tsv",
              ["gene", "chrom", "transcript", "n_exons", "exonic_bp_padded", "gene_start", "gene_end"])
    total = sum(r["exonic_bp_padded"] for r in rows)
    write_json({"genes_requested": len(symbols), "genes_resolved": len(rows), "genes_unresolved": missing,
                "intervals": len(bed), "approx_target_mb": round(total / 1e6, 3), "padding_bp": a.pad,
                "method": "Exons of the MANE Select (or canonical protein-coding) transcript from Ensembl REST, "
                          "padded to include splice sites.",
                "caveat": "This approximates the assay's design: the manufacturer targets selected exons of some "
                          "genes rather than all of them, so the true territory is smaller. Use the observed "
                          "coverage (mosdepth) for the TMB denominator, and this BED for calling intervals."},
               a.out + "_summary.json")
    print(f"{len(rows)}/{len(symbols)} genes -> {len(bed)} intervals, ~{total/1e6:.2f} Mb ({a.out}.bed)",
          file=sys.stderr)
    if missing:
        print("unresolved: " + ", ".join(missing), file=sys.stderr)


if __name__ == "__main__":
    main()
