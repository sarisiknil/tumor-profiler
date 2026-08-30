#!/usr/bin/env python3
"""Map the sample's alterations onto canonical oncogenic signalling pathways.

Uses the ten pathways curated by Sanchez-Vega et al. 2018 (Cell 173:321) plus two supplementary sets
(DNA-damage repair/MMR and chromatin regulators) that the original ten deliberately exclude.

IMPORTANT (and stated in the output): this is a *membership* summary — which pathways contain an altered gene —
not a statistical enrichment. Enrichment testing is meaningless for a targeted panel because the gene universe
is chosen by the assay, and for a single sample there is no null distribution.
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

def read_gmt(p):
    sets = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            sets[f[0]] = {"description": f[1], "genes": [g for g in f[2:] if g]}
    return sets

def read_tsv(p):
    if not Path(p).exists():
        return []
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True, help="tiered or filtered variants TSV")
    ap.add_argument("--fusions", help="fusion summary TSV (columns gene1, gene2)")
    ap.add_argument("--gmt", default=str(Path(__file__).resolve().parents[1] / "resources" / "oncogenic_pathways.gmt"))
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    sets = read_gmt(a.gmt)
    alt = {}
    for v in read_tsv(a.variants):
        if v.get("class") and v["class"] != "SOMATIC_LIKELY":
            continue
        g = v.get("gene")
        if g:
            alt.setdefault(g, []).append(f"{v.get('hgvsp') or v.get('consequence','variant')} "
                                         f"(VAF {v.get('vaf','?')}, tier {v.get('amp_tier','-')})")
    if a.fusions and Path(a.fusions).exists():
        for f in read_tsv(a.fusions):
            for g in (f.get("gene1"), f.get("gene2")):
                if g:
                    alt.setdefault(g, []).append(f"fusion {f.get('gene1')}-{f.get('gene2')} "
                                                 f"({f.get('confidence','?')})")
    rows = []
    for name, s in sets.items():
        hits = [g for g in s["genes"] if g in alt]
        rows.append({"pathway": name, "n_genes_in_set": len(s["genes"]), "n_altered": len(hits),
                     "altered_genes": ",".join(hits),
                     "alterations": " | ".join(f"{g}: {'; '.join(alt[g])}" for g in hits),
                     "source": s["description"]})
    rows.sort(key=lambda r: -r["n_altered"])
    write_tsv(rows, a.out + ".tsv", ["pathway","n_altered","n_genes_in_set","altered_genes","alterations","source"])
    write_json({"pathways_hit": [r["pathway"] for r in rows if r["n_altered"]],
                "genes_altered": sorted(alt),
                "unassigned_genes": sorted(g for g in alt if not any(g in s["genes"] for s in sets.values())),
                "method": "Membership mapping against Sanchez-Vega et al. 2018 (Cell) curated pathway gene lists "
                          "plus supplementary DNA-repair and chromatin sets.",
                "caveat": "This is not an enrichment analysis: the gene universe is defined by the panel, and a "
                          "single sample provides no null distribution. Absence of a pathway hit may simply mean "
                          "the panel does not cover that pathway's genes."},
               a.out + "_summary.json")
    print("\n".join(f"{r['pathway']}: {r['n_altered']} altered ({r['altered_genes']})" for r in rows if r["n_altered"]),
          file=sys.stderr)

if __name__ == "__main__":
    main()
