#!/usr/bin/env python3
"""Turn an Arriba fusions.tsv into an interpretable, panel-aware fusion report.

What it adds on top of Arriba's own output:
  * keeps Arriba's confidence (high/medium/low) but re-ranks by clinical relevance: is either partner a
    known targetable kinase (ALK, ROS1, RET, NTRK1/2/3, MET, FGFR1-3, BRAF, NRG1, EGFR, PDGFRA/B, ABL1)?
  * reading frame, retained protein domains, breakpoint coordinates and supporting read counts
  * a per-fusion note on what would be needed to call it clinically (orthogonal confirmation)
  * flags fusions where only one partner is on the panel — expected for anchored multiplex PCR, which uses a
    single gene-specific primer and discovers the partner (Zheng et al. 2014); this is a feature, not an error.

Arriba's filters are tuned for whole-transcriptome data; on AMP/FusionPlex panels open-source callers recover
fewer fusions than the vendor pipeline (Capone et al. 2022: Arriba 86 % of Archer calls in lung, 57 % in sarcoma),
so low-confidence calls are reported, not dropped.
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

TARGETABLE = {"ALK","ROS1","RET","NTRK1","NTRK2","NTRK3","MET","FGFR1","FGFR2","FGFR3","BRAF","NRG1","EGFR",
              "PDGFRA","PDGFRB","ABL1","ERBB2","AKT3","ESR1","THADA","FGR","MAML2","NUTM1","PRKACA","RAF1"}
RECURRENT_PARTNERS = {"EML4","KIF5B","CD74","SLC34A2","EZR","TPM3","SDC4","KLC1","STRN","CCDC6","NCOA4",
                      "TFG","ETV6","LMNA","BCAN","SQSTM1","PAX8","TRIM24","GOPC","SLC45A3","TMPRSS2"}

def read_arriba(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arriba", required=True, help="Arriba fusions.tsv")
    ap.add_argument("--discarded", help="Arriba fusions.discarded.tsv (optional; scanned for targetable genes)")
    ap.add_argument("--panel-genes", help="one gene symbol per line: genes targeted by the RNA panel")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    panel = set()
    if a.panel_genes and Path(a.panel_genes).exists():
        panel = {l.strip().split("\t")[0] for l in open(a.panel_genes) if l.strip()}
    rows = []
    for f in read_arriba(a.arriba):
        g1 = (f.get("gene1") or "").split(",")[0]
        g2 = (f.get("gene2") or "").split(",")[0]
        support = 0
        for k in ("split_reads1", "split_reads2", "discordant_mates"):
            try: support += int(f.get(k, 0) or 0)
            except ValueError: pass
        targetable = sorted({g for g in (g1, g2) if g in TARGETABLE})
        notes = []
        if targetable:
            notes.append(f"targetable kinase partner: {','.join(targetable)}")
        if (g1 in RECURRENT_PARTNERS) or (g2 in RECURRENT_PARTNERS):
            notes.append("recurrent fusion partner")
        if panel:
            on = [g for g in (g1, g2) if g in panel]
            if len(on) == 1:
                notes.append(f"only {on[0]} is a panel target — expected for anchored multiplex PCR "
                             "(partner discovered without prior knowledge)")
            elif not on:
                notes.append("neither partner is a panel target — treat with caution (possible artefact)")
        if (f.get("reading_frame") or "") == "in-frame":
            notes.append("in-frame")
        rows.append({"fusion": f"{g1}--{g2}", "gene1": g1, "gene2": g2,
                     "confidence": f.get("confidence", ""), "reading_frame": f.get("reading_frame", ""),
                     "type": f.get("type", ""), "site1": f.get("site1", ""), "site2": f.get("site2", ""),
                     "breakpoint1": f.get("breakpoint1", ""), "breakpoint2": f.get("breakpoint2", ""),
                     "split_reads": f.get("split_reads1", ""), "split_reads2": f.get("split_reads2", ""),
                     "discordant_mates": f.get("discordant_mates", ""), "support_total": support,
                     "coverage1": f.get("coverage1", ""), "coverage2": f.get("coverage2", ""),
                     "retained_domains": f.get("retained_protein_domains", "")[:200],
                     "targetable": ",".join(targetable), "notes": "; ".join(notes),
                     "filters": f.get("filters", "")})
    conf_rank = {"high": 0, "medium": 1, "low": 2, "": 3}
    rows.sort(key=lambda r: (0 if r["targetable"] else 1, conf_rank.get(r["confidence"], 3), -r["support_total"]))
    cols = ["fusion","confidence","reading_frame","targetable","support_total","split_reads","split_reads2",
            "discordant_mates","breakpoint1","breakpoint2","site1","site2","type","retained_domains","notes",
            "gene1","gene2","coverage1","coverage2","filters"]
    write_tsv(rows, a.out + ".tsv", cols)
    rescued = []
    if a.discarded and Path(a.discarded).exists():
        for f in read_arriba(a.discarded):
            g1 = (f.get("gene1") or "").split(",")[0]; g2 = (f.get("gene2") or "").split(",")[0]
            if g1 in TARGETABLE or g2 in TARGETABLE:
                rescued.append({"fusion": f"{g1}--{g2}", "filters": f.get("filters", ""),
                                "split_reads": f.get("split_reads1", ""), "confidence": f.get("confidence", "")})
        write_tsv(rescued[:100], a.out + "_discarded_targetable.tsv", ["fusion","confidence","split_reads","filters"])
    write_json({"n_fusions": len(rows),
                "by_confidence": {c: sum(1 for r in rows if r["confidence"] == c) for c in {r["confidence"] for r in rows}},
                "targetable_fusions": [r["fusion"] for r in rows if r["targetable"]],
                "n_discarded_with_targetable_gene": len(rescued),
                "caveats": ["Arriba is tuned for whole-transcriptome data; on anchored-multiplex-PCR panels it "
                            "recovers ~86 % (lung) / 57 % (sarcoma) of the vendor pipeline's calls "
                            "(Capone et al. 2022), so low-confidence calls are kept and shown.",
                            "Arriba does not detect intragenic exon-skipping events such as MET exon 14 skipping "
                            "— use exon_skipping.py for those.",
                            "Any fusion used for a treatment decision needs orthogonal confirmation "
                            "(FISH, IHC or RT-PCR)."]},
               a.out + "_summary.json")
    print(f"{len(rows)} fusions ({sum(1 for r in rows if r['targetable'])} with a targetable partner)", file=sys.stderr)

if __name__ == "__main__":
    main()
