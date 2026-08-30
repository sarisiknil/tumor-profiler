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

# Paralogue families: reads and PCR products move between these genes easily, so a "fusion" joining two members
# of one family is far more likely to be a mapping or chemistry artefact than biology.
PARALOGUE_FAMILIES = [
    {"FGFR1", "FGFR2", "FGFR3", "FGFR4"},
    {"NTRK1", "NTRK2", "NTRK3"},
    {"ERBB2", "ERBB3", "ERBB4", "EGFR"},
    {"PIK3CA", "PIK3CB", "PIK3CD"},
    {"AKT1", "AKT2", "AKT3"},
    {"RAF1", "BRAF", "ARAF"},
    {"HRAS", "KRAS", "NRAS"},
]

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
    def same_family(g1, g2):
        return any({g1, g2} <= fam for fam in PARALOGUE_FAMILIES)

    raw = read_arriba(a.arriba)
    seen_pairs = {(f.get("gene1", "").split(",")[0], f.get("gene2", "").split(",")[0]) for f in raw}
    # a gene that turns up with many different partners is an artefact hub, not a driver
    partner_count = {}
    for f in raw:
        x, y = f.get("gene1", "").split(",")[0], f.get("gene2", "").split(",")[0]
        partner_count.setdefault(x, set()).add(y)
        partner_count.setdefault(y, set()).add(x)
    rows = []
    for f in raw:
        g1 = (f.get("gene1") or "").split(",")[0]
        g2 = (f.get("gene2") or "").split(",")[0]
        support = 0
        for k in ("split_reads1", "split_reads2", "discordant_mates"):
            try: support += int(f.get(k, 0) or 0)
            except ValueError: pass
        targetable = sorted({g for g in (g1, g2) if g in TARGETABLE})
        # ---- artefact screening, in the order that matters most for this chemistry -------------------
        risk = []
        if (g2, g1) in seen_pairs and g1 != g2:
            risk.append("reciprocal call: the same pair is reported in both orientations, which a real fusion "
                        "does not produce")
        if same_family(g1, g2):
            risk.append(f"{g1} and {g2} are paralogues: reads and PCR products cross between them easily")
        if panel and g1 in panel and g2 in panel:
            risk.append("both partners are panel targets: anchored multiplex PCR primes one gene and discovers "
                        "the partner, so a join between two primed genes is a classic PCR chimera")
        if support <= 2:
            risk.append(f"only {support} supporting read(s)")
        for g in (g1, g2):
            if len(partner_count.get(g, ())) >= 3:
                risk.append(f"{g} is reported with {len(partner_count[g])} different partners in this sample: "
                            "an artefact hub rather than a driver")
                break
        unnamed = [g for g in (g1, g2) if g.startswith("ENSG") or g in (".", "")]
        if unnamed:
            risk.append(f"partner {unnamed[0]} has no gene symbol: an uncharacterised or novel locus, "
                        "commonly a mapping artefact")
        if "5'-5'" in (f.get("type") or ""):
            risk.append("5'-5' orientation: the two 5' ends are joined, which cannot produce a functional "
                        "fusion protein")

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
                     "artefact_risk": "high" if len(risk) >= 2 else ("possible" if risk else "not flagged"),
                     "artefact_reasons": "; ".join(risk),
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
    risk_rank = {"not flagged": 0, "possible": 1, "high": 2}
    # artefact risk outranks everything: a flagged call must not sit at the top of the table just because one
    # partner happens to be a druggable kinase
    rows.sort(key=lambda r: (risk_rank.get(r["artefact_risk"], 0), 0 if r["targetable"] else 1,
                             conf_rank.get(r["confidence"], 3), -r["support_total"]))
    cols = ["fusion","artefact_risk","confidence","reading_frame","targetable","support_total","split_reads",
            "split_reads2","discordant_mates","artefact_reasons","breakpoint1","breakpoint2","site1","site2",
            "type","retained_domains","notes","gene1","gene2","coverage1","coverage2","filters"]
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
                "n_flagged_high_artefact_risk": sum(1 for r in rows if r["artefact_risk"] == "high"),
                "n_not_flagged": sum(1 for r in rows if r["artefact_risk"] == "not flagged"),
                "by_confidence": {c: sum(1 for r in rows if r["confidence"] == c) for c in {r["confidence"] for r in rows}},
                "targetable_fusions": [r["fusion"] for r in rows if r["targetable"]],
                "n_discarded_with_targetable_gene": len(rescued),
                "caveats": ["Arriba is tuned for whole-transcriptome data; on anchored-multiplex-PCR panels it "
                            "recovers ~86 % (lung) / 57 % (sarcoma) of the vendor pipeline's calls "
                            "(Capone et al. 2022), so low-confidence calls are kept and shown.",
                            "Arriba does not detect intragenic exon-skipping events such as MET exon 14 skipping "
                            "— use exon_skipping.py for those.",
                            "Any fusion used for a treatment decision needs orthogonal confirmation "
                            "(FISH, IHC or RT-PCR).",
                            "Artefact screening flags reciprocal calls, paralogue pairs, joins between two "
                            "primed panel genes, and calls with minimal read support. On anchored multiplex "
                            "PCR these patterns dominate: the assay amplifies its target genes in one tube, so "
                            "partially extended products readily cross-prime between them."]},
               a.out + "_summary.json")
    print(f"{len(rows)} fusions ({sum(1 for r in rows if r['targetable'])} with a targetable partner)", file=sys.stderr)

if __name__ == "__main__":
    main()
