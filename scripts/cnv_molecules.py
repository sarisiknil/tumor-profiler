#!/usr/bin/env python3
"""Coverage-based copy-number hints for a targeted panel, tumour-only, without a normal reference cohort.

Rationale: in anchored multiplex PCR the number of *unique molecules* (UMI families) captured per target is
roughly proportional to the number of template copies, so amplifications show up as targets with far more
molecules than the panel median (ArcherDX technical note PN-MKT-0032 demonstrates PDGFRA/KIT amplification this
way at only 2 M reads). We compute a robust per-gene log2 ratio versus the panel median, using mosdepth
per-region depth (or a per-target molecule count) as input.

Hard limitations, printed with the output and to be quoted in the report:
  * no matched normal and no pooled panel-of-normals -> capture/GC/primer-efficiency bias is NOT corrected;
    a gene can look amplified merely because its primers work better
  * tumour purity and ploidy are unknown, so log2 ratios cannot be converted to absolute copy number
  * single-exon or focal events are unreliable with few targets per gene
  * CNVkit with a flat reference is the next step up, and a real cohort reference is the proper solution
"""
import argparse, statistics, sys, gzip
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

AMPLIFIED_ONCOGENES = {"EGFR","ERBB2","MET","MYC","MYCN","CCND1","CCNE1","FGFR1","FGFR2","KRAS","PDGFRA","KIT",
                       "CDK4","MDM2","AR","ALK","BRAF","TERT","AKT2","RICTOR"}
DELETED_TSG = {"CDKN2A","CDKN2B","PTEN","RB1","SMAD4","STK11","TP53","NF1","NF2","BRCA1","BRCA2","ATM","APC"}

def read_regions(path):
    """mosdepth *.regions.bed(.gz): chrom start end [name] mean_depth"""
    op = gzip.open if str(path).endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                depth = float(f[-1]); start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            name = f[3] if len(f) >= 5 else f"{f[0]}:{start}-{end}"
            out.append({"chrom": f[0], "start": start, "end": end, "gene": name.split("_")[0], "depth": depth})
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--regions", required=True, help="mosdepth *.regions.bed.gz over the panel targets")
    ap.add_argument("--min-depth", type=float, default=50, help="ignore targets below this depth")
    ap.add_argument("--gain-log2", type=float, default=1.0)
    ap.add_argument("--loss-log2", type=float, default=-1.0)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    regs = [r for r in read_regions(a.regions) if r["depth"] >= a.min_depth]
    if not regs:
        write_tsv([], a.out + ".tsv"); write_json({"error": "no targets above min depth"}, a.out + "_summary.json"); return
    med = statistics.median(r["depth"] for r in regs)
    by_gene = {}
    for r in regs:
        by_gene.setdefault(r["gene"], []).append(r["depth"])
    import math
    rows = []
    for g, ds in by_gene.items():
        gm = statistics.median(ds)
        log2 = math.log2(gm / med) if med > 0 and gm > 0 else 0.0
        call = "neutral"
        if log2 >= a.gain_log2:
            call = "gain (possible amplification)" if g in AMPLIFIED_ONCOGENES else "gain"
        elif log2 <= a.loss_log2:
            call = "loss (possible deletion)" if g in DELETED_TSG else "loss"
        rows.append({"gene": g, "n_targets": len(ds), "median_depth": round(gm, 1),
                     "panel_median_depth": round(med, 1), "log2_ratio": round(log2, 3), "call": call,
                     "known_amplified_oncogene": g in AMPLIFIED_ONCOGENES,
                     "known_deleted_tsg": g in DELETED_TSG})
    rows.sort(key=lambda r: -abs(r["log2_ratio"]))
    write_tsv(rows, a.out + ".tsv",
              ["gene","call","log2_ratio","median_depth","panel_median_depth","n_targets",
               "known_amplified_oncogene","known_deleted_tsg"])
    write_json({"panel_median_depth": round(med, 1), "n_targets_used": len(regs), "n_genes": len(rows),
                "gains": [r["gene"] for r in rows if r["call"].startswith("gain")],
                "losses": [r["gene"] for r in rows if r["call"].startswith("loss")],
                "thresholds": {"gain_log2": a.gain_log2, "loss_log2": a.loss_log2, "min_depth": a.min_depth},
                "limitations": ["No matched normal or panel of normals: primer-efficiency, GC and capture bias "
                                "are uncorrected, so a 'gain' may be a technical artefact.",
                                "Tumour purity and ploidy are unknown; log2 ratios cannot be converted to "
                                "absolute copy number.",
                                "Genes with few targets give noisy estimates; interpret only large, "
                                "consistent shifts.",
                                "For a defensible CNV call use CNVkit with a pooled reference of >= 10-20 "
                                "samples processed identically (Talevich et al. 2016)."]},
               a.out + "_summary.json")
    top = [r for r in rows if r["call"] != "neutral"][:10]
    print(f"panel median depth {med:.0f}; {len(top)} non-neutral genes: "
          + ", ".join(f"{r['gene']}({r['log2_ratio']:+.2f})" for r in top), file=sys.stderr)

if __name__ == "__main__":
    main()
