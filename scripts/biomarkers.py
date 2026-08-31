#!/usr/bin/env python3
"""Panel-level biomarkers from tumour-only targeted DNA data: TMB, mutation spectrum, FFPE damage indicator.

TMB: nonsynonymous somatic variants per Mb of *assessable* territory (regions with depth >= config
min_depth_assessable, taken from mosdepth). Panel TMB is only comparable across assays after harmonisation
(Merino et al. 2020): small panels give wide confidence intervals and, in tumour-only mode, residual germline
calls inflate the value. We therefore report the count, the denominator, a 95 % Poisson interval and the caveats.

FFPE indicator: fraction of C>T/G>A among somatic SNVs (Do & Dobrovic 2015; Diossy et al. 2021) — a high value
concentrated at low VAF suggests formalin deamination rather than biology.

MSI and mutational signatures are deliberately NOT computed; the reason is recorded in the output so the report
can state it explicitly.
"""
import argparse, json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_json, load_config

NONSYN = ("missense_variant", "stop_gained", "stop_lost", "start_lost", "frameshift_variant",
          "inframe_insertion", "inframe_deletion", "splice_acceptor_variant", "splice_donor_variant",
          "protein_altering_variant", "coding_sequence_variant")

def read_tsv(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

def assessable_mb(mosdepth_regions, min_depth):
    """Sum the length of regions whose mean depth >= min_depth, from a mosdepth *.regions.bed(.gz)."""
    import gzip
    total = 0
    op = gzip.open if str(mosdepth_regions).endswith(".gz") else open
    with op(mosdepth_regions, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                depth = float(f[-1]); start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            if depth >= min_depth:
                total += end - start
    return total / 1e6

def _norm_inv(p):
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02,
         -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00,
         4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
    pl, pu = 0.02425, 1 - 0.02425
    if p < pl:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > pu:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5; r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)

def _chi2_inv(p, df):
    if df <= 0:
        return 0.0
    z = _norm_inv(p)
    return df * (1 - 2 / (9 * df) + z * math.sqrt(2 / (9 * df))) ** 3

def poisson_ci(n, exposure):
    """95 % confidence interval for the rate n/exposure (Garwood/chi-square)."""
    if exposure <= 0:
        return (None, None)
    lo = 0.0 if n == 0 else 0.5 * _chi2_inv(0.025, 2 * n)
    hi = 0.5 * _chi2_inv(0.975, 2 * (n + 1))
    return (round(lo / exposure, 2), round(hi / exposure, 2))

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True, help="filtered variants TSV")
    ap.add_argument("--mosdepth-regions", help="mosdepth *.regions.bed.gz for the assessable denominator")
    ap.add_argument("--panel-size-mb", type=float, help="fallback panel size in Mb")
    ap.add_argument("--config", default=None)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    cfg = load_config(a.config)
    min_depth = cfg["panel"]["min_depth_assessable"]
    rows = read_tsv(a.variants)
    som_all = [r for r in rows if r.get("class") == "SOMATIC_LIKELY"]
    # TMB must be computed on the calls one actually believes. Variants carrying an artefact flag - positional
    # clustering, or membership of a dominant substitution class - are excluded from the numerator, and both
    # figures are reported so the effect of the screen is visible rather than assumed.
    som = [r for r in som_all if r.get("artefact_risk", "not flagged") == "not flagged"]
    nonsyn = [r for r in som if any(c in (r.get("consequence") or "") for c in NONSYN)]
    nonsyn_all = [r for r in som_all if any(c in (r.get("consequence") or "") for c in NONSYN)]
    if a.mosdepth_regions and Path(a.mosdepth_regions).exists():
        mb = assessable_mb(a.mosdepth_regions, min_depth)
        denom_src = f"mosdepth regions with mean depth >= {min_depth}x"
    else:
        mb = a.panel_size_mb or 0.0
        denom_src = "panel size given on the command line (mosdepth output not available)"
    tmb = round(len(nonsyn) / mb, 2) if mb else None
    tmb_unscreened = round(len(nonsyn_all) / mb, 2) if mb else None
    lo, hi = poisson_ci(len(nonsyn), mb) if mb else (None, None)
    snvs = [r for r in som_all if len(r.get("ref", "")) == 1 and len(r.get("alt", "")) == 1]
    fold = {("G", "A"): ("C", "T"), ("T", "C"): ("A", "G"), ("G", "T"): ("C", "A"),
            ("G", "C"): ("C", "G"), ("A", "T"): ("T", "A"), ("A", "G"): ("T", "C"), ("A", "C"): ("T", "G")}
    spectrum = {}
    for r in snvs:
        pair = fold.get((r["ref"], r["alt"]), (r["ref"], r["alt"]))
        k = f"{pair[0]}>{pair[1]}"
        spectrum[k] = spectrum.get(k, 0) + 1
    ct = spectrum.get("C>T", 0)
    low_vaf_ct = sum(1 for r in snvs if (r["ref"], r["alt"]) in (("C", "T"), ("G", "A"))
                     and float(r.get("vaf") or 0) < 0.10)
    out = {
        "tmb": {"nonsynonymous_somatic_variants": len(nonsyn),
                "nonsynonymous_before_artefact_screening": len(nonsyn_all),
                "tmb_per_mb_before_screening": tmb_unscreened,
                "assessable_mb": round(mb, 3),
                "denominator_source": denom_src, "tmb_per_mb": tmb, "ci95_per_mb": [lo, hi],
                "reference_threshold": "TMB-high is commonly defined as >= 10 mutations/Mb (FDA pembrolizumab "
                                       "tissue-agnostic indication, measured with FoundationOne CDx)",
                "meaningful": mb >= 0.1,
                "screening_note": (f"{len(som_all) - len(som)} of {len(som_all)} somatic candidates carried an "
                                   "artefact flag and were excluded from the numerator."),
                "caveats": ([] if mb >= 0.1 else
                            [f"Assessable territory is only {mb:.3f} Mb: far below the ~1 Mb that panel TMB "
                             "estimation assumes. The value is reported for completeness but must not be "
                             "interpreted as a tumour mutational burden."]) +
                           ["Panel TMB is not directly comparable between assays without harmonisation "
                            "(Merino et al. 2020, J Immunother Cancer 8:e000147).",
                            "Tumour-only calling leaves residual germline variants, which inflate TMB.",
                            "With a small panel the Poisson confidence interval is wide; always report it."]},
        "mutation_spectrum": spectrum,
        "ffpe_indicator": {"C>T_or_G>A_fraction_of_snvs": round(ct / len(snvs), 3) if len(snvs) >= 10 else None,
                           "n_snvs": len(snvs),
                           "interpretable": len(snvs) >= 10,
                           "low_vaf_C>T_count": low_vaf_ct,
                           "note": ("" if len(snvs) >= 10 else
                                    f"Only {len(snvs)} somatic SNV(s): the C>T fraction is not interpretable "
                                    "and is therefore not reported. ") +
                                   "Formalin fixation deaminates cytosine, producing C>T/G>A changes at low VAF "
                                   "(Do & Dobrovic 2015). A high fraction argues for artefact filtering "
                                   "(Mutect2 read-orientation model, SOBDetector), not for a biological signature."},
        "not_computed": {"MSI": "MSIsensor2/-pro require a panel-specific site model (or >= 20 matched normals) "
                                "and >= 50 covered microsatellite loci; not satisfied by this panel.",
                         "mutational_signatures": f"SBS refitting needs >= 100-200 SNVs (PCGR default 200); only "
                                                  f"{len(snvs)} somatic SNVs were called."},
        "counts": {"all_variants": len(rows), "somatic_likely": len(som), "somatic_snvs": len(snvs),
                   "nonsynonymous": len(nonsyn)},
    }
    write_json(out, a.out)
    print(json.dumps(out["tmb"], indent=1), file=sys.stderr)

if __name__ == "__main__":
    main()
