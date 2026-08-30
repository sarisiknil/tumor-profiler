#!/usr/bin/env python3
"""Assign clinical-significance tiers to somatic alterations.

Two complementary systems, both implemented as transparent, rule-based mappings from the evidence table:

AMP/ASCO/CAP (Li et al. 2017, J Mol Diagn 19:4-23)
  Tier I   strong clinical significance   – FDA-approved therapy / professional guidelines for this tumour type
                                            (OncoKB level 1-2, CIViC A-B predictive evidence in this disease)
  Tier II  potential clinical significance – therapy approved in another tumour type, or clinical-trial evidence
                                            (OncoKB 3A-4/R1, CIViC C-D)
  Tier III variant of unknown significance – in a cancer gene, no clinical evidence
  Tier IV  benign / likely benign
ESCAT (Mateo et al. 2018, Ann Oncol 29:1895)
  I-A/I-B ready for routine use, II investigational with clinical data, III benefit in another tumour type,
  IV preclinical, V evidence of actionability but no clinical benefit, X no evidence.

Mapping tables are printed into the output so a reader can audit exactly why a variant got its tier;
tiers are a decision aid for education, NOT a clinical report.
"""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

# Evidence in the PATIENT'S tumour type.  AMP/ASCO/CAP is explicit that Tier I means an approved therapy or
# guideline *for this tumour type*; the same evidence in a different tumour type is Tier II (level C).
ONCOKB_TIER = {"LEVEL_1": ("I", "I-A"), "LEVEL_2": ("I", "I-B"), "LEVEL_3A": ("II", "II"),
               "LEVEL_3B": ("II", "III"), "LEVEL_4": ("II", "IV"), "LEVEL_R1": ("I", "I-A"),
               "LEVEL_R2": ("II", "II")}
CIVIC_TIER = {"A": ("I", "I-B"), "B": ("I", "I-B"), "C": ("II", "II"), "D": ("II", "IV"), "E": ("III", "V")}
# The same evidence when the disease does not match: strong evidence elsewhere is Tier II, weak stays where it is.
OTHER_DISEASE = {"I": ("II", "III"), "II": ("II", "IV"), "III": ("III", "V"), "IV": ("IV", "X")}
BENIGN = ("benign", "likely_benign")
PATHOGENIC = ("pathogenic", "likely_pathogenic")
HIGH_IMPACT = ("stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost")

def read_tsv(p):
    if not Path(p).exists():
        return []
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

def rank(tier):
    return {"I": 0, "II": 1, "III": 2, "IV": 3}.get(tier, 4)


def disease_matches(evidence_disease: str, cancer_type: str) -> bool:
    """Loose match between an evidence item's disease and the patient's diagnosis.

    Deliberately generous (substring / shared significant word) rather than an ontology lookup: the aim is to
    separate 'evidence in this tumour type' from 'evidence in a different tumour type', and a missed match
    only ever makes the tier more conservative."""
    if not cancer_type:
        return True                     # no diagnosis configured: do not pretend to judge disease context
    d, c = (evidence_disease or "").lower(), cancer_type.lower()
    if not d:
        return False
    if d in c or c in d:
        return True
    stop = {"cancer", "carcinoma", "tumor", "tumour", "neoplasm", "malignant", "advanced", "metastatic", "the"}
    dw = {w for w in d.replace(",", " ").split() if len(w) > 4 and w not in stop}
    cw = {w for w in c.replace(",", " ").split() if len(w) > 4 and w not in stop}
    return bool(dw & cw)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--evidence", required=True, help="evidence TSV from annotate_evidence.py")
    ap.add_argument("--cancer-genes", default=str(Path(__file__).resolve().parents[1] / "resources" / "cancer_genes.txt"))
    ap.add_argument("--cancer-type", default=None,
                    help="the patient's diagnosis; evidence from other tumour types is capped at Tier II")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    if a.cancer_type is None:
        try:
            import yaml
            a.cancer_type = yaml.safe_load(
                open(Path(__file__).resolve().parents[1] / "config" / "config.yaml")).get("cancer_type", "")
        except Exception:
            a.cancer_type = ""
    ev_by_key = {}
    for e in read_tsv(a.evidence):
        ev_by_key.setdefault(e["key"], []).append(e)
    cancer_genes = set()
    if Path(a.cancer_genes).exists():
        cancer_genes = {l.strip().split("\t")[0] for l in open(a.cancer_genes) if l.strip() and not l.startswith("#")}
    rows = []
    for v in read_tsv(a.variants):
        if v.get("class") != "SOMATIC_LIKELY":
            continue
        evs = ev_by_key.get(v["key"], [])
        best_amp, best_escat, why = "III", "X", []
        clin = (v.get("clinvar_sig") or "").lower()
        matched_disease = False
        for e in evs:
            lvl = (e.get("level") or "").upper()
            src = e.get("source", "")
            if src == "OncoKB" and lvl in ONCOKB_TIER:
                amp, esc = ONCOKB_TIER[lvl]
                note = f"OncoKB {lvl} ({e.get('therapies','')})"
            elif src == "CIViC" and lvl and lvl[0] in CIVIC_TIER and (e.get("type") or "").lower().startswith("predict"):
                amp, esc = CIVIC_TIER[lvl[0]]
                note = f"CIViC {lvl} {e.get('significance','')} ({e.get('therapies','')})"
            else:
                continue
            same = disease_matches(e.get("disease", ""), a.cancer_type)
            matched_disease |= same
            if not same:
                amp, esc = OTHER_DISEASE.get(amp, (amp, esc))
                note += f" — in {e.get('disease','another tumour type')}, not {a.cancer_type}"
            else:
                note += f" — in {e.get('disease','')}"
            why.append(note)
            if rank(amp) < rank(best_amp):
                best_amp, best_escat = amp, esc
        if not evs and v.get("gene") in cancer_genes:
            why.append("cancer-gene, no clinical evidence found")
        if any(b in clin for b in BENIGN) and not any(p in clin for p in PATHOGENIC):
            best_amp, best_escat = "IV", "X"; why.append("ClinVar benign/likely benign")
        elif best_amp == "III" and (any(h in (v.get("consequence") or "") for h in HIGH_IMPACT)
                                    and v.get("gene") in cancer_genes):
            why.append("truncating variant in a cancer gene (potential loss of function)")
        rows.append({**v, "amp_tier": best_amp, "escat": best_escat,
                     "disease_matched_evidence": "yes" if matched_disease else "no",
                     "tier_rationale": "; ".join(why[:6]) or "no evidence found",
                     "n_evidence_items": len(evs),
                     "therapies": ";".join(sorted({e.get("therapies", "") for e in evs if e.get("therapies")}))})
    rows.sort(key=lambda r: (rank(r["amp_tier"]), -float(r.get("vaf") or 0)))
    cols = ["key","gene","hgvsp","hgvsc","consequence","amp_tier","escat","disease_matched_evidence",
            "therapies","tier_rationale","vaf","depth","n_callers","gnomad_af_max","clinvar_sig",
            "n_evidence_items"]
    write_tsv(rows, a.out + ".tsv", cols)
    counts = {}
    for r in rows: counts[r["amp_tier"]] = counts.get(r["amp_tier"], 0) + 1
    write_json({"tier_counts": counts, "n_somatic": len(rows), "cancer_type": a.cancer_type,
                "variants_with_evidence_in_this_tumour_type":
                    sum(1 for r in rows if r["disease_matched_evidence"] == "yes"),
                "systems": {"AMP/ASCO/CAP": "Li MM et al. 2017 J Mol Diagn 19(1):4-23",
                            "ESCAT": "Mateo J et al. 2018 Ann Oncol 29(9):1895-1902"},
                "mapping_used": {"OncoKB": ONCOKB_TIER, "CIViC": CIVIC_TIER},
                "disease_context_rule": "AMP/ASCO/CAP Tier I requires evidence in the patient's own tumour "
                                        "type; identical evidence from a different tumour type is capped at "
                                        "Tier II. Free knowledge bases cover tumour types unevenly, so a "
                                        "variant can be genuinely Tier I under a guideline (e.g. NCCN) while "
                                        "the open evidence available here only supports Tier II.",
                "disclaimer": "Educational output. Tier assignment here is automated and does not replace "
                              "expert review or a clinical molecular tumour board."},
               a.out + "_summary.json")
    print(f"tiers: {counts}", file=sys.stderr)

if __name__ == "__main__":
    main()
