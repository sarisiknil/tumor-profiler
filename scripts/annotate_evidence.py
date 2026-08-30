#!/usr/bin/env python3
"""Attach clinical evidence to somatic variants from free knowledge bases, and build a therapy/trial table.

Sources (all free for academic use; see docs/licensing.md):
  * CIViC   – public GraphQL API, data CC0, may be redistributed          https://civicdb.org/api/graphql
  * OncoKB  – free academic API with a personal token (env ONCOKB_TOKEN); annotations may NOT be redistributed,
              so the pipeline stores them only under results/ (git-ignored)
  * DGIdb   – drug–gene interactions, GraphQL                             https://dgidb.org/api/graphql
  * ClinicalTrials.gov API v2 – recruiting trials matching gene/drug terms
COSMIC is deliberately NOT queried: its licence forbids redistribution and excludes for-profit settings.

Usage: annotate_evidence.py --variants results/dna/variants_filtered.tsv --cancer-type "Lung Adenocarcinoma"
                            -o results/dna/evidence [--no-oncokb] [--no-trials]
"""
import argparse, json, os, sys, time, urllib.request, urllib.error
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

UA = "tumor-profiler/1.0 (educational pipeline)"

AA3TO1 = {"Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q","Glu":"E","Gly":"G","His":"H","Ile":"I",
          "Leu":"L","Lys":"K","Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W","Tyr":"Y","Val":"V",
          "Ter":"*", "Sec":"U", "Xaa":"X"}

def short_protein_change(hgvsp: str) -> str:
    """`ENSP00000288602.6:p.Val600Glu` -> `V600E`; knowledge bases index variants in one-letter notation."""
    if not hgvsp:
        return ""
    p = hgvsp.split(":")[-1]
    if p.startswith("p."):
        p = p[2:]
    for three, one in AA3TO1.items():
        p = p.replace(three, one)
    return p

def http(url, data=None, headers=None, tries=3):
    hdr = {"User-Agent": UA, "Accept": "application/json"}
    if data is not None:
        hdr["Content-Type"] = "application/json"
        data = json.dumps(data).encode()
    hdr.update(headers or {})
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, data=data, headers=hdr), timeout=90) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503) and i < tries - 1:
                time.sleep(3 * (i + 1)); continue
            return {"_error": f"HTTP {e.code} {e.reason}", "_url": url}
        except Exception as e:
            if i < tries - 1:
                time.sleep(2 * (i + 1)); continue
            return {"_error": str(e), "_url": url}

# ---------------------------------------------------------------- CIViC
# CIViC's GraphQL schema attaches evidence to variants through molecular profiles; the `evidenceItems`
# root field accepts a variantId directly, which is the shortest reliable path (schema checked 2026-08).
CIVIC_VARIANT_Q = """
query($name: String, $feature: String) {
  browseVariants(variantName: $name, featureName: $feature, first: 8) {
    nodes { id name featureName evidenceItemCount }
  }
}"""
CIVIC_EVIDENCE_Q = """
query($vid: Int!) {
  evidenceItems(variantId: $vid, first: 25, status: ACCEPTED) {
    nodes { id evidenceType evidenceLevel evidenceDirection significance
            therapies { name } disease { name } link description }
  }
}"""

def civic_for(gene, protein_change):
    """Look up a variant by gene + one-letter protein change, then pull its accepted evidence items."""
    if not gene or not protein_change:
        return []
    r = http("https://civicdb.org/api/graphql",
             {"query": CIVIC_VARIANT_Q, "variables": {"name": protein_change, "feature": gene}})
    nodes = (((r or {}).get("data") or {}).get("browseVariants") or {}).get("nodes") or []
    if not nodes:
        return []
    exact = [n for n in nodes if n["name"].upper() == protein_change.upper()] or nodes[:1]
    node = exact[0]
    ev = http("https://civicdb.org/api/graphql",
              {"query": CIVIC_EVIDENCE_Q, "variables": {"vid": int(node["id"])}})
    items = (((ev or {}).get("data") or {}).get("evidenceItems") or {}).get("nodes") or []
    out = []
    for it in items:
        out.append({"source": "CIViC", "variant": node["name"], "type": it.get("evidenceType"),
                    "level": it.get("evidenceLevel"), "significance": it.get("significance"),
                    "direction": it.get("evidenceDirection"),
                    "therapies": ",".join(t["name"] for t in (it.get("therapies") or [])),
                    "disease": (it.get("disease") or {}).get("name", ""),
                    "url": it.get("link") or f"https://civicdb.org/evidence/{it.get('id')}",
                    "description": (it.get("description") or "")[:300]})
    return out

# ---------------------------------------------------------------- OncoKB
def oncokb_for(gene, protein_change, cancer_type, token):
    if not token or not gene or not protein_change:
        return []
    alt = protein_change.replace("p.", "")
    url = ("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange"
           f"?hugoSymbol={urllib.parse.quote(gene)}&alteration={urllib.parse.quote(alt)}"
           f"&tumorType={urllib.parse.quote(cancer_type or '')}&referenceGenome=GRCh38")
    r = http(url, headers={"Authorization": f"Bearer {token}"})
    if not isinstance(r, dict) or r.get("_error"):
        return [{"source": "OncoKB", "error": (r or {}).get("_error", "no response")}]
    rows = []
    for t in (r.get("treatments") or []):
        rows.append({"source": "OncoKB", "level": t.get("level", ""),
                     "therapies": ", ".join(d.get("drugName", "") for d in (t.get("drugs") or [])),
                     "disease": ((t.get("levelAssociatedCancerType") or {}).get("mainType") or {}).get("name", ""),
                     "significance": r.get("oncogenic", ""), "type": "Predictive",
                     "description": (t.get("description") or "")[:300],
                     "url": "https://www.oncokb.org/gene/" + gene})
    if not rows:
        rows.append({"source": "OncoKB", "significance": r.get("oncogenic", "Unknown"),
                     "type": "Oncogenicity", "level": "", "therapies": "", "disease": "",
                     "description": (r.get("mutationEffect") or {}).get("description", "")[:300],
                     "url": "https://www.oncokb.org/gene/" + gene})
    return rows

# ---------------------------------------------------------------- DGIdb
DGIDB_Q = """
query($names: [String!]) {
  genes(names: $names) {
    nodes { name interactions { drug { name conceptId } interactionScore interactionTypes { type } sources { sourceDbName } } }
  }
}"""

def dgidb_for(genes):
    if not genes:
        return {}
    r = http("https://dgidb.org/api/graphql", {"query": DGIDB_Q, "variables": {"names": sorted(genes)}})
    out = {}
    for n in ((((r or {}).get("data") or {}).get("genes") or {}).get("nodes") or []):
        rows = []
        for it in (n.get("interactions") or [])[:15]:
            rows.append({"drug": (it.get("drug") or {}).get("name", ""),
                         "score": it.get("interactionScore"),
                         "types": ",".join(t.get("type", "") for t in (it.get("interactionTypes") or [])),
                         "sources": ",".join(s.get("sourceDbName", "") for s in (it.get("sources") or []))})
        out[n["name"]] = sorted(rows, key=lambda x: -(x["score"] or 0))
    return out

# ---------------------------------------------------------------- ClinicalTrials.gov v2
def trials_for(term, cancer_type, max_n=8):
    """Recruiting trials mentioning the gene. ClinicalTrials.gov's relevance ranking is loose, so results whose
    title/conditions/interventions do not actually mention the gene are dropped."""
    q = f"{term} {cancer_type}".strip()
    url = ("https://clinicaltrials.gov/api/v2/studies?query.term=" + urllib.parse.quote(q) +
           "&filter.overallStatus=RECRUITING&pageSize=" + str(max_n) +
           "&fields=NCTId,BriefTitle,OverallStatus,Phase,Condition,LocationCountry")
    r = http(url)
    rows = []
    for s in (r or {}).get("studies", [])[: max_n * 4]:
        p = s.get("protocolSection", {})
        blob = json.dumps(p).upper()
        if term.upper() not in blob:          # keep only trials that really mention the gene
            continue
        if len(rows) >= max_n:
            break
        rows.append({"nct_id": p.get("identificationModule", {}).get("nctId", ""),
                     "title": p.get("identificationModule", {}).get("briefTitle", "")[:160],
                     "phase": ",".join(p.get("designModule", {}).get("phases", []) or []),
                     "status": p.get("statusModule", {}).get("overallStatus", ""),
                     "conditions": ",".join(p.get("conditionsModule", {}).get("conditions", []) or [])[:120],
                     "url": "https://clinicaltrials.gov/study/" + p.get("identificationModule", {}).get("nctId", "")})
    return rows

def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]

def main():
    import urllib.parse as _  # ensure imported name exists at module level for helpers
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True, help="filtered variants TSV (needs gene, hgvsp, class)")
    ap.add_argument("--fusions", help="fusion summary TSV (gene1,gene2) to include in the therapy search")
    ap.add_argument("--cancer-type", default="")
    ap.add_argument("--classes", default="SOMATIC_LIKELY", help="comma-separated classes to annotate")
    ap.add_argument("--no-oncokb", action="store_true"); ap.add_argument("--no-trials", action="store_true")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    keep = set(a.classes.split(","))
    variants = [v for v in read_tsv(a.variants) if v.get("class", "SOMATIC_LIKELY") in keep]
    token = os.environ.get("ONCOKB_TOKEN", "")
    evidence, genes = [], set()
    for v in variants:
        gene = v.get("gene", "")
        pchange = short_protein_change(v.get("hgvsp", ""))
        if gene: genes.add(gene)
        for e in civic_for(gene, pchange):
            evidence.append({**e, "key": v["key"], "gene": gene, "alteration": pchange})
        if not a.no_oncokb:
            for e in oncokb_for(gene, pchange, a.cancer_type, token):
                evidence.append({**e, "key": v["key"], "gene": gene, "alteration": pchange})
        time.sleep(0.2)
    fusion_genes = []
    if a.fusions and Path(a.fusions).exists():
        for f in read_tsv(a.fusions):
            for g in (f.get("gene1", ""), f.get("gene2", "")):
                if g: genes.add(g); fusion_genes.append(g)
    drugs = dgidb_for(genes)
    trials = []
    if not a.no_trials:
        for g in sorted(genes)[:12]:
            for t in trials_for(g, a.cancer_type):
                trials.append({**t, "matched_gene": g})
            time.sleep(0.3)
    cols_e = ["key","gene","alteration","source","type","level","significance","direction","therapies","disease","url","description"]
    write_tsv(evidence, a.out + "_evidence.tsv", cols_e)
    write_tsv([{"gene": g, **d} for g, rows in drugs.items() for d in rows],
              a.out + "_druggability.tsv", ["gene","drug","score","types","sources"])
    write_tsv(trials, a.out + "_trials.tsv", ["matched_gene","nct_id","title","phase","status","conditions","url"])
    write_json({"n_variants_annotated": len(variants), "n_evidence_items": len(evidence),
                "genes_queried": sorted(genes), "oncokb_token_present": bool(token),
                "cancer_type": a.cancer_type,
                "note": "OncoKB annotations are licensed for academic use and must not be redistributed; "
                        "they are written under results/ which is git-ignored."},
               a.out + "_summary.json")
    print(f"evidence items: {len(evidence)}; genes: {len(genes)}; trials: {len(trials)}", file=sys.stderr)

if __name__ == "__main__":
    import urllib.parse
    main()
