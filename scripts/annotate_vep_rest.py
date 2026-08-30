#!/usr/bin/env python3
"""Annotate variants with Ensembl VEP through its REST API (no 20 GB offline cache needed).

Why REST: a targeted panel yields hundreds, not millions, of variants. rest.ensembl.org allows ~55,000
requests/hour and its POST endpoint takes 200 variants at a time, so a whole panel annotates in seconds and
the pipeline stays runnable on a laptop. For WES/WGS scale you would install the offline cache instead.

Input : merged variant TSV from merge_callers.py (columns chrom,pos,ref,alt) or a VCF.
Output: <out>.tsv (one row per variant, most severe consequence in the canonical/MANE transcript) and <out>.json
        (full VEP payload, kept for reproducibility/debugging).

Fields added: gene symbol, HGVSc/HGVSp, consequence, impact, canonical/MANE transcript, exon, protein position,
              gnomAD popmax AF (germline filtering), ClinVar significance, SIFT/PolyPhen, variant class.
"""
import argparse, json, sys, time, urllib.request, urllib.error
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_vcf, write_tsv, write_json

SERVER = "https://rest.ensembl.org"
CHUNK = 180

def post(server, endpoint, payload, tries=4):
    data = json.dumps(payload).encode()
    req = urllib.request.Request(server + endpoint, data=data,
                                 headers={"Content-Type": "application/json", "Accept": "application/json",
                                          "User-Agent": "tumor-profiler/1.0"})
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 429:                      # rate limited: obey Retry-After
                wait = int(e.headers.get("Retry-After", "5")) + 1
                print(f"  rate-limited, sleeping {wait}s", file=sys.stderr); time.sleep(wait); continue
            if e.code >= 500 and attempt < tries - 1:
                time.sleep(3 * (attempt + 1)); continue
            raise
        except Exception:
            if attempt < tries - 1:
                time.sleep(3 * (attempt + 1)); continue
            raise
    raise RuntimeError("VEP REST request failed after retries")

def to_vep_notation(chrom, pos, ref, alt):
    """Ensembl 'region' input format: 'chr pos . ref alt . . .' (VCF-like, 1-based)."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return f"{c} {pos} . {ref} {alt} . . ."

def read_input(path):
    rows = []
    if str(path).endswith((".vcf", ".vcf.gz")):
        for rec in read_vcf(path):
            rows.append({"chrom": rec["chrom"], "pos": rec["pos"], "ref": rec["ref"],
                         "alt": rec["alt"].split(",")[0], "key": f"{rec['chrom']}:{rec['pos']}:{rec['ref']}>{rec['alt'].split(',')[0]}"})
    else:
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if not f.get("chrom"):
                    continue
                f["pos"] = int(f["pos"])
                rows.append(f)
    return rows

def pick_transcript(tcs):
    """Prefer MANE Select, then canonical, then the most severe consequence, then the first."""
    if not tcs:
        return {}
    for key in ("mane_select", "mane"):
        m = [t for t in tcs if t.get(key)]
        if m:
            return m[0]
    can = [t for t in tcs if t.get("canonical")]
    return (can or tcs)[0]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="merged variants TSV or VCF")
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument("--server", default=SERVER)
    ap.add_argument("--assembly", default="GRCh38", choices=["GRCh38", "GRCh37"])
    a = ap.parse_args()
    server = a.server if a.assembly == "GRCh38" else "https://grch37.rest.ensembl.org"
    variants = read_input(a.input)
    if not variants:
        write_tsv([], a.out + ".tsv"); write_json({"n": 0}, a.out + ".json")
        print("no variants to annotate", file=sys.stderr); return
    print(f"annotating {len(variants)} variants via {server}", file=sys.stderr)
    raw = []
    for i in range(0, len(variants), CHUNK):
        batch = variants[i:i + CHUNK]
        payload = {"variants": [to_vep_notation(v["chrom"], v["pos"], v["ref"], v["alt"]) for v in batch],
                   "canonical": 1, "hgvs": 1, "mane": 1, "numbers": 1, "domains": 1,
                   "sift": "b", "polyphen": "b", "af_gnomadg": 1, "af_gnomade": 1, "af": 1,
                   "variant_class": 1, "check_existing": 1, "clinvar": 1, "pick_order": "mane_select,canonical,rank"}
        raw.extend(post(server, "/vep/homo_sapiens/region", payload))
        print(f"  {min(i+CHUNK, len(variants))}/{len(variants)}", file=sys.stderr)
        time.sleep(0.4)
    by_input = {r.get("input", ""): r for r in raw}
    rows = []
    for v in variants:
        r = by_input.get(to_vep_notation(v["chrom"], v["pos"], v["ref"], v["alt"]), {})
        t = pick_transcript(r.get("transcript_consequences", []))
        cons = t.get("consequence_terms") or r.get("most_severe_consequence", "")
        col = r.get("colocated_variants", []) or []
        # population frequency: max over gnomAD genome/exome population maxima reported by VEP
        af_fields, clinvar, dbsnp = [], set(), ""
        for c in col:
            for k, val in (c.get("frequencies", {}) or {}).items():
                for pop, af in (val or {}).items():
                    if pop.startswith("gnomad") and isinstance(af, (int, float)):
                        af_fields.append(af)
            for k in ("gnomadg", "gnomade", "gnomadg_af", "af"):
                if isinstance(c.get(k), (int, float)):
                    af_fields.append(c[k])
            for cs in (c.get("clin_sig") or []):
                clinvar.add(cs)
            if str(c.get("id", "")).startswith("rs") and not dbsnp:
                dbsnp = c["id"]
        rows.append({**{k: v[k] for k in ("key", "chrom", "pos", "ref", "alt") if k in v},
                     "gene": t.get("gene_symbol", ""), "transcript": t.get("transcript_id", ""),
                     "mane": bool(t.get("mane_select")), "canonical": bool(t.get("canonical")),
                     "consequence": ",".join(cons) if isinstance(cons, list) else cons,
                     "impact": t.get("impact", ""), "exon": t.get("exon", ""), "intron": t.get("intron", ""),
                     "hgvsc": t.get("hgvsc", ""), "hgvsp": t.get("hgvsp", ""),
                     "protein_position": t.get("protein_start", ""), "amino_acids": t.get("amino_acids", ""),
                     "variant_class": r.get("variant_class", ""),
                     "gnomad_af_max": round(max(af_fields), 6) if af_fields else 0.0,
                     "clinvar_sig": ",".join(sorted(clinvar)), "dbsnp": dbsnp,
                     "sift": t.get("sift_prediction", ""), "polyphen": t.get("polyphen_prediction", ""),
                     "most_severe_consequence": r.get("most_severe_consequence", "")})
    cols = ["key","chrom","pos","ref","alt","gene","consequence","impact","hgvsc","hgvsp","protein_position",
            "amino_acids","exon","intron","transcript","mane","canonical","variant_class","gnomad_af_max",
            "clinvar_sig","dbsnp","sift","polyphen","most_severe_consequence"]
    write_tsv(rows, a.out + ".tsv", cols)
    write_json({"n_variants": len(rows), "server": server, "assembly": a.assembly, "vep_raw": raw}, a.out + ".json")
    print(f"wrote {a.out}.tsv ({len(rows)} rows)", file=sys.stderr)

if __name__ == "__main__":
    main()
