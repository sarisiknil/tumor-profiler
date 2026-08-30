#!/usr/bin/env python3
"""Name the genes a panel targets, starting from the primers recovered out of Read 2.

Two routes, because the fast one needs an aligner:

  **--bed** (recommended)  The primer FASTA is aligned to GRCh38 - on Galaxy with BWA-MEM, or locally when a
  reference is at hand - and the resulting BED is annotated here through the Ensembl REST `/overlap/region`
  endpoint, which is fast, free and unrestricted. This is the default route because alignment is already part
  of the workflow.

  **--blast**  No aligner available: the primers are sent to NCBI BLAST's public URL API against human RefSeq
  RNA. Correct, but NCBI's public queue frequently quotes hours for a multi-hundred-sequence job, so treat it
  as a fallback and submit a small subset.

Outputs `<out>.tsv` (primer -> gene), `<out>_genes.tsv` (gene -> number of primers, i.e. how densely each gene
is covered) and `<out>_summary.json`.

Usage:
  python3 scripts/map_primers.py --bed results/primers/gsp2_RNA_mapped.bed -o results/primers/rna_panel
  python3 scripts/map_primers.py --fasta results/primers/gsp2_RNA.fa --blast --max 50 -o results/primers/rna_panel
"""
import argparse, json, re, sys, time, urllib.parse, urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv, write_json

BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "tumor-profiler/1.0 (educational bioinformatics pipeline)"}


def post(url, data):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode(), headers=UA)
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read().decode(errors="replace")


def get(url):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r:
        return r.read().decode(errors="replace")


def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def blast_submit(fasta_text, database="refseq_select_rna", program="blastn", word_size=16):
    body = post(BLAST, {"CMD": "Put", "PROGRAM": program, "MEGABLAST": "on", "DATABASE": database,
                        "QUERY": fasta_text, "ENTREZ_QUERY": "txid9606[ORGN]", "WORD_SIZE": str(word_size),
                        "HITLIST_SIZE": "3", "EXPECT": "1", "FILTER": "L",
                        "TOOL": "tumor-profiler", "FORMAT_TYPE": "Tabular"})
    m = re.search(r"RID = (\S+)", body)
    if not m:
        raise RuntimeError("NCBI BLAST did not return a request id; response began: " + body[:300])
    rtoe = re.search(r"RTOE = (\d+)", body)
    return m.group(1), int(rtoe.group(1)) if rtoe else 30


def blast_wait(rid, poll=20, max_wait=1800):
    waited = 0
    while waited < max_wait:
        body = get(f"{BLAST}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}")
        if "Status=READY" in body:
            return "ThereAreHits=yes" in body
        if "Status=FAILED" in body or "Status=UNKNOWN" in body:
            raise RuntimeError(f"BLAST job {rid} failed")
        time.sleep(poll); waited += poll
        print(f"  waiting for BLAST ({waited}s)", file=sys.stderr)
    raise TimeoutError(f"BLAST job {rid} did not finish within {max_wait}s")


def blast_results(rid):
    txt = get(f"{BLAST}?CMD=Get&RID={rid}&FORMAT_TYPE=Tabular&ALIGNMENTS=3&DESCRIPTIONS=3")
    best = {}
    for line in txt.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, s, pid, alen, bits = f[0], f[1], float(f[2]), int(f[3]), float(f[11])
        if q not in best or bits > best[q]["bits"]:
            best[q] = {"subject": s, "identity": pid, "aln_len": alen, "bits": bits}
    return best


def accession_titles(accessions, chunk=180):
    """Map RefSeq accessions to titles (which contain the gene symbol in parentheses)."""
    titles = {}
    accs = sorted({a.split("|")[-1] if "|" in a else a for a in accessions})
    for i in range(0, len(accs), chunk):
        ids = ",".join(accs[i:i + chunk])
        url = f"{EUTILS}/esummary.fcgi?db=nuccore&retmode=json&id={urllib.parse.quote(ids)}&tool=tumor-profiler"
        try:
            d = json.loads(get(url)).get("result", {})
        except Exception as e:
            print(f"  esummary failed: {e}", file=sys.stderr); continue
        for uid in d.get("uids", []):
            rec = d[uid]
            titles[rec.get("accessionversion", uid)] = rec.get("title", "")
            titles[rec.get("caption", uid)] = rec.get("title", "")
        time.sleep(0.4)
    return titles


def gene_from_title(title):
    m = re.search(r"\(([A-Z0-9orf\-]{2,12})\)[,\s]", title or "")
    return m.group(1) if m else ""


def read_bed(path):
    rows = []
    for line in open(path):
        if line.startswith(("#", "track", "browser")) or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        rows.append({"chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                     "primer": f[3] if len(f) > 3 else f"{f[0]}:{f[1]}",
                     "strand": f[5] if len(f) > 5 else "."})
    return rows


def genes_at(chrom, start, end, cache):
    """Gene symbols overlapping a region, via the Ensembl REST API (GRCh38)."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    key = (c, start // 10000)
    if key in cache:
        return cache[key]
    url = (f"https://rest.ensembl.org/overlap/region/human/{c}:{max(1, start)}-{end}"
           "?feature=gene;content-type=application/json")
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={**UA, "Accept": "application/json"}),
                                    timeout=60) as r:
            d = json.load(r)
    except Exception as e:
        print(f"    Ensembl overlap failed for {chrom}:{start}-{end}: {e}", file=sys.stderr)
        d = []
    gs = sorted({g.get("external_name") for g in d
                 if isinstance(g, dict) and g.get("external_name") and g.get("biotype") == "protein_coding"})
    cache[key] = gs
    return gs


def run_bed(bed_path, out):
    rows, cache = [], {}
    bed = read_bed(bed_path)
    print(f"annotating {len(bed)} mapped primers via Ensembl REST", file=sys.stderr)
    for i, b in enumerate(bed, start=1):
        gs = genes_at(b["chrom"], b["start"], b["end"], cache)
        rows.append({"primer": b["primer"], "gene": gs[0] if gs else "", "all_genes": ",".join(gs),
                     "chrom": b["chrom"], "start": b["start"], "end": b["end"], "strand": b["strand"]})
        if i % 100 == 0:
            print(f"  {i}/{len(bed)}", file=sys.stderr)
        time.sleep(0.08)
    write_tsv(rows, out + ".tsv", ["primer", "gene", "all_genes", "chrom", "start", "end", "strand"])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta")
    ap.add_argument("--bed", help="BED of primers already aligned to GRCh38 (e.g. Galaxy BWA-MEM output)")
    ap.add_argument("--blast", action="store_true", help="use the NCBI BLAST fallback instead of --bed")
    ap.add_argument("--max", type=int, default=600, help="map at most this many primers (ranked by read count)")
    ap.add_argument("--database", default="refseq_select_rna")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    if a.bed:
        rows = run_bed(a.bed, a.out)
        counts = Counter(r["gene"] for r in rows if r["gene"])
        write_tsv([{"gene": g, "n_primers": n} for g, n in counts.most_common()],
                  a.out + "_genes.tsv", ["gene", "n_primers"])
        write_json({"primers_annotated": len(rows), "primers_with_a_gene": sum(1 for r in rows if r["gene"]),
                    "n_genes": len(counts), "genes": sorted(counts),
                    "top_genes_by_primer_count": counts.most_common(60),
                    "method": "Primer coordinates from alignment; gene symbols from the Ensembl REST "
                              "/overlap/region endpoint (GRCh38, protein-coding genes).",
                    "caveats": ["Primer candidates are read-2 prefixes, so sequencing errors create "
                                "near-duplicates and the primer count per gene is an upper bound.",
                                "A primer in a paralogous region can overlap several genes; all are listed in "
                                "all_genes and the first is used for counting.",
                                "The recovered gene list approximates the panel design; it does not reproduce "
                                "the manufacturer's target file."]},
                   a.out + "_summary.json")
        print(f"{sum(1 for r in rows if r['gene'])}/{len(rows)} primers assigned to {len(counts)} genes",
              file=sys.stderr)
        print("top: " + ", ".join(f"{g}({n})" for g, n in counts.most_common(30)), file=sys.stderr)
        return
    if not a.fasta:
        ap.error("give --bed (recommended) or --fasta with --blast")
    primers = read_fasta(a.fasta)[: a.max]
    fasta_text = "".join(f">{n}\n{s}\n" for n, s in primers)
    print(f"submitting {len(primers)} primers to NCBI BLAST ({a.database})", file=sys.stderr)
    rid, rtoe = blast_submit(fasta_text, a.database)
    print(f"  request id {rid}, estimated {rtoe}s", file=sys.stderr)
    time.sleep(min(rtoe, 30))
    blast_wait(rid)
    best = blast_results(rid)
    titles = accession_titles({b["subject"] for b in best.values()})
    rows = []
    for name, seq in primers:
        b = best.get(name)
        if not b:
            rows.append({"primer": name, "sequence": seq, "gene": "", "transcript": "", "identity": "",
                         "aln_len": "", "title": ""})
            continue
        acc = b["subject"].split("|")[-1]
        title = titles.get(acc) or titles.get(acc.split(".")[0], "")
        rows.append({"primer": name, "sequence": seq, "gene": gene_from_title(title), "transcript": acc,
                     "identity": b["identity"], "aln_len": b["aln_len"], "title": title[:120]})
    write_tsv(rows, a.out + ".tsv", ["primer", "gene", "transcript", "identity", "aln_len", "title", "sequence"])
    counts = Counter(r["gene"] for r in rows if r["gene"])
    write_tsv([{"gene": g, "n_primers": n} for g, n in counts.most_common()],
              a.out + "_genes.tsv", ["gene", "n_primers"])
    mapped = sum(1 for r in rows if r["gene"])
    write_json({"primers_submitted": len(rows), "primers_with_a_gene": mapped, "n_genes": len(counts),
                "genes": sorted(counts), "top_genes_by_primer_count": counts.most_common(60),
                "blast_request_id": rid, "database": a.database,
                "method": "NCBI BLAST URL API (megablast vs human RefSeq RNA); the gene symbol is read from the "
                          "matching transcript's title.",
                "caveats": ["Primer candidates are read-2 prefixes: sequencing errors create near-duplicates, so "
                            "the primer count per gene is an upper bound.",
                            "A primer inside a paralogous or repetitive region can match several transcripts; "
                            "only the best hit is kept.",
                            "The recovered gene list approximates the panel design; it does not reproduce the "
                            "manufacturer's target file."]},
               a.out + "_summary.json")
    print(f"{mapped}/{len(rows)} primers assigned to {len(counts)} genes -> {a.out}_genes.tsv", file=sys.stderr)
    print("top: " + ", ".join(f"{g}({n})" for g, n in counts.most_common(30)), file=sys.stderr)


if __name__ == "__main__":
    main()
