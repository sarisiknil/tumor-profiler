#!/usr/bin/env python3
"""Characterize the read structure of an Illumina library directly from FASTQ (no alignment, no dependencies).

Answers, from the first N reads:
  * read-length distribution, N content
  * per-position base composition -> detects a random UMI prefix (uniform bases) followed by a fixed
    "common region"/anchor (one base >90 %), as in Archer AMP MBC adapters
  * adapter read-through rate and position (proxy for short inserts / FFPE fragmentation)
  * how concentrated the read starts are (distinct k-mers needed to cover 25/50/80 % of reads):
    a handful of k-mers covering most reads = reads start at gene-specific primers (targeted panel)
  * index sequences from the header (i7+i5 lengths are a kit fingerprint, e.g. 8+10 = IDT Liquid P5 MBC + P7)

Usage: read_structure.py --r1 R1.fastq.gz [--r2 R2.fastq.gz] [-n 200000] -o out_prefix
Writes <out_prefix>.json and <out_prefix>.md
"""
import argparse, gzip, json, sys, collections, statistics

ADAPTER = "AGATCGGAAGAGC"   # Illumina TruSeq/Nextera universal adapter stem (both reads)

def open_fq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def iter_reads(path, n):
    with open_fq(path) as fh:
        i = 0
        while i < n:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip("\n"); fh.readline(); q = fh.readline().rstrip("\n")
            yield h.rstrip("\n"), s, q
            i += 1

def analyze(path, n, k=20, npos=40):
    comp = [collections.Counter() for _ in range(npos)]
    lengths = collections.Counter()
    adapter_pos = collections.Counter()
    start_kmers = collections.Counter()
    index_seqs = collections.Counter()
    n_reads = 0; n_with_N = 0; n_adapter = 0
    instrument = None
    for h, s, q in iter_reads(path, n):
        n_reads += 1
        lengths[len(s)] += 1
        if "N" in s: n_with_N += 1
        for i, b in enumerate(s[:npos]):
            comp[i][b] += 1
        p = s.find(ADAPTER)
        if p >= 0:
            n_adapter += 1; adapter_pos[p] += 1
        start_kmers[s[:k]] += 1
        # header: @INSTR:RUN:FLOWCELL:LANE:TILE:X:Y READ:FILTER:0:INDEX(+INDEX2)
        if instrument is None:
            instrument = h[1:].split(":")[0]
        parts = h.split(" ")
        if len(parts) > 1 and parts[1].count(":") >= 3:
            index_seqs[parts[1].split(":")[3]] += 1
    # positional consensus: fraction of the most common base at each position
    consensus = []
    for i, c in enumerate(comp):
        tot = sum(c.values()) or 1
        b, cnt = c.most_common(1)[0]
        consensus.append({"pos": i + 1, "base": b, "frac": round(cnt / tot, 4),
                          "A": round(c["A"]/tot, 3), "C": round(c["C"]/tot, 3), "G": round(c["G"]/tot, 3), "T": round(c["T"]/tot, 3)})
    # Structure inference: a random UMI is a leading run of positions where no base dominates (<50 %),
    # immediately followed by a fixed "common region" (>=90 % one base). If no such fixed region follows,
    # there is no inline UMI (e.g. Read 2 of a one-sided panel, which starts at a gene-specific primer):
    # in that case report umi_len = 0 and no anchor rather than the length of the scanned window.
    lead = 0
    while lead < npos and consensus[lead]["frac"] < 0.5:
        lead += 1
    anchor = ""
    j = lead
    while j < npos and consensus[j]["frac"] >= 0.9:
        anchor += consensus[j]["base"]; j += 1
    anchor_tail = ""   # positions with 0.75-0.9 consensus right after the anchor (e.g. ligation T-overhang)
    while j < npos and consensus[j]["frac"] >= 0.75:
        anchor_tail += consensus[j]["base"]; j += 1
    if len(anchor) < 6 or lead >= npos:      # no credible fixed common region -> no inline UMI
        umi_len, anchor, anchor_tail = 0, "", ""
    else:
        umi_len = lead
    # k-mer concentration at read start
    tot = sum(start_kmers.values()) or 1
    cum = 0; q25 = q50 = q80 = None
    for rank, (_, c) in enumerate(start_kmers.most_common(), start=1):
        cum += c
        if q25 is None and cum >= 0.25 * tot: q25 = rank
        if q50 is None and cum >= 0.50 * tot: q50 = rank
        if q80 is None and cum >= 0.80 * tot: q80 = rank; break
    top_kmers = [{"kmer": km, "count": c, "frac": round(c / tot, 4)} for km, c in start_kmers.most_common(10)]
    idx = index_seqs.most_common(1)[0][0] if index_seqs else ""
    idx_lens = [len(x) for x in idx.split("+")] if idx else []
    return {
        "file": path, "reads_analyzed": n_reads, "instrument": instrument,
        "index_sequence": idx, "index_lengths": idx_lens,
        "read_length": {"min": min(lengths), "max": max(lengths),
                         "mode": lengths.most_common(1)[0][0],
                         "frac_full_length": round(lengths.most_common(1)[0][1] / n_reads, 4)},
        "frac_reads_with_N": round(n_with_N / n_reads, 4),
        "adapter": {"stem": ADAPTER, "frac_reads_with_adapter": round(n_adapter / n_reads, 4),
                    "median_adapter_position": (statistics.median(adapter_pos.elements()) if n_adapter else None)},
        "umi_len_inferred": umi_len, "anchor_inferred": anchor, "anchor_tail_inferred": anchor_tail,
        "insert_starts_at_base": umi_len + len(anchor) + len(anchor_tail) + 1,
        "start_kmer_k": k, "distinct_start_kmers": len(start_kmers),
        "start_kmers_covering": {"25pct": q25, "50pct": q50, "80pct": q80},
        "top_start_kmers": top_kmers, "positional_consensus": consensus,
    }

def md_summary(res, label):
    r = res
    lines = [f"### {label}: `{r['file'].split('/')[-1]}`", "",
             f"- Instrument `{r['instrument']}`; index read `{r['index_sequence']}` (lengths {r['index_lengths']})",
             f"- Reads analyzed: {r['reads_analyzed']:,}; read length mode {r['read_length']['mode']} "
             f"({r['read_length']['frac_full_length']*100:.1f} % full length; min {r['read_length']['min']})",
             f"- Adapter stem found in {r['adapter']['frac_reads_with_adapter']*100:.1f} % of reads "
             f"(median position {r['adapter']['median_adapter_position']})",
             (f"- Inferred UMI length **{r['umi_len_inferred']}**, fixed common region **`{r['anchor_inferred']}`**"
              + (f" + `{r['anchor_tail_inferred']}`" if r['anchor_tail_inferred'] else "")
              + f" → insert starts at base {r['insert_starts_at_base']}")
             if r['umi_len_inferred'] else
             "- No inline UMI / fixed common region detected (read starts directly in target sequence, "
             "as expected for Read 2 of a one-sided gene-specific-primer panel)",
             f"- Read-start {r['start_kmer_k']}-mers: {r['distinct_start_kmers']:,} distinct; "
             f"{r['start_kmers_covering']['25pct']} / {r['start_kmers_covering']['50pct']} / {r['start_kmers_covering']['80pct']} "
             f"k-mers cover 25 / 50 / 80 % of reads", ""]
    return "\n".join(lines)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True); ap.add_argument("--r2")
    ap.add_argument("-n", type=int, default=200000); ap.add_argument("-k", type=int, default=20)
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    a = ap.parse_args()
    out = {"R1": analyze(a.r1, a.n, a.k)}
    if a.r2: out["R2"] = analyze(a.r2, a.n, a.k)
    json.dump(out, open(a.out + ".json", "w"), indent=1)
    with open(a.out + ".md", "w") as fh:
        for lab, res in out.items(): fh.write(md_summary(res, lab) + "\n")
    print(open(a.out + ".md").read())

if __name__ == "__main__":
    main()
