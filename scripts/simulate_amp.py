#!/usr/bin/env python3
"""Generate a small anchored-multiplex-PCR (AMP) dataset on a REAL (mini) GRCh38 reference, so the whole
pipeline — including live annotation against VEP/ClinVar/gnomAD/OncoKB/CIViC — can be run, taught and
CI-tested WITHOUT any patient data.

Reproduces the read architecture of an Archer-style library:
    Read 1: [12-nt random UMI][fixed common region][T][insert ...]
    Read 2: [gene-specific primer][insert ...]
Spikes real driver hotspots at real GRCh38 coordinates and writes a truth set so tests can assert that the
pipeline recovers them. PCR duplicates are simulated (several reads per UMI family) so that UMI deduplication
has something to do, and a configurable sequencing error rate makes consensus calling meaningful.

Prerequisite: examples/reference/mini_GRCh38.fa (scripts/fetch_example_reference.py).
Usage: simulate_amp.py -o examples/simulated [--depth 400] [--seed 7]
"""
import argparse, gzip, json, random
from pathlib import Path

DNA_COMMON = "AGTCGTCTCGAAG"
RNA_COMMON = "CTGGATAGTACGCT"

# gene, genomic chrom, 1-based position, ref, alt, target VAF, protein change (for the truth set)
HOTSPOTS = [
    ("KRAS",   "chr12", 25245351, "C", "A", 0.32, "p.Gly12Cys"),  # c.34G>T on the minus strand
    ("BRAF",   "chr7",  140753336, "A", "T", 0.18, "p.Val600Glu"),
    ("TP53",   "chr17", 7675088,  "C", "T", 0.55, "p.Arg175His"),
    ("EGFR",   "chr7",  55191822, "T", "G", 0.07, "p.Leu858Arg"),   # subclonal, tests low-VAF sensitivity
    ("PIK3CA", "chr3",  179234297, "A", "G", 0.24, "p.His1047Arg"),
    ("PTEN",   "chr10", 87933147, "C", "T", 0.40, "p.Arg130Ter"),
]

def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]

def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

def read_regions(path):
    out = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                out.append({"contig": f[0], "chrom": f[1], "start": int(f[2]), "end": int(f[3]), "gene": f[4]})
    return out

def sequencing_errors(rng, seq, rate):
    if rate <= 0:
        return seq
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--outdir", default="examples/simulated")
    ap.add_argument("--reference", default="examples/reference/mini_GRCh38.fa")
    ap.add_argument("--regions", default="examples/reference/regions.tsv")
    ap.add_argument("--depth", type=int, default=400, help="unique molecules per amplicon")
    ap.add_argument("--dup-rate", type=float, default=2.5, help="mean reads per UMI family (PCR duplicates)")
    ap.add_argument("--error-rate", type=float, default=0.002)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    ref = read_fasta(a.reference)
    regions = read_regions(a.regions)
    by_gene = {}
    for r in regions:
        by_gene.setdefault(r["gene"], []).append(r)

    r1, r2, primer_recs, truth, bed = [], [], [], [], []
    rid = 0
    for gene, regs in by_gene.items():
        reg = regs[0]
        contig, seq = reg["contig"], ref[reg["contig"]]
        spikes = [h for h in HOTSPOTS if h[0] == gene]
        # one amplicon centred on each hotspot, plus one background amplicon
        centres = [(h[2] - reg["start"]) for h in spikes] or [len(seq) // 2]
        centres.append(min(len(seq) - 400, max(400, centres[0] + 3000)))
        for ai, centre in enumerate(centres):
            p_off = max(0, centre - 60)                        # primer sits ~60 bp 5' of the hotspot
            primer = seq[p_off:p_off + 25]
            if "N" in primer or len(primer) < 25:
                continue
            primer_recs.append((f"{gene}_amp{ai+1}", primer, reg["chrom"], reg["start"] + p_off))
            bed.append((contig, p_off, min(len(seq), p_off + a.read_len), f"{gene}_amp{ai+1}"))
            spike = next((h for h in spikes if p_off <= (h[2] - reg["start"]) < p_off + a.read_len), None)
            for mol in range(a.depth):
                template = seq[p_off:p_off + a.read_len + 40]
                carries = False
                if spike:
                    off = (spike[2] - reg["start"]) - p_off
                    if rng.random() < spike[5]:
                        template = template[:off] + spike[4] + template[off + len(spike[3]):]
                        carries = True
                umi = rand_seq(rng, 12)
                n_dups = max(1, int(rng.expovariate(1 / a.dup_rate)))
                for d in range(n_dups):
                    rid += 1
                    ins = sequencing_errors(rng, template, a.error_rate)
                    read2 = ins[:a.read_len]
                    read1 = (umi + DNA_COMMON + "T" + rc(ins)[:a.read_len - 26])[:a.read_len]
                    nm = (f"SIM:1:EXAMPLEFC:1:1101:{rid}:{mol} "
                          "{}:N:0:ACGTACGT+ACGTACGTAC")
                    r1.append((nm.format(1), read1, "F" * len(read1)))
                    r2.append((nm.format(2), read2, "F" * len(read2)))
                if spike and carries:
                    pass
        for h in spikes:
            truth.append({"gene": h[0], "chrom": h[1], "pos": h[2], "ref": h[3], "alt": h[4],
                          "expected_vaf": h[5], "protein_change": h[6],
                          "contig": reg["contig"], "contig_pos": h[2] - reg["start"] + 1})

    with gzip.open(out / "sim_DNA_R1.fastq.gz", "wt") as f1, gzip.open(out / "sim_DNA_R2.fastq.gz", "wt") as f2:
        for (n, s, q), (n2, s2, q2) in zip(r1, r2):
            f1.write(f"@{n}\n{s}\n+\n{q}\n"); f2.write(f"@{n2}\n{s2}\n+\n{q2}\n")

    # ---- RNA: an EML4-ALK-like fusion plus a MET-exon-14-skipping-like event on synthetic transcripts -----
    rr1, rr2 = [], []
    def rna_reads(tx, n, tag):
        for i in range(n):
            pos = rng.randint(0, max(0, len(tx) - a.read_len))
            ins = sequencing_errors(rng, tx[pos:pos + a.read_len], a.error_rate)
            umi = rand_seq(rng, 12)
            nm = f"SIM:2:EXAMPLEFC:1:1101:{rng.randint(1,10**6)}:{i}_{tag}"
            read1 = (umi + RNA_COMMON + rc(ins)[:a.read_len - 26])[:a.read_len]
            rr1.append((nm + " 1:N:0:ACGTACGT+ACGTACGTAC", read1, "F" * len(read1)))
            rr2.append((nm + " 2:N:0:ACGTACGT+ACGTACGTAC", ins, "F" * len(ins)))
    donor, acceptor = rand_seq(rng, 300), rand_seq(rng, 300)
    rna_reads(donor + acceptor, 150, "FUSION")
    rna_reads(donor + rand_seq(rng, 300), 400, "WT5p")
    ex13, ex14, ex15 = rand_seq(rng, 200), rand_seq(rng, 141), rand_seq(rng, 200)
    rna_reads(ex13 + ex15, 120, "SKIP")
    rna_reads(ex13 + ex14 + ex15, 130, "CANON")
    with gzip.open(out / "sim_RNA_R1.fastq.gz", "wt") as f1, gzip.open(out / "sim_RNA_R2.fastq.gz", "wt") as f2:
        for (n, s, q), (n2, s2, q2) in zip(rr1, rr2):
            f1.write(f"@{n}\n{s}\n+\n{q}\n"); f2.write(f"@{n2}\n{s2}\n+\n{q2}\n")

    with open(out / "sim_primers.fa", "w") as fh:
        for nm, s, c, p in primer_recs:
            fh.write(f">{nm} {c}:{p}\n{s}\n")
    with open(out / "sim_targets.bed", "w") as fh:
        for c, s, e, n in sorted(bed):
            fh.write(f"{c}\t{s}\t{e}\t{n}\n")
    (out / "truth.json").write_text(json.dumps(
        {"description": "Synthetic anchored-multiplex-PCR dataset on a real mini-GRCh38 reference. "
                        "Contains no patient data.",
         "reference": str(a.reference), "regions": str(a.regions),
         "read_structure": {"read1": f"12-nt UMI + {DNA_COMMON} + T + insert",
                            "read2": "gene-specific primer + insert",
                            "rna_read1": f"12-nt UMI + {RNA_COMMON} + insert"},
         "dna_read_pairs": len(r1), "rna_read_pairs": len(rr1),
         "parameters": {"molecules_per_amplicon": a.depth, "mean_reads_per_umi": a.dup_rate,
                        "sequencing_error_rate": a.error_rate, "seed": a.seed},
         "expected_variants": truth,
         "expected_fusion": {"supporting_reads": 150, "note": "synthetic transcripts; tests the fusion module "
                                                              "plumbing, not real ALK biology"},
         "expected_exon_skipping": {"skipping_reads": 120, "canonical_reads": 130,
                                    "approx_ratio": round(120 / 250, 2)}}, indent=1))
    print(f"DNA: {len(r1):,} read pairs over {len(primer_recs)} amplicons; RNA: {len(rr1):,} read pairs")
    print(f"truth set: {len(truth)} spiked hotspots -> {out}/truth.json")

if __name__ == "__main__":
    main()
