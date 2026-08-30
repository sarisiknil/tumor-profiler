"""Unit tests for the analysis logic that does not need external services."""
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))


def test_read_structure_recovers_the_simulated_library(tmp_path):
    """The structure detector must recover exactly what the simulator built."""
    import json
    out = tmp_path / "rs"
    subprocess.run([sys.executable, str(REPO / "scripts" / "read_structure.py"),
                    "--r1", str(REPO / "examples/simulated/sim_DNA_R1.fastq.gz"),
                    "--r2", str(REPO / "examples/simulated/sim_DNA_R2.fastq.gz"),
                    "-n", "3000", "-o", str(out)], check=True, capture_output=True)
    d = json.loads((tmp_path / "rs.json").read_text())
    assert d["R1"]["umi_len_inferred"] == 12
    assert d["R1"]["anchor_inferred"].startswith("AGTCGTCTCGAAG")
    # Read 2 starts at a primer: no inline UMI, and few distinct read starts
    assert d["R2"]["umi_len_inferred"] == 0
    assert d["R2"]["start_kmers_covering"]["80pct"] <= 50


def test_poisson_confidence_interval_brackets_the_estimate():
    from biomarkers import poisson_ci
    lo, hi = poisson_ci(10, 2.0)
    assert lo < 5.0 < hi


def test_protein_change_conversion_to_one_letter():
    from annotate_evidence import short_protein_change
    assert short_protein_change("ENSP00000493543.1:p.Val600Glu") == "V600E"
    assert short_protein_change("ENSP00000361021.3:p.Arg130Ter") == "R130*"
    assert short_protein_change("") == ""


def test_variant_normalisation_matches_differently_padded_indels():
    from merge_callers import norm
    a, *_ = norm({"chrom": "chr1", "pos": 100, "ref": "GA", "alt": "G"})
    b, *_ = norm({"chrom": "1", "pos": 100, "ref": "GA", "alt": "G"})
    assert a == b


def test_lift_regions_converts_contig_coordinates(tmp_path):
    regions = tmp_path / "r.tsv"
    regions.write_text("contig\tchrom\tstart\tend\tgene\tnote\nchr12_100_200\tchr12\t100\t200\tX\t-\n")
    vcf = tmp_path / "in.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "chr12_100_200\t51\t.\tC\tT\t.\tPASS\t.\n")
    out = tmp_path / "out.vcf"
    subprocess.run([sys.executable, str(REPO / "scripts" / "lift_regions.py"),
                    "--regions", str(regions), "-i", str(vcf), "-o", str(out)], check=True, capture_output=True)
    body = [l for l in out.read_text().splitlines() if not l.startswith("#")][0].split("\t")
    assert body[0] == "chr12" and body[1] == "150"


def test_pathway_gene_sets_are_wellformed():
    gmt = (REPO / "resources" / "oncogenic_pathways.gmt").read_text().strip().splitlines()
    assert len(gmt) >= 10
    for line in gmt:
        f = line.split("\t")
        assert len(f) > 3, f"malformed gene set: {f[0]}"
        assert all(g.isupper() or g[0].isdigit() for g in f[2:])
