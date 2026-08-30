"""End-to-end assertions on the example run: the pipeline must recover what was spiked into the data.

These tests are the reason the synthetic dataset exists. They are the difference between "the workflow ran"
and "the workflow is right".
"""
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RES = REPO / "results_example"
TRUTH = json.loads((REPO / "examples" / "simulated" / "truth.json").read_text())

pytestmark = pytest.mark.skipif(not (RES / "summary.json").exists(),
                                reason="run `snakemake -c4` first to produce results_example/")


def read_tsv(p):
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def test_all_spiked_variants_are_called():
    called = {r["key"] for r in read_tsv(RES / "dna" / "variants_filtered.tsv")}
    missing = [f"{v['gene']} {v['chrom']}:{v['pos']}{v['ref']}>{v['alt']}"
               for v in TRUTH["expected_variants"]
               if f"{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}" not in called]
    assert not missing, f"spiked variants not recovered: {missing}"


def test_spiked_variants_are_classified_somatic():
    rows = {r["key"]: r for r in read_tsv(RES / "dna" / "variants_filtered.tsv")}
    for v in TRUTH["expected_variants"]:
        r = rows[f"{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}"]
        assert r["class"] == "SOMATIC_LIKELY", f"{v['gene']} classified {r['class']}: {r['reasons']}"


def test_variant_allele_fractions_are_close_to_the_truth():
    rows = {r["key"]: r for r in read_tsv(RES / "dna" / "variants_filtered.tsv")}
    for v in TRUTH["expected_variants"]:
        r = rows[f"{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}"]
        obs, exp = float(r["vaf"]), v["expected_vaf"]
        assert abs(obs - exp) < 0.15, f"{v['gene']} VAF {obs:.3f} vs expected {exp}"


def test_annotation_assigns_the_expected_gene_and_protein_change():
    rows = {r["key"]: r for r in read_tsv(RES / "dna" / "vep.tsv")}
    for v in TRUTH["expected_variants"]:
        r = rows[f"{v['chrom']}:{v['pos']}:{v['ref']}>{v['alt']}"]
        assert r["gene"] == v["gene"], f"expected {v['gene']}, VEP said {r['gene']}"
        expected_p = v["protein_change"].replace("p.", "")
        assert expected_p in r["hgvsp"], f"{v['gene']}: expected {expected_p}, got {r['hgvsp']}"


def test_known_actionable_variants_reach_tier_I_or_II():
    """KRAS G12C, BRAF V600E, EGFR L858R and PIK3CA H1047R all have level-A CIViC evidence."""
    rows = {r["gene"]: r for r in read_tsv(RES / "dna" / "variants_tiered.tsv")}
    for gene in ("KRAS", "BRAF", "EGFR", "PIK3CA"):
        assert gene in rows, f"{gene} missing from the tiered table"
        assert rows[gene]["amp_tier"] in ("I", "II"), \
            f"{gene} got tier {rows[gene]['amp_tier']} ({rows[gene]['tier_rationale']})"


def test_pathway_mapping_finds_the_rtk_ras_pathway():
    rows = {r["pathway"]: r for r in read_tsv(RES / "pathways" / "pathway_hits.tsv")}
    assert int(rows["RTK_RAS"]["n_altered"]) >= 3
    assert int(rows["PI3K"]["n_altered"]) >= 1


def test_summary_and_report_exist_and_carry_the_disclaimer():
    s = json.loads((RES / "summary.json").read_text())
    assert "disclaimer" in s and "not a clinical report" in s["disclaimer"].lower()
    assert (RES / "report.md").exists()


def test_no_patient_identifiers_in_shareable_outputs():
    """Nothing in the example results may contain the lab's sample identifier."""
    for p in RES.rglob("*"):
        if p.is_file() and p.suffix in (".json", ".tsv", ".md"):
            assert "REDACTED" not in p.read_text(errors="ignore"), f"identifier leaked into {p}"
