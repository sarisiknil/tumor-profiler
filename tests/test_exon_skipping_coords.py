"""Re-derive the splice-junction coordinates in scripts/exon_skipping.py from Ensembl and fail if they drift.

These constants are the kind that go wrong quietly: a wrong intron boundary makes the tool report
"no skipping detected" for a targetable event, which looks exactly like a true negative.
"""
import json
import sys
import urllib.request
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
UA = {"User-Agent": "tumor-profiler-tests/1.0", "Accept": "application/json"}


def _exons(symbol, transcript_id):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?expand=1;content-type=application/json"
    d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90))
    tx = next((t for t in d["Transcript"] if t["id"] == transcript_id), None)
    assert tx, f"{transcript_id} is no longer a transcript of {symbol}"
    return sorted(tx["Exon"], key=lambda e: e["start"])


@pytest.mark.network
def test_met_exon14_junctions():
    from exon_skipping import EVENTS
    ev = EVENTS["MET_exon14_skipping"]
    try:
        ex = _exons("MET", ev["transcript"])
    except Exception as e:
        pytest.skip(f"Ensembl unreachable: {e}")
    e13, e14, e15 = ex[12], ex[13], ex[14]
    assert e14["end"] - e14["start"] + 1 == 141, "MET exon 14 should be 141 bp"
    assert ev["skipping_junction"] == (e13["end"] + 1, e15["start"] - 1)
    assert (e13["end"] + 1, e14["start"] - 1) in ev["canonical_junctions"]
    assert (e14["end"] + 1, e15["start"] - 1) in ev["canonical_junctions"]


@pytest.mark.network
def test_egfrviii_junctions():
    from exon_skipping import EVENTS
    ev = EVENTS["EGFRvIII_exon2-7_deletion"]
    try:
        ex = _exons("EGFR", ev["transcript"])
    except Exception as e:
        pytest.skip(f"Ensembl unreachable: {e}")
    e1, e2, e8 = ex[0], ex[1], ex[7]
    assert ev["skipping_junction"] == (e1["end"] + 1, e8["start"] - 1)
    assert (e1["end"] + 1, e2["start"] - 1) in ev["canonical_junctions"]
