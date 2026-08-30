"""Render every dashboard page headlessly and fail if any of them raises.

Streamlit's AppTest actually executes the page script, so this catches missing columns, bad JSON keys and
plotting errors that a plain HTTP check would miss.
"""
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PAGES = [REPO / "dashboard" / "app.py"] + sorted((REPO / "dashboard" / "pages").glob("*.py"))


@pytest.mark.parametrize("page", PAGES, ids=lambda p: p.name)
def test_page_renders(page, monkeypatch):
    from streamlit.testing.v1 import AppTest
    monkeypatch.chdir(REPO)
    at = AppTest.from_file(str(page), default_timeout=120)
    at.run()
    assert not at.exception, f"{page.name} raised: {at.exception}"
