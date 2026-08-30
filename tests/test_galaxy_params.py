"""Check that every Galaxy parameter path in scripts/galaxy_run.py still exists on the server.

Galaxy tools get updated and parameters get renamed. This test turns "the run silently used a default" into a
red build. It needs network access, so it is skipped when Galaxy is unreachable.
"""
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SERVER = "https://usegalaxy.eu"
UA = {"User-Agent": "tumor-profiler-tests/1.0", "Accept": "application/json"}


def _get(url):
    return json.load(urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=90))


def _flatten(inputs, prefix=""):
    paths = set()
    for i in inputs or []:
        full = f"{prefix}{i.get('name', '')}"
        t = i.get("type")
        if t == "conditional":
            paths.add(full)
            tp = (i.get("test_param") or {}).get("name")
            if tp:
                paths.add(f"{full}|{tp}")
            for c in i.get("cases", []) or []:
                paths |= _flatten(c.get("inputs"), f"{full}|")
        elif t in ("section", "repeat"):
            paths |= _flatten(i.get("inputs"), f"{full}|")
        else:
            paths.add(full)
    return paths


def _blocks():
    src = (REPO / "scripts" / "galaxy_run.py").read_text()
    tools = dict(re.findall(r'"(\w+)":\s*"(toolshed[^"]+)"', src))
    out = []
    for key, label, body in re.findall(r'r\.run\("(\w+)",\s*"([^"]+)",\s*\{(.*?)\n    \},', src, re.S):
        out.append((key, tools[key], label, re.findall(r'^\s*"([^"]+)":', body, re.M)))
    return out


@pytest.mark.network
@pytest.mark.parametrize("key,tool_id,label,params", _blocks(), ids=lambda x: x if isinstance(x, str) else "")
def test_parameters_exist(key, tool_id, label, params):
    try:
        tool = _get(f"{SERVER}/api/tools/{urllib.parse.quote(tool_id, safe='')}?io_details=true")
    except Exception as e:                                  # offline or Galaxy down
        pytest.skip(f"cannot reach {SERVER}: {e}")
    legal = _flatten(tool.get("inputs"))
    bad = [p for p in params if p not in legal]
    assert not bad, f"{label}: Galaxy no longer has {bad} (tool {tool_id})"
