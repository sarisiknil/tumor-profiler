#!/usr/bin/env python3
"""Upload the run's input files to Galaxy over the API, so the browser does not have to stay open.

Creates (or reuses) a history, uploads the FASTQ files and panel BEDs with the right datatype and genome build,
fetches the gnomAD germline resource server-side from its URL, and waits until everything is ready.

  export GALAXY_API_KEY=<from User -> Preferences -> Manage API Key>
  python3 scripts/galaxy_upload.py --check                     # verify the key and show the quota
  python3 scripts/galaxy_upload.py --history "TUMOR01 tumour profiling"

If anything goes wrong, the browser upload still works — see galaxy/UPLOAD_CHECKLIST.md.
"""
import argparse, json, os, sys, time
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("This script needs `requests`: conda activate tp-py, or pip install requests")

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SERVER = "https://usegalaxy.eu"
GNOMAD = "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"


def api(server, key, path, method="GET", **kw):
    url = f"{server}/api/{path}"
    kw.setdefault("params", {})["key"] = key
    r = requests.request(method, url, timeout=kw.pop("timeout", 300), **kw)
    if r.status_code >= 400:
        raise SystemExit(f"Galaxy API {method} {path} failed ({r.status_code}): {r.text[:400]}")
    return r.json() if r.content else {}


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def upload_file(server, key, history_id, path, ext, dbkey="hg38"):
    """Upload one local file with the fetch API (streams from disk; nothing is held in memory)."""
    path = Path(path)
    payload = {
        "history_id": history_id,
        "targets": json.dumps([{
            "destination": {"type": "hdas"},
            "elements": [{"src": "files", "name": path.name, "ext": ext, "dbkey": dbkey,
                          "auto_decompress": False}],
        }]),
        "auto_decompress": "false",
    }
    with open(path, "rb") as fh:
        files = {"files_0|file_data": (path.name, fh)}
        r = requests.post(f"{server}/api/tools/fetch", data=payload, files=files,
                          params={"key": key}, timeout=7200)
    if r.status_code >= 400:
        raise SystemExit(f"upload of {path.name} failed ({r.status_code}): {r.text[:400]}")
    out = r.json()
    return [d["id"] for d in out.get("outputs", [])]


def fetch_url(server, key, history_id, url, name, ext, dbkey="hg38"):
    """Ask Galaxy to download a URL itself — nothing passes through this machine."""
    payload = {"history_id": history_id,
               "targets": json.dumps([{"destination": {"type": "hdas"},
                                       "elements": [{"src": "url", "url": url, "name": name,
                                                     "ext": ext, "dbkey": dbkey}]}]),
               "auto_decompress": "false"}
    r = requests.post(f"{server}/api/tools/fetch", data=payload, params={"key": key}, timeout=600)
    if r.status_code >= 400:
        raise SystemExit(f"URL fetch failed ({r.status_code}): {r.text[:400]}")
    return [d["id"] for d in r.json().get("outputs", [])]


def wait(server, key, history_id, poll=20):
    print("\nwaiting for Galaxy to finish ingesting the datasets ...")
    while True:
        h = api(server, key, f"histories/{history_id}")
        counts = h.get("state_details", {}) or {}
        busy = counts.get("queued", 0) + counts.get("running", 0) + counts.get("new", 0) + counts.get("upload", 0)
        bad = counts.get("error", 0)
        print(f"  ok={counts.get('ok', 0)}  running/queued={busy}  error={bad}")
        if busy == 0:
            return bad == 0
        time.sleep(poll)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--server", default=os.environ.get("GALAXY_SERVER", DEFAULT_SERVER))
    ap.add_argument("--history", default="TUMOR01 tumour profiling")
    ap.add_argument("--config", default=str(REPO / "config" / "config.yaml"))
    ap.add_argument("--check", action="store_true", help="verify the API key and show the quota, then exit")
    ap.add_argument("--skip-gnomad", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="show what would be uploaded")
    a = ap.parse_args()
    key = os.environ.get("GALAXY_API_KEY", "")
    if not key:
        sys.exit("Set GALAXY_API_KEY first:  export GALAXY_API_KEY=<your key>\n"
                 "Get one at " + a.server + "/user/api_key")

    me = api(a.server, key, "users/current")
    print(f"signed in as {me.get('username')} ({me.get('email')}) on {a.server}")
    print(f"disk usage {me.get('nice_total_disk_usage', '?')} of quota {me.get('quota', '?')}")
    if a.check:
        return 0

    import yaml
    cfg = yaml.safe_load(open(a.config))
    d = cfg["data"]
    uploads = [(REPO / d["dna_r1"], "fastqsanger.gz"), (REPO / d["dna_r2"], "fastqsanger.gz"),
               (REPO / d["rna_r1"], "fastqsanger.gz"), (REPO / d["rna_r2"], "fastqsanger.gz"),
               (REPO / "resources" / "panel_dna.bed", "bed"),
               (REPO / "resources" / "panel_rna.bed", "bed")]
    missing = [str(p) for p, _ in uploads if not p.exists()]
    if missing:
        sys.exit("cannot find: " + ", ".join(missing))
    total = sum(p.stat().st_size for p, _ in uploads)
    print(f"\n{len(uploads)} files, {human(total)} total:")
    for p, ext in uploads:
        print(f"  {p.name:48s} {human(p.stat().st_size):>10s}  as {ext}")
    if not a.skip_gnomad:
        print(f"  + gnomAD germline resource, fetched by Galaxy itself from storage.googleapis.com")
    if a.dry_run:
        return 0

    hist = api(a.server, key, "histories", method="POST", json={"name": a.history})
    hid = hist["id"]
    print(f"\ncreated history {a.history!r} ({hid})")
    print(f"watch it at {a.server}/histories/view?id={hid}\n")

    for p, ext in uploads:
        t0 = time.time()
        print(f"uploading {p.name} ({human(p.stat().st_size)}) ...", flush=True)
        upload_file(a.server, key, hid, p, ext)
        print(f"  sent in {time.time()-t0:.0f}s")
    if not a.skip_gnomad:
        print("asking Galaxy to fetch the gnomAD germline resource ...")
        fetch_url(a.server, key, hid, GNOMAD, "af-only-gnomad.hg38.vcf.gz", "vcf_bgzip")

    ok = wait(a.server, key, hid)
    print("\nall datasets ready." if ok else "\nsome datasets failed — open the history and check them.")
    print(f"history id: {hid}\nnext: run the DNA and RNA workflows (galaxy/UPLOAD_CHECKLIST.md steps 4-5)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
