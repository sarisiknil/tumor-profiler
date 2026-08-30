#!/usr/bin/env python3
"""Upload the run's input files to Galaxy, resumably, so the browser does not have to stay open.

Uses the TUS resumable-upload protocol through BioBlend (the official Galaxy Python client). A single-shot
multipart POST of a 600 MB FASTQ is fragile — a brief network hiccup aborts the whole transfer with an
HTTP 499 — whereas TUS sends the file in chunks and can pick up where it left off.

Files are uploaded under an alias (TUMOR01_DNA_R1.fastq.gz ...), so the laboratory's sample code never reaches
a public server: nothing downstream needs it, and the project's own rule is that identifiers stay local.

  export GALAXY_API_KEY=<from User -> Preferences -> Manage API Key>
  python3 scripts/galaxy_upload.py --check                        # verify the key and quota
  python3 scripts/galaxy_upload.py --history "TUMOR01 tumour profiling"
  python3 scripts/galaxy_upload.py --history "..." --only fastq   # resume just the FASTQ files

Interrupted transfers resume automatically: the TUS session URLs are kept in logs/galaxy_tus_state.json.
"""
import argparse, os, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

try:
    from bioblend.galaxy import GalaxyInstance
except ImportError:                     # the shell's `python3` may not be the environment's python
    import subprocess
    for cand in [os.environ.get("CONDA_PREFIX", "") + "/bin/python",
                 "/opt/miniconda3/envs/tp-py/bin/python",
                 str(Path.home() / "miniconda3/envs/tp-py/bin/python")]:
        if cand and Path(cand).exists() and cand != sys.executable:
            if subprocess.run([cand, "-c", "import bioblend"], capture_output=True).returncode == 0:
                print(f"[re-running with {cand}]", file=sys.stderr)
                os.execv(cand, [cand, __file__] + sys.argv[1:])
    sys.exit("This script needs BioBlend:\n"
             "  /opt/miniconda3/envs/tp-py/bin/python -m pip install bioblend")

DEFAULT_SERVER = "https://usegalaxy.eu"
GNOMAD = "https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--server", default=os.environ.get("GALAXY_SERVER", DEFAULT_SERVER))
    ap.add_argument("--history", default="TUMOR01 tumour profiling")
    ap.add_argument("--to-history", help="upload into this history id instead of creating/reusing by name")
    ap.add_argument("--config", default=str(REPO / "config" / "config.yaml"))
    ap.add_argument("--only", help="restrict to 'fastq' and/or 'beds' (comma-separated)")
    ap.add_argument("--check", action="store_true", help="verify the key and quota, then exit")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-gnomad", action="store_true")
    ap.add_argument("--chunk-mb", type=int, default=10,
                    help="usegalaxy.eu's proxy rejects chunks above ~10 MB with HTTP 413")
    ap.add_argument("--retries", type=int, default=4)
    a = ap.parse_args()

    key = os.environ.get("GALAXY_API_KEY", "")
    if not key:
        sys.exit(f"Set GALAXY_API_KEY first:  export GALAXY_API_KEY=<your key>\n"
                 f"Get one at {a.server}/user/api_key")
    gi = GalaxyInstance(url=a.server, key=key)
    me = gi.users.get_current_user()
    print(f"signed in as {me.get('username')} ({me.get('email')}) on {a.server}")
    print(f"disk usage {me.get('nice_total_disk_usage')} of quota {me.get('quota')}")
    if a.check:
        return 0

    import yaml
    cfg = yaml.safe_load(open(a.config))
    d = cfg["data"]
    alias = cfg.get("sample_id", "TUMOR01")
    uploads = [(REPO / d["dna_r1"], "fastqsanger.gz", f"{alias}_DNA_R1.fastq.gz"),
               (REPO / d["dna_r2"], "fastqsanger.gz", f"{alias}_DNA_R2.fastq.gz"),
               (REPO / d["rna_r1"], "fastqsanger.gz", f"{alias}_RNA_R1.fastq.gz"),
               (REPO / d["rna_r2"], "fastqsanger.gz", f"{alias}_RNA_R2.fastq.gz"),
               (REPO / "resources" / "panel_dna.bed", "bed", "panel_dna.bed"),
               (REPO / "resources" / "panel_rna.bed", "bed", "panel_rna.bed")]
    if a.only:
        keep = set(x.strip() for x in a.only.split(","))
        uploads = [u for u in uploads if ("beds" in keep and u[1] == "bed")
                   or ("fastq" in keep and u[1] != "bed")]
    missing = [str(p) for p, _, _ in uploads if not p.exists()]
    if missing:
        sys.exit("cannot find: " + ", ".join(missing))
    total = sum(p.stat().st_size for p, _, _ in uploads)
    print(f"\n{len(uploads)} files, {human(total)}:")
    for p, ext, nm in uploads:
        print(f"  {p.name[:40]:42s} -> {nm:26s} {human(p.stat().st_size):>10s}  as {ext}")
    if a.dry_run:
        return 0

    if a.to_history:
        hid = a.to_history
        print(f"\nusing history {hid}")
    else:
        existing = [h for h in gi.histories.get_histories(name=a.history)]
        hid = existing[0]["id"] if existing else gi.histories.create_history(name=a.history)["id"]
        print(f"\n{'reusing' if existing else 'created'} history {a.history!r} ({hid})")
    print(f"watch it at {a.server}/histories/view?id={hid}\n")

    # names already present and healthy in the history are not re-uploaded
    present = {}
    for ds in gi.histories.show_history(hid, contents=True, deleted=False):
        present.setdefault(ds.get("name"), []).append(ds.get("state"))
    # tuspy's storage is a FILE that records resume URLs, not a directory
    storage = REPO / "logs" / "galaxy_tus_state.json"
    storage.parent.mkdir(parents=True, exist_ok=True)

    for p, ext, nm in uploads:
        if any(s in ("ok", "running", "queued") for s in present.get(nm, [])):
            print(f"skip {nm} — already in the history")
            continue
        size = p.stat().st_size
        for attempt in range(1, a.retries + 1):
            t0 = time.time()
            print(f"uploading {nm} ({human(size)}) attempt {attempt}/{a.retries} ...", flush=True)
            try:
                gi.tools.upload_file(str(p), hid, file_name=nm, file_type=ext, dbkey="hg38",
                                     auto_decompress=False, to_posix_lines=(ext == "bed"),
                                     chunk_size=a.chunk_mb * 1000 * 1000,
                                     storage=str(storage))
                dt = time.time() - t0
                print(f"  done in {dt/60:.1f} min ({human(size/max(dt,1))}/s)", flush=True)
                break
            except Exception as e:
                print(f"  attempt {attempt} failed: {type(e).__name__}: {str(e)[:200]}", flush=True)
                if attempt == a.retries:
                    print(f"  giving up on {nm}; rerun with --only fastq to resume", flush=True)
                else:
                    time.sleep(15 * attempt)

    if not a.skip_gnomad and not any(s in ("ok", "running", "queued")
                                     for s in present.get("af-only-gnomad.hg38.vcf.gz", [])):
        print("asking Galaxy to fetch the gnomAD germline resource from Google ...", flush=True)
        try:
            gi.tools.put_url(GNOMAD, hid, file_name="af-only-gnomad.hg38.vcf.gz", file_type="vcf_bgzip")
        except Exception as e:
            print(f"  URL fetch failed ({e}); add it by hand with Paste/Fetch data:\n  {GNOMAD}")

    print("\nwaiting for Galaxy to finish ingesting ...", flush=True)
    while True:
        h = gi.histories.show_history(hid)
        c = h.get("state_details", {}) or {}
        busy = c.get("queued", 0) + c.get("running", 0) + c.get("new", 0) + c.get("upload", 0)
        print(f"  ok={c.get('ok', 0)}  busy={busy}  error={c.get('error', 0)}", flush=True)
        if busy == 0:
            break
        time.sleep(30)
    print(f"\nhistory id: {hid}")
    print("next: DNA and RNA workflows — galaxy/UPLOAD_CHECKLIST.md steps 3 and 4")
    return 0


if __name__ == "__main__":
    sys.exit(main())
