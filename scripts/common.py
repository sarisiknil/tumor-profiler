"""Shared helpers: config loading, VCF reading (no pysam needed), output paths, provenance stamping."""
from __future__ import annotations
import gzip, json, os, subprocess, sys, datetime, hashlib
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

def load_config(path: str | os.PathLike | None = None) -> dict:
    import yaml
    p = Path(path) if path else REPO / "config" / "config.yaml"
    with open(p) as fh:
        return yaml.safe_load(fh)

def open_maybe_gzip(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)

def read_vcf(path):
    """Yield dicts for each VCF record. Minimal, dependency-free; keeps INFO and the first sample's FORMAT."""
    samples = []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line[1:].split("\t")
                samples = cols[9:] if len(cols) > 9 else []
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            info = {}
            for kv in f[7].split(";"):
                if not kv:
                    continue
                k, _, v = kv.partition("=")
                info[k] = v if v else True
            rec = {"chrom": f[0], "pos": int(f[1]), "id": f[2], "ref": f[3], "alt": f[4],
                   "qual": f[5], "filter": f[6], "info": info, "fmt": {}}
            if len(f) > 9 and samples:
                keys = f[8].split(":")
                vals = f[9].split(":")
                rec["fmt"] = dict(zip(keys, vals))
            yield rec

def vaf_and_depth(rec):
    """Best-effort VAF/DP extraction across Mutect2 (AF/AD/DP), VarDict (AF/VD/DP), LoFreq (INFO AF/DP4/DP)."""
    fmt, info = rec["fmt"], rec["info"]
    dp = None; alt = None; vaf = None
    if "DP" in fmt:
        try: dp = int(fmt["DP"])
        except ValueError: pass
    if dp is None and "DP" in info:
        try: dp = int(info["DP"])
        except (ValueError, TypeError): pass
    if "AD" in fmt:
        parts = fmt["AD"].split(",")
        if len(parts) >= 2:
            try:
                ref_c, alt = int(parts[0]), int(parts[1])
                if dp is None: dp = ref_c + alt
            except ValueError: pass
    if alt is None and "VD" in fmt:
        try: alt = int(fmt["VD"])
        except ValueError: pass
    if "AF" in fmt:
        try: vaf = float(fmt["AF"].split(",")[0])
        except ValueError: pass
    if vaf is None and "AF" in info and info["AF"] is not True:
        try: vaf = float(str(info["AF"]).split(",")[0])
        except ValueError: pass
    if vaf is None and alt is not None and dp:
        vaf = alt / dp
    if alt is None and vaf is not None and dp:
        alt = round(vaf * dp)
    return vaf, dp, alt

def variant_key(rec) -> str:
    return f"{rec['chrom']}:{rec['pos']}:{rec['ref']}>{rec['alt']}"

def provenance(extra: dict | None = None) -> dict:
    def sh(cmd):
        try: return subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=REPO).stdout.strip()
        except Exception: return ""
    p = {"generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
         "git_commit": sh("git rev-parse --short HEAD"),
         "git_dirty": bool(sh("git status --porcelain")),
         "python": sys.version.split()[0],
         "argv": " ".join(sys.argv)}
    if extra: p.update(extra)
    return p

def write_json(obj, path, add_provenance=True):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    if add_provenance and isinstance(obj, dict):
        obj = {**obj, "_provenance": provenance()}
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
    return path

def write_tsv(rows, path, columns=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        Path(path).write_text("\t".join(columns or []) + "\n"); return path
    columns = columns or list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)).replace("\t", " ") for c in columns) + "\n")
    return path
