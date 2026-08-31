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
    """Observed variant allele fraction, depth and alternate-read count, across caller conventions.

    Order matters. Allele *counts* come first because they are unambiguous, and only then the AF fields:
    LoFreq's INFO/AF is the observed frequency, but FreeBayes' INFO/AF is a genotype-model estimate that is
    0.5 for every heterozygous call. Reading that as a VAF silently replaces every FreeBayes measurement with
    0.5 — which is exactly the kind of wrong number that looks plausible in a table.
    """
    fmt, info = rec["fmt"], rec["info"]
    dp = alt = vaf = None

    def as_int(x):
        try:
            return int(x)
        except (TypeError, ValueError):
            return None

    # 1. allele depths: AD=ref,alt (GATK, VarDict, FreeBayes) -- the most direct measurement
    if "AD" in fmt:
        parts = fmt["AD"].split(",")
        if len(parts) >= 2:
            ref_c, alt_c = as_int(parts[0]), as_int(parts[1])
            if ref_c is not None and alt_c is not None and (ref_c + alt_c) > 0:
                alt, dp, vaf = alt_c, ref_c + alt_c, alt_c / (ref_c + alt_c)
    # 2. FreeBayes also gives reference/alternate observation counts explicitly
    if vaf is None and "RO" in fmt and "AO" in fmt:
        ref_c, alt_c = as_int(fmt["RO"]), as_int(str(fmt["AO"]).split(",")[0])
        if ref_c is not None and alt_c is not None and (ref_c + alt_c) > 0:
            alt, dp, vaf = alt_c, ref_c + alt_c, alt_c / (ref_c + alt_c)
    # 3. VarDict reports the alternate count as VD alongside DP
    if vaf is None and "VD" in fmt and "DP" in fmt:
        alt_c, d = as_int(fmt["VD"]), as_int(fmt["DP"])
        if alt_c is not None and d:
            alt, dp, vaf = alt_c, d, alt_c / d
    # 4. an explicit per-sample AF (Mutect2) is a real measurement
    if vaf is None and "AF" in fmt:
        try:
            vaf = float(fmt["AF"].split(",")[0])
        except ValueError:
            pass
    # 5. INFO/DP4 = ref-fwd, ref-rev, alt-fwd, alt-rev (LoFreq, bcftools)
    if vaf is None and isinstance(info.get("DP4"), str):
        parts = [as_int(x) for x in info["DP4"].split(",")]
        if len(parts) == 4 and all(p is not None for p in parts):
            ref_c, alt_c = parts[0] + parts[1], parts[2] + parts[3]
            if ref_c + alt_c > 0:
                alt, dp, vaf = alt_c, ref_c + alt_c, alt_c / (ref_c + alt_c)
    # 6. last resort: INFO/AF. Correct for LoFreq, a genotype estimate for FreeBayes, so it is only used
    #    when nothing above produced a value and never for a record that carried allele counts.
    if vaf is None and info.get("AF") not in (None, True):
        try:
            vaf = float(str(info["AF"]).split(",")[0])
        except ValueError:
            pass
    if dp is None:
        dp = as_int(fmt.get("DP")) or as_int(info.get("DP"))
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
