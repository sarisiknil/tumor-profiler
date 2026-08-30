#!/usr/bin/env python3
"""Verify that every external tool the workflow needs is present AND executable on this machine.

Why this exists: Bioconda occasionally resolves a Linux build into an osx-arm64 environment, which then fails at
run time with a confusing "cannot execute binary file" (exit 126) in the middle of a workflow. Checking up front
turns that into a clear message. Run: python3 scripts/check_tools.py
"""
import shutil, subprocess, sys

TOOLS = {
    "bwa":          ["bwa"],                       # prints usage to stderr and exits 1 - handled below
    "samtools":     ["samtools", "--version"],
    "bcftools":     ["bcftools", "--version"],
    "umi_tools":    ["umi_tools", "--version"],
    "freebayes":    ["freebayes", "--version"],
    "lofreq":       ["lofreq", "version"],
    "vardict-java": ["vardict-java", "-h"],
    "mosdepth":     ["mosdepth", "--version"],
    "fastqc":       ["fastqc", "--version"],
    "fastp":        ["fastp", "--version"],
}
OPTIONAL = {"gatk": ["gatk", "--list"], "fgbio": ["fgbio", "--version"], "STAR": ["STAR", "--version"],
            "arriba": ["arriba", "-h"], "featureCounts": ["featureCounts", "-v"]}

def check(name, cmd):
    path = shutil.which(cmd[0])
    if not path:
        return False, "not found on PATH"
    try:
        r = subprocess.run(cmd, capture_output=True, timeout=120)
    except OSError as e:
        return False, f"{path}: {e}"
    except subprocess.TimeoutExpired:
        return False, f"{path}: timed out"
    if r.returncode == 126 or b"cannot execute binary file" in r.stderr:
        return False, f"{path}: wrong architecture for this machine (conda resolved a foreign build)"
    return True, path

def main():
    bad = []
    for name, cmd in TOOLS.items():
        ok, info = check(name, cmd)
        print(f"{'OK  ' if ok else 'FAIL'} {name:14s} {info}")
        if not ok:
            bad.append(name)
    print("--- optional ---")
    for name, cmd in OPTIONAL.items():
        ok, info = check(name, cmd)
        print(f"{'OK  ' if ok else '--  '} {name:14s} {info}")
    if bad:
        print(f"\nMissing or broken: {', '.join(bad)}. "
              f"Activate the right environment (environment/create_envs.sh) or reinstall those packages.",
              file=sys.stderr)
        return 1
    print("\nAll required tools are present and executable.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
