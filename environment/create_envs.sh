#!/usr/bin/env bash
# Creates the three conda environments (Bioconda / conda-forge). Idempotent: skips envs that exist.
set -u
cd "$(dirname "$0")/.."
for e in qc dna rna py; do
  name="tp-$e"
  if conda env list | awk '{print $1}' | grep -qx "$name"; then echo "[skip] $name exists"; continue; fi
  echo "[create] $name  ($(date))"
  conda env create -f environment/$e.yaml --yes > logs/env_$e.log 2>&1 && echo "[ok] $name" || { echo "[FAIL] $name -> see logs/env_$e.log"; tail -20 logs/env_$e.log; }
done
conda env list
