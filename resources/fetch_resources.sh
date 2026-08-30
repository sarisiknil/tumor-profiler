#!/usr/bin/env bash
# Downloads reference resources that are free to obtain but that we do NOT redistribute in this repository
# (see docs/licensing.md). Run once after cloning:  bash resources/fetch_resources.sh
set -euo pipefail
cd "$(dirname "$0")"
UA="tumor-profiler/1.0 (educational pipeline)"

# 1) OncoKB Cancer Gene List — downloadable without a licence; not redistributed here (git-ignored).
if [ ! -s cancer_genes.tsv ]; then
  echo "[fetch] OncoKB cancer gene list"
  curl -fsSL -A "$UA" "https://www.oncokb.org/api/v1/utils/cancerGeneList.txt" -o cancer_genes.tsv
fi
awk -F'\t' 'NR>1 && $1!="" {print $1"\t"$7}' cancer_genes.tsv > cancer_genes.txt
echo "[ok] $(wc -l < cancer_genes.txt) cancer genes -> resources/cancer_genes.txt (symbol, gene type)"

# 2) ClinVar VCF for GRCh38 (public domain, ~90 MB) — only needed for offline annotation; the default
#    pipeline reads ClinVar significance through the Ensembl VEP REST API instead.
if [ "${WITH_CLINVAR:-0}" = "1" ] && [ ! -s clinvar_GRCh38.vcf.gz ]; then
  echo "[fetch] ClinVar GRCh38 VCF"
  curl -fsSL -A "$UA" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz" -o clinvar_GRCh38.vcf.gz
fi

# 3) Cancer Hotspots (Chang et al. 2016/2018), free for research use.
if [ "${WITH_HOTSPOTS:-0}" = "1" ] && [ ! -s hotspots.tsv ]; then
  echo "[fetch] cancerhotspots.org v2 single-residue hotspots"
  curl -fsSL -A "$UA" "https://www.cancerhotspots.org/files/hotspots_v2.xls" -o hotspots_v2.xls || true
fi
echo "[done]"
