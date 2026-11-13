#!/usr/bin/env bash
# Download and prepare 1KG Phase 3 hg38 EUR pfile.
# Source: https://www.cog-genomics.org/plink/2.0/resources
set -euo pipefail

OUT_DIR="${HOME}/db/shared_data/ref_geno/1kg_hg38"
mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

# --- Download ---
echo "==> Downloading 1KG hg38 files..."

curl -L -o all_hg38.pgen.zst \
  "https://www.dropbox.com/s/y6ytfoybz6ruoq0/all_hg38.pgen.zst?dl=1"

curl -L -o all_hg38_rs_noannot.pvar.zst \
  "https://www.dropbox.com/scl/fi/id642dpdd858uy41og8qi/all_hg38_rs_noannot.pvar.zst?rlkey=sskyiyam1bsqweujjmxqv1h55&dl=1"

curl -L -o hg38_corrected.psam \
  "https://www.dropbox.com/s/zj8d14vv0s6bi82/hg38_corrected.psam?dl=1"

# --- Decompress ---
echo "==> Decompressing..."
zstd -d all_hg38.pgen.zst -o all_hg38.pgen
zstd -d all_hg38_rs_noannot.pvar.zst -o all_hg38.pvar

# --- Extract EUR, autosomes only ---
echo "==> Extracting EUR samples (chr1-22)..."
awk 'NR>1 && $5=="EUR" {print $1}' hg38_corrected.psam > EUR_samples.txt

plink2 --pfile all_hg38 \
       --psam hg38_corrected.psam \
       --keep EUR_samples.txt \
       --chr chr1-chr22 \
       --allow-extra-chr \
       --make-pgen \
       --out EUR_hg38

echo "==> Generating EUR_hg38.bim (required for preimpute QC)..."
plink2 --pfile EUR_hg38 \
       --make-just-bim \
       --out EUR_hg38

echo "==> Cleaning up decompressed full-sample files..."
rm -f all_hg38.pgen all_hg38.pvar

echo ""
echo "=== Done ==="
echo "Output: ${OUT_DIR}/EUR_hg38.{pgen,pvar,psam,bim}"
echo "$(awk 'NR>1' EUR_hg38.psam | wc -l | tr -d ' ') EUR samples"
