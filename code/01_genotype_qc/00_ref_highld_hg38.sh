#!/usr/bin/env bash
# Convert high-LD exclusion regions from hg19 to hg38 using liftOver.
# Output: ~/db/shared_data/ref_geno/1kg_v3/highld_hg38.range
set -euo pipefail

LIFTOVER="${HOME}/soft/liftOver"
CHAIN="${HOME}/db/shared_data/ref_geno/liftover/hg19ToHg38.over.chain.gz"
INPUT="${HOME}/db/shared_data/ref_geno/1kg_v3/highld_hg19.range"
OUT="${HOME}/db/shared_data/ref_geno/1kg_v3/highld_hg38.range"
TMP_BED="/tmp/highld_hg19.bed"
TMP_OUT="/tmp/highld_hg38.bed"
TMP_UNMAP="/tmp/highld_hg38_unmapped.bed"

# convert to 0-based BED with chr prefix
awk 'BEGIN{OFS="\t"} {print "chr"$1, $2-1, $3, $4}' "${INPUT}" > "${TMP_BED}"

# liftOver with relaxed minMatch to handle large regions spanning assembly gaps
"${LIFTOVER}" "${TMP_BED}" "${CHAIN}" "${TMP_OUT}" "${TMP_UNMAP}" -minMatch=0.1 || true

unmapped=$(grep -c '^[^#]' "${TMP_UNMAP}" || true)
unmapped="${unmapped:-0}"
echo "==> Unmapped regions: ${unmapped}"
if [ "${unmapped}" -gt 0 ] 2>/dev/null; then
  echo "==> Unmapped regions (review manually):"
  grep -v '^#' "${TMP_UNMAP}"
fi

# convert back to plink range format (1-based, no chr prefix)
awk 'BEGIN{OFS="\t"} {c=$1; sub(/^chr/,"",c); print c, $2+1, $3, $4}' \
  "${TMP_OUT}" > "${OUT}"

echo "==> Written: ${OUT}"
cat "${OUT}"
