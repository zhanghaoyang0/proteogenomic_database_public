#!/bin/bash
set -euo pipefail

# Download GTEx V8 cis-eQTL full summary statistics (all 49 tissues) from the
# eQTL Catalogue (EBI reprocessing of GTEx V8 — free, anonymous, rsIDs included).
#   Catalogue: https://www.ebi.ac.uk/eqtl/
#   FTP:       https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/{QTS}/{QTD}/{QTD}.all.tsv.gz
#
# Datasets to fetch are listed in info/gtex/datasets.csv (dataset_id=QTD, study_id=QTS,
# sample_group, tissue_label, sample_size), generated from the eQTL Catalogue API.
# We pull gene-expression (ge) cis-eQTL "all" files (full cis sumstats, hg38) + .tbi.
#
# Why not GTEx's own bucket: the official GTEx_Analysis_v8_eQTL_all_associations now
# lives in a requester-pays GCP bucket (gs://gtex-resources). eQTL Catalogue is the
# free equivalent. Run on the login node (NOT via sbatch). wget -c resumes partials.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${SCRIPT_DIR}/../../info/gtex/datasets.csv" # from https://www.ebi.ac.uk/eqtl/api/v2/datasets
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"
RAW_DIR="${HOME}/db/shared_data/gwas/gtex/raw"
mkdir -p "${RAW_DIR}"
cd "${RAW_DIR}"

[ -f "${MANIFEST}" ] || { echo "Missing manifest: ${MANIFEST}" >&2; exit 1; }

# skip header; iterate datasets
tail -n +2 "${MANIFEST}" | while IFS=, read -r qtd qts sample_group tissue_label sample_size; do
  url="${BASE_URL}/${qts}/${qtd}/${qtd}.all.tsv.gz"
  echo "==> ${sample_group} (${qtd})"
  wget -c "${url}"
  wget -c "${url}.tbi"
done

echo "Done. Per-dataset {QTD}.all.tsv.gz under ${RAW_DIR}"
