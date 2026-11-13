#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DATASET="bf2_hg38"
INPUT_TYPE="bfile"
INPUT_PREFIX="${HOME}/db/secure_data/bf/array/bf2/raw/UB-2844_240229_PLINK_PCF_PLUS/1_BF-2_Genetic_Data"
OUT_DIR="${HOME}/db/secure_data/bf/array/bf2/for_impute"

source "${SCRIPT_DIR}/lib/qc_hg38.sh"
run_preimpute_qc_hg38
