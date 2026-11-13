#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DATASET="bf1_hg19"
# Raw 2018 delivery, 1680 samples: the lab's QC'd raw/GRch37 pfile (1546) plus
# the 134 it marked QC FAIL. Our QC decides who goes; until 2026-11-02 the input
# was that pfile (run deleted 2026-11-03; its key is a subset of the current one).
INPUT_TYPE="bfile"
INPUT_PREFIX="${HOME}/db/secure_data/bf/array/bf1/raw/BioFINDER1RAWDATA/GSP_2018-054-biofinder"
OUT_DIR="${HOME}/db/secure_data/bf/array/bf1/for_impute"

source "${SCRIPT_DIR}/lib/qc_hg19.sh"
run_preimpute_qc
