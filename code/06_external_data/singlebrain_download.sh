#!/bin/bash
# Download SingleBrain single-nucleus cis-eQTL summary stats (7 major cell types).
#   Paper : Nat Genet 2026, doi:10.1038/s41588-026-02541-x (Raj lab)
#   Data  : https://zenodo.org/records/14908182 (CC-BY 4.0), 983 EUR donors, 4 cohorts
#   Build : hg38, GENCODE v38; variant_id is already an rsID.
# Subtypes (Ast1-4, MG1-4, ...) and MiGA3 are NOT downloaded -- add to CELL_TYPES to include.
# Run on the login node: bash singlebrain_download.sh
set -euo pipefail

RECORD=14908182
RAW_DIR="$HOME/db/shared_data/gwas/singlebrain/raw"
CELL_TYPES=(Ast End Ext IN MG OD OPC)

mkdir -p "$RAW_DIR"
cd "$RAW_DIR"

# Parallel, not a sequential loop: each full_assoc is 4-10 GB and one-at-a-time
# leaves the link idle waiting on a single stream. wget -c so a re-run resumes.
# setsid + </dev/null detaches each from the caller's process group, so they survive
# the parent shell going away.
PAR="${PAR:-4}"
for ct in "${CELL_TYPES[@]}"; do
  for assoc in full top; do
    f="${ct}_eqtl_${assoc}_assoc.tsv.gz"
    while [[ $(jobs -rp | wc -l) -ge "$PAR" ]]; do sleep 5; done
    echo "[singlebrain] start $f"
    setsid wget -q -c "https://zenodo.org/records/${RECORD}/files/${f}?download=1" -O "$f" </dev/null &
  done
done
wait

# Fail loudly: every expected file must exist and be non-empty.
missing=0
for ct in "${CELL_TYPES[@]}"; do
  for assoc in full top; do
    f="${RAW_DIR}/${ct}_eqtl_${assoc}_assoc.tsv.gz"
    [[ -s "$f" ]] || { echo "MISSING/EMPTY: $f" >&2; missing=1; }
  done
done
[[ "$missing" -eq 0 ]] || { echo "[singlebrain] download incomplete" >&2; exit 1; }

echo "[singlebrain] done: $(ls -1 "$RAW_DIR" | wc -l) files, $(du -sh "$RAW_DIR" | cut -f1)"
