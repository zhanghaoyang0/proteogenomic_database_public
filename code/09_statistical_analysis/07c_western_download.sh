#!/usr/bin/env bash
#
# Downloads the Western et al. 2026 CSF Olink summary statistics that 07c_western_replication.R
# listed (tmp/proteogenomic_database/07c_western_download.tsv: accession, url) and keeps only the rows whose rs_id
# is one of our CSF credible-set variants (tmp/proteogenomic_database/07c_western_rsids.txt). Each file (~190 MB)
# is streamed and filtered, never stored whole; the extract is written to a temporary name
# and renamed when complete, so a rerun skips finished studies and redoes broken ones.
#
# Usage (repo root, after the first pass of 07c_western_replication.R):
#   sbatch code/09_statistical_analysis/07c_western_download.sh
#   PILOT=2 bash code/09_statistical_analysis/07c_western_download.sh   # first 2 only
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

TMP_DIR=tmp/proteogenomic_database
LIST=$TMP_DIR/07c_western_download.tsv
RSIDS=$TMP_DIR/07c_western_rsids.txt
LOOKUP_DIR=$HOME/db/shared_data/pqtl_ref/western_csf_olink/lookup
JOBS=${SLURM_CPUS_PER_TASK:-4}
for f in "$LIST" "$RSIDS"; do [[ -s $f ]] || { echo "missing $f: run 07c_western_replication.R first" >&2; exit 1; }; done
mkdir -p "$LOOKUP_DIR"

fetch_one() {
  set -euo pipefail
  local acc=$1 url=$2 out=$LOOKUP_DIR/$1.tsv.gz
  [[ -s $out ]] && return 0
  # column 10 must be rs_id (GWAS-SSF layout of these files); stop otherwise
  curl -sS --fail --retry 3 --retry-delay 20 "$url" | zcat |
    awk -F'\t' -v acc="$acc" 'NR == FNR { keep[$1]; next }
      FNR == 1 { if ($10 != "rs_id") { print acc ": column 10 is " $10 > "/dev/stderr"; exit 2 }; print; next }
      $10 in keep' "$RSIDS" - | gzip > "$out.part" || return 1
  # set -e above is off here: fetch_one runs left of || in the xargs call below
  mv "$out.part" "$out"
}
export -f fetch_one
export LOOKUP_DIR RSIDS

n=$(wc -l < "$LIST")
[[ -n ${PILOT:-} ]] && n=$PILOT
echo "$(date) fetching $n studies with $JOBS parallel streams -> $LOOKUP_DIR"
head -n "$n" "$LIST" | xargs -P "$JOBS" -n 2 bash -c 'fetch_one "$1" "$2" || { echo "FAILED $1" >&2; exit 255; }' _
done_n=$(head -n "$n" "$LIST" | cut -f1 | while read -r a; do [[ -s $LOOKUP_DIR/$a.tsv.gz ]] && echo; done | wc -l)
echo "$(date) $done_n of $n extracts present"
[[ $done_n -eq $n ]] || { echo "some extracts are missing" >&2; exit 1; }
echo "07c_western_download done"
