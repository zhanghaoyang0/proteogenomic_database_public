#!/bin/bash
# Part 1 only: filter clean GWAS files to clump input (p < PTHR).
# IDX_START, IDX_END (1-based), COHORT, TYPE, PTHR passed via env.
# An input older than its clean file is redone; output goes via .tmp + mv, so a
# killed task never leaves a partial file that a re-run would keep.
set -euo pipefail

COHORT="${COHORT:-bf2}"
TYPE="${TYPE:-}"
PTHR="${PTHR:-5e-8}"
IDX_START="${IDX_START:-1}"
IDX_END="${IDX_END:-400}"

source ../lib/common.sh   # GWAS_BASE (run from clump/ by its sbatch)
GWAS_DIR="${GWAS_BASE}/${COHORT}"

# BF1/BF2: clean/<TYPE>/. FinnGen/other: flat clean/.
if [[ "$COHORT" == "finngen" || "$COHORT" == "other" ]]; then
  TYPE=""
  CLEAN_DIR="${GWAS_DIR}/clean"
  CLUMP_INPUT_DIR="${GWAS_DIR}/clump/input_${PTHR}Pthres"
else
  TYPE="${TYPE:-protein}"
  CLEAN_DIR="${GWAS_DIR}/clean/${TYPE}"
  CLUMP_INPUT_DIR="${GWAS_DIR}/clump/${TYPE}/input_${PTHR}Pthres"
fi

mkdir -p "$CLUMP_INPUT_DIR"

mapfile -t FILES < <(ls "$CLEAN_DIR"/*.tsv.gz 2>/dev/null)
TOTAL=${#FILES[@]}
IDX_END=$(( IDX_END < TOTAL ? IDX_END : TOTAL ))

echo "COHORT=$COHORT TYPE=$TYPE PTHR=$PTHR idx=${IDX_START}-${IDX_END} total=${TOTAL}"

for (( i=IDX_START-1; i<IDX_END; i++ )); do
  f="${FILES[$i]}"
  [[ -f "$f" ]] || continue
  protein=$(basename "$f" .tsv.gz)
  out="$CLUMP_INPUT_DIR/${protein}.tsv"
  if [[ -e "$out" && ! "$f" -nt "$out" ]]; then echo "  skip $protein"; continue; fi
  echo "  processing $protein"
  tmp="${out}.tmp"
  rm -f "$tmp"
  zcat "$f" | awk -F'\t' -v OFS="\t" -v pfile="$tmp" -v pthr="$PTHR" \
    'NR==1{for(i=1;i<=NF;i++){if($i=="P")pcol=i} if(!pcol){print "no P column" > "/dev/stderr"; exit 1} header=$0; next}
     # A non-numeric P (NA, empty) would become 0 under +0 and pass as significant;
     # a cleaned GWAS should have none, so stop instead of guessing.
     $pcol !~ /^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/ {print "non-numeric P \"" $pcol "\" at line " NR > "/dev/stderr"; exit 1}
     {p=$pcol+0; if(pthr+0>=1||p<=pthr+0) {if(!hprinted){print header>pfile; hprinted=1} print>pfile}}'
  # 0-byte file = no significant SNPs (05_clump.sh logs that, not a failure)
  [[ -e "$tmp" ]] || : > "$tmp"
  mv "$tmp" "$out"
done

echo "Done idx ${IDX_START}-${IDX_END}"
