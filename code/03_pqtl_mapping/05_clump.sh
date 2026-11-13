#!/bin/bash
# Genome-wide PLINK clump, chunk mode.
# Required env: COHORT, TYPE, IDX_START, IDX_END
# A trait is done when its .clumped exists, or its .log says there was nothing
# to clump; both must be at least as new as the input. Anything else (a failed
# plink run included) is redone, and any failure makes the job exit 1.
#
# No high-LD region is excluded: highld_hg19.range is a PCA pruning list, and
# applying it here hid every signal in 136 Mb (402 BF2 assays could not get a cis
# lead). The MHC alone is special-cased: its leads collapse to the minimum-P one,
# so each trait has at most one MHC locus (rule: ../lib/mhc.R).
set -euo pipefail

COHORT="${COHORT:-bf2}"
TYPE="${TYPE:-}"
IDX_START="${IDX_START:-1}"
IDX_END="${IDX_END:-99999}"

source ../lib/common.sh   # GWAS_BASE, LD_REF (run from clump/ by its sbatch)
GWAS_DIR="${GWAS_BASE}/${COHORT}"
PTHR="${PTHR:-5e-8}"
CLUMP_KB="${CLUMP_KB:-1000}"
CLUMP_R2="${CLUMP_R2:-0.01}"

REF_GENO="$LD_REF"
MHC_RANGE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../lib" && pwd)/mhc_hg19.range"
[[ -s "$MHC_RANGE" ]] || { echo "ERROR: missing $MHC_RANGE" >&2; exit 1; }
read -r MHC_CHR MHC_START MHC_END _ < "$MHC_RANGE"

# Keep only the minimum-P lead inside the MHC (.clumped columns: CHR F SNP BP P ...).
collapse_mhc() {
  local f=$1
  awk -v c="$MHC_CHR" -v s="$MHC_START" -v e="$MHC_END" '
    function mhc() { return FNR > 1 && $1 == c && $4 >= s && $4 <= e }
    NR == FNR { if (mhc() && (best == "" || $5 + 0 < bestp)) { best = FNR; bestp = $5 + 0 }; next }
    !mhc() || FNR == best' "$f" "$f" > "$f.tmp"
  mv "$f.tmp" "$f"
}

# PLINK_BIN skips the module: aurora has none, but the module's static binary runs
# there from a copy in $HOME, e.g. PLINK_BIN=$HOME/bin/plink-1.9b_6.21 -p aurora.
if [[ -z "${PLINK_BIN:-}" ]]; then
  set +u   # Lmod reads unset vars
  module load PLINK/1.9b_6.21-x86_64
  set -u
  PLINK_BIN=$(command -v plink || true)
fi
[[ -x "$PLINK_BIN" ]] || { echo "ERROR: plink not found (module or PLINK_BIN)" >&2; exit 1; }

if [[ -n "$TYPE" ]]; then
  CLUMP_INPUT_DIR="${GWAS_DIR}/clump/${TYPE}/input_${PTHR}Pthres"
  CLUMP_OUTPUT_DIR="${GWAS_DIR}/clump/${TYPE}/output_${PTHR}Pthres_${CLUMP_KB}kb_${CLUMP_R2}r2"
else
  CLUMP_INPUT_DIR="${GWAS_DIR}/clump/input_${PTHR}Pthres"
  CLUMP_OUTPUT_DIR="${GWAS_DIR}/clump/output_${PTHR}Pthres_${CLUMP_KB}kb_${CLUMP_R2}r2"
fi
mkdir -p "$CLUMP_OUTPUT_DIR"

mapfile -t all_files < <(ls -r "$CLUMP_INPUT_DIR"/*.tsv 2>/dev/null)
n_total=${#all_files[@]}
idx_end=$(( IDX_END > n_total ? n_total : IDX_END ))
start=$(( IDX_START - 1 ))
count=$(( idx_end - IDX_START + 1 ))
files=("${all_files[@]:$start:$count}")

echo "Processing traits ${IDX_START}-${idx_end} (n=${#files[@]}) of $n_total"

if awk "BEGIN{exit !($PTHR+0>=1)}"; then CLUMP_P1=1; else CLUMP_P1=$PTHR; fi

processed=0; failed=0; skipped=0
for f in "${files[@]}"; do
  [[ -f "$f" ]] || continue
  protein=$(basename "$f" .tsv)
  outpref="$CLUMP_OUTPUT_DIR/$protein"
  if [[ -s "${outpref}.clumped" && ! "$f" -nt "${outpref}.clumped" ]]; then skipped=$((skipped+1)); continue; fi
  if [[ -s "${outpref}.log" && ! "$f" -nt "${outpref}.log" ]] &&
     grep -qE "^No significant SNPs in input|No significant --clump results" "${outpref}.log"; then
    skipped=$((skipped+1)); continue
  fi
  rm -f "${outpref}.clumped"
  if [[ ! -s "$f" ]]; then
    echo "No significant SNPs in input." > "${outpref}.log"
    skipped=$((skipped+1)); continue
  fi
  echo "Clump $protein..."
  if "$PLINK_BIN" --bfile "$REF_GENO" \
       --clump "$f" --clump-snp-field SNP --clump-field P \
       --clump-p1 "$CLUMP_P1" --clump-p2 1 --clump-r2 "$CLUMP_R2" --clump-kb "$CLUMP_KB" \
       --out "$outpref" > "${outpref}.log" 2>&1; then
    [[ -s "${outpref}.clumped" ]] && collapse_mhc "${outpref}.clumped"
    processed=$((processed+1)); echo "  [OK] $protein"
  else
    failed=$((failed+1)); echo "  [FAIL] $protein (see ${outpref}.log)"
  fi
done
rm -f "$CLUMP_OUTPUT_DIR"/*.nosex
echo "Done: processed=$processed skipped=$skipped failed=$failed"
[[ "$failed" -eq 0 ]] || { echo "ERROR: $failed trait(s) failed to clump" >&2; exit 1; }
