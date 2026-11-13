#!/bin/bash
# Baseline protein GWAS: one plink2 --glm pass per panel (no age grouping).
# Skips any protein whose .glm.linear.zst already exists in the output dir,
# so re-runs are idempotent and previously-finished panels add zero work.
#
# All panels are inverse-normal-transformed (--pheno-quantile-normalize) by
# plink2 before --glm, so every platform shares a common phenotype scale.
#
# Panel list = baseline (Visit==0) panels with gwas sample > MIN_GWAS_SAMPLES.
# Env overrides: ONLY_COHORT, ONLY_PANEL, CHUNK_SIZE, THREADS, MEM_MB,
#                MIN_GWAS_SAMPLES, QC_SUMMARY, PANEL_MANIFEST,
#                NO_MEAN_PROTEIN (set to 1 to drop mean_protein_{panel} covar).

set -uo pipefail

BF_ROOT="$HOME/db/secure_data/bf"
PROTEIN_DIR="$BF_ROOT/array/gwas_input/protein"
GWAS_DIR="$HOME/db/shared_data/gwas"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$SCRIPT_DIR/log"

PFILE_BF1="${PFILE_BF1:-$BF_ROOT/array/bf1/imputed_qc/bf1_hg19}"
PFILE_BF2="${PFILE_BF2:-$BF_ROOT/array/bf2/imputed_qc/bf2_hg19}"
CHUNK_SIZE="${CHUNK_SIZE:-200}"
THREADS="${THREADS:-6}"
MEM_MB="${MEM_MB:-10000}"
ONLY_COHORT="${ONLY_COHORT:-}"
ONLY_PANEL="${ONLY_PANEL:-}"
NO_MEAN_PROTEIN="${NO_MEAN_PROTEIN:-0}"
MIN_GWAS_SAMPLES="${MIN_GWAS_SAMPLES:-1000}"
QC_SUMMARY="${QC_SUMMARY:-$SCRIPT_DIR/../../info/bf_protein/qc_panel_summary.csv}"
PANEL_MANIFEST="${PANEL_MANIFEST:-}"
FAILED=0

mkdir -p "$LOG_DIR"

join_by_comma() {
    local IFS=,
    echo "$*"
}

has_glm_zst() {
    # Treat protein as already-run if either the zst output OR a sibling .done
    # sentinel exists. The sentinel lets us delete the (large) zst after it
    # has been archived elsewhere while preserving the skip logic.
    local out_dir=$1 out_tag=$2 protein_col=$3
    compgen -G "$out_dir/${out_tag}.${protein_col}.glm.linear.zst" >/dev/null \
        || compgen -G "$out_dir/${out_tag}.${protein_col}.glm.linear.zst.done" >/dev/null
}

run_panel() {
    local cohort=$1 tag=$2 pfile=$3
    if [ -n "$ONLY_COHORT" ] && [ "$cohort" != "$ONLY_COHORT" ]; then return; fi
    if [ -n "$ONLY_PANEL" ] && [ "$tag" != "$ONLY_PANEL" ]; then return; fi

    local pheno="$PROTEIN_DIR/${cohort}_${tag}.phe"
    local covar="$PROTEIN_DIR/${cohort}_covar.phe"
    local out_dir="$GWAS_DIR/$cohort/raw/protein"
    mkdir -p "$out_dir"

    # Every listed panel needs its inputs (01 writes a .phe for each one it lists).
    if [ ! -f "$pheno" ]; then echo "ERROR: $cohort $tag: missing $pheno"; FAILED=1; return; fi
    if [ ! -f "$covar" ]; then echo "ERROR: $cohort $tag: missing $covar"; FAILED=1; return; fi

    # Protein columns = .phe header minus FID, IID.
    local proteins=()
    while IFS= read -r p; do proteins+=("$p"); done < <(
        head -1 "$pheno" | tr '\t' '\n' | tail -n +3
    )
    local n=${#proteins[@]}
    [ "$n" -eq 0 ] && { echo "ERROR: $cohort $tag: no protein columns in $pheno"; FAILED=1; return; }

    # Covariate names: age, sex, pc1..pc10, plus mean_protein_{tag} when present.
    local covar_names="age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10"
    if [ "$NO_MEAN_PROTEIN" != "1" ]; then
        if head -1 "$covar" | tr '\t' '\n' | grep -qx "mean_protein_${tag}"; then
            # A present-but-empty covariate column passes every existence check and
            # then silently drops every sample at --glm time. Count real values.
            local mp_n
            mp_n=$(awk -F'\t' -v col="mean_protein_${tag}" '
                NR==1 { for (i = 1; i <= NF; i++) if ($i == col) k = i; next }
                k && $k != "NA" && $k != "" && $k != "nan" { n++ }
                END { print n + 0 }
            ' "$covar")
            if [ "$mp_n" -eq 0 ]; then
                echo "ERROR: $cohort $tag: mean_protein_${tag} is present in $covar but has no values" >&2
                FAILED=1
                return 1
            fi
            local n_pheno_samples
            n_pheno_samples=$(($(wc -l < "$pheno") - 1))
            if [ "$mp_n" -lt $((n_pheno_samples / 2)) ]; then
                echo "WARNING: $cohort $tag: mean_protein_${tag} covers only $mp_n of $n_pheno_samples phenotyped samples" >&2
            fi
            covar_names="${covar_names},mean_protein_${tag}"
        fi
    fi

    echo "=============================="
    echo "$cohort $tag: $n proteins | chunk=$CHUNK_SIZE threads=$THREADS mem=$MEM_MB"
    echo "covar_names=$covar_names"
    echo "=============================="

    local chunk_n=0 start=0
    while [ "$start" -lt "$n" ]; do
        local chunk=()
        local j=0
        while [ "$j" -lt "$CHUNK_SIZE" ] && [ $((start + j)) -lt "$n" ]; do
            local p="${proteins[$((start + j))]}"
            if ! has_glm_zst "$out_dir" "$tag" "$p"; then
                chunk+=("$p")
            fi
            j=$((j + 1))
        done
        start=$((start + j))
        [ "${#chunk[@]}" -eq 0 ] && continue
        chunk_n=$((chunk_n + 1))

        local pheno_names
        pheno_names=$(join_by_comma "${chunk[@]}")
        local out_prefix="$out_dir/${tag}"
        local last_idx=$((${#chunk[@]} - 1))
        echo "  Chunk $chunk_n: ${#chunk[@]} proteins (${chunk[0]} ... ${chunk[$last_idx]})"

        # INT (rank-based inverse-normal) all panel phenotypes via plink2, so
        # every platform (Olink, NULISA, SomaLogic) is on a common scale.
        # --no-input-missing-phenotype: plink2 otherwise reads -9 as missing,
        # and log-scale protein values can be exactly -9. Missing is written as
        # NA by 01/02 (na_rep='NA'), so nothing relies on the -9 code.
        plink2 --pfile "$pfile" \
               --no-input-missing-phenotype \
               --pheno "$pheno" \
               --pheno-name "$pheno_names" \
               --pheno-quantile-normalize \
               --covar "$covar" \
               --covar-name "$covar_names" \
               --covar-variance-standardize \
               --glm zs hide-covar cols=chrom,pos,omitted,a1freq,nobs,beta,se,p omit-ref \
               --snps-only just-acgt \
               --maf 0.01 --geno 0.05 dosage --max-alleles 2 \
               --threads "$THREADS" \
               --memory "$MEM_MB" \
               --out "$out_prefix"
        local status=$?
        local raw_log="$out_prefix.log"
        local code_log="$LOG_DIR/${cohort}_${tag}_chunk_${chunk_n}.plink2.log"
        [ -f "$raw_log" ] && mv "$raw_log" "$code_log"
        if [ "$status" -ne 0 ]; then
            echo "  plink2 failed (status=$status): $code_log"
            FAILED=1
        fi
    done
}

run_panels_from_summary() {
    if [ ! -f "$QC_SUMMARY" ]; then
        echo "missing QC summary: $QC_SUMMARY"
        FAILED=1
        return
    fi
    if ! command -v python3 >/dev/null 2>&1; then
        echo "python3 not found; cannot parse $QC_SUMMARY"
        FAILED=1
        return
    fi

    echo "Panel source: $QC_SUMMARY"
    echo "Panel filter: gwas sample > $MIN_GWAS_SAMPLES and gwas protein > 0"

    local n_listed=0
    while IFS=$'\t' read -r cohort tag n_samples n_proteins; do
        n_listed=$((n_listed + 1))
        case "$cohort" in
            bf1) run_panel "$cohort" "$tag" "$PFILE_BF1" ;;
            bf2) run_panel "$cohort" "$tag" "$PFILE_BF2" ;;
            *) echo "skip unknown cohort in summary: $cohort $tag" ;;
        esac
    done < <(python3 - "$QC_SUMMARY" "$MIN_GWAS_SAMPLES" <<'PY'
import csv
import re
import sys

summary_path = sys.argv[1]
min_samples = int(sys.argv[2])
gwas_col = "gwas (sample x protein)"

with open(summary_path, newline="") as fh:
    for row in csv.DictReader(fh):
        match = re.match(r"\s*(\d+)\s*x\s*(\d+)\s*$", row.get(gwas_col, "") or "")
        if not match:
            continue
        n_samples = int(match.group(1))
        n_proteins = int(match.group(2))
        if n_samples > min_samples and n_proteins > 0:
            print(f"{row['cohort']}\t{row['panel']}\t{n_samples}\t{n_proteins}")
PY
    )
    # A failed parse lists nothing, which would otherwise end as a clean exit 0.
    if [ "$n_listed" -eq 0 ]; then
        echo "ERROR: no GWAS panel listed from $QC_SUMMARY"
        FAILED=1
    fi
}

run_panels_from_manifest() {
    if [ ! -f "$PANEL_MANIFEST" ]; then
        echo "missing panel manifest: $PANEL_MANIFEST"
        FAILED=1
        return
    fi

    echo "Panel source: $PANEL_MANIFEST"
    while IFS=$'\t' read -r cohort tag _rest; do
        [ -z "$cohort" ] && continue
        [ "${cohort#\#}" != "$cohort" ] && continue
        [ "$cohort" = "cohort" ] && continue
        case "$cohort" in
            bf1) run_panel "$cohort" "$tag" "$PFILE_BF1" ;;
            bf2) run_panel "$cohort" "$tag" "$PFILE_BF2" ;;
            *) echo "invalid cohort in panel manifest: $cohort $tag"; FAILED=1 ;;
        esac
    done < "$PANEL_MANIFEST"
}

if [ -n "$PANEL_MANIFEST" ]; then
    run_panels_from_manifest
else
    run_panels_from_summary
fi

exit "$FAILED"
