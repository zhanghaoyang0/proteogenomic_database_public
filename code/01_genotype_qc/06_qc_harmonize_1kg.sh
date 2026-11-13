#!/usr/bin/env bash
# Harmonize TOPMed-imputed data to a 1KG/hg19 PRS/MR target.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/lib/qc_common.sh"

GENO_DIR="${GENO_DIR:-$HOME/db/secure_data/bf/array/bf2}"
IMPUTED_DIR="${IMPUTED_DIR:-${GENO_DIR}/imputed_qc/raw/01_topmed_qc}"
INPUT_PREFIX="${INPUT_PREFIX:-${IMPUTED_DIR}/raw_imputed}"
INPUT_BUILD="${INPUT_BUILD:-hg38}"
OUT_DIR="${OUT_DIR:-${GENO_DIR}/imputed_qc}"
DATASET="${DATASET:-$(basename "${GENO_DIR}" | sed 's/_geno$//')}"
OUT_PREFIX="${OUT_PREFIX:-${OUT_DIR}/${DATASET}_hg19}"

REF_BIM_1KG="${REF_BIM_1KG:-$HOME/db/shared_data/ref_geno/1kg_v3/EUR.bim}"
HIGH_LD_RANGE="${HIGH_LD_RANGE:-$HOME/db/shared_data/ref_geno/1kg_v3/highld_hg19.range}"
CHAIN_HG38_TO_HG19="${CHAIN_HG38_TO_HG19:-$HOME/db/shared_data/ref_geno/liftover/hg38ToHg19.over.chain.gz}"
LIFTOVER="${LIFTOVER:-$HOME/soft/liftOver}"
KEY_CSV="${KEY_CSV:-$HOME/db/secure_data/bf/pheno/clean/${DATASET}_sid_iid.csv}"
MAKE_KEY_SCRIPT="${MAKE_KEY_SCRIPT:-${SCRIPT_DIR}/03_make_genetic_key.R}"
SEX_EXCLUDE="${SEX_EXCLUDE:-$HOME/db/secure_data/bf/array/sex_check/${DATASET}_sex_exclude.txt}"
NO_PHENO_EXCLUDE="${NO_PHENO_EXCLUDE:-$HOME/db/secure_data/bf/pheno/clean/${DATASET}_iid_no_pheno.txt}"

THREADS="${THREADS:-4}"
MEM_MB="${MEM_MB:-12000}"
PAL_MAF_MAX="${PAL_MAF_MAX:-0.4}"
KEEP_STRAND_FLIP="${KEEP_STRAND_FLIP:-1}"
RUN_MAKE_KEY="${RUN_MAKE_KEY:-0}"
RUN_KEY_FILTER="${RUN_KEY_FILTER:-1}"
RUN_SAMPLE_QC="${RUN_SAMPLE_QC:-1}"
GENO_MAX="${GENO_MAX:-0.02}"
MIND_MAX="${MIND_MAX:-0.02}"
SAMPLE_QC_MAF_MIN="${SAMPLE_QC_MAF_MIN:-0.01}"
OUTPUT_MAF_MIN="${OUTPUT_MAF_MIN:-0.01}"
HWE_MIN="${HWE_MIN:-1e-6}"
HET_SD="${HET_SD:-3}"
KING_CUTOFF="${KING_CUTOFF:-0.2}"
RUN_PCA="${RUN_PCA:-1}"
PCA_COUNT="${PCA_COUNT:-10}"
RUN_APOE="${RUN_APOE:-0}"
APOE_OUT="${APOE_OUT:-$HOME/db/secure_data/bf/pheno/clean/apoe}"
RSCRIPT_BIN="${RSCRIPT_BIN:-}"

RAW_DIR="${OUT_DIR}/raw"
# Scratch, cleaned by 05_cleanup.sh (which keeps only log/). Deliberately not
# hidden: a dot-prefixed dir is skipped by `ls` and `du -sh *`, which is how the
# previous .qc_tmp/ reached 142 GB without anyone noticing.
QC_TMP_DIR="${OUT_DIR}/tmp/02_qc_harmonize"
TMP_DIR="${QC_TMP_DIR}/tmp"
SAMPLE_QC_DIR="${QC_TMP_DIR}/sample_qc"
RAW_REPORT_DIR="${OUT_DIR}/report"
REPORT_DIR="${OUT_DIR}/report"
LOG_DIR="${QC_TMP_DIR}/log"
RAW_HARMONIZED_PREFIX="${SAMPLE_QC_DIR}/hg19_raw"
KEYED_PREFIX="${SAMPLE_QC_DIR}/hg19_key0"
PCA_PRUNE_PREFIX="${RAW_DIR}/02_to_1kg/hg19_prune"

mkdir -p "${TMP_DIR}" "${SAMPLE_QC_DIR}" "${RAW_REPORT_DIR}" "${REPORT_DIR}" "${LOG_DIR}"


check_inputs() {
  need_cmd awk
  need_cmd plink2
  setup_rscript

  case "${INPUT_BUILD}" in
    hg19|GRCh37|grch37|b37) INPUT_BUILD_NORM="hg19" ;;
    hg38|GRCh38|grch38|b38) INPUT_BUILD_NORM="hg38" ;;
    *)
      echo "ERROR: INPUT_BUILD must be hg38 or hg19, got: ${INPUT_BUILD}" >&2
      exit 1
      ;;
  esac

  if [[ "${INPUT_BUILD_NORM}" = "hg19" && "${INPUT_PREFIX}" == */raw/01_topmed_qc/raw_imputed ]]; then
    echo "ERROR: raw TOPMed post-imputation input is expected to be hg38; do not run ${INPUT_PREFIX} with INPUT_BUILD=hg19" >&2
    exit 1
  fi

  if [[ "${INPUT_BUILD_NORM}" = "hg38" ]]; then
    [[ -x "${LIFTOVER}" ]] || {
      echo "ERROR: liftOver executable not found or not executable: ${LIFTOVER}" >&2
      exit 1
    }
    [[ -s "${CHAIN_HG38_TO_HG19}" ]] || {
      echo "ERROR: missing hg38ToHg19 chain: ${CHAIN_HG38_TO_HG19}" >&2
      exit 1
    }
  fi

  for ext in pgen pvar psam; do
    [[ -s "${INPUT_PREFIX}.${ext}" ]] || {
      echo "ERROR: missing ${INPUT_PREFIX}.${ext}" >&2
      exit 1
    }
  done

  [[ -s "${REF_BIM_1KG}" ]] || {
    echo "ERROR: missing 1KG hg19 BIM: ${REF_BIM_1KG}" >&2
    exit 1
  }

  [[ -s "${HIGH_LD_RANGE}" ]] || {
    echo "ERROR: missing high-LD range file: ${HIGH_LD_RANGE}" >&2
    exit 1
  }

  if [[ "${RUN_MAKE_KEY}" = "1" && ! -s "${MAKE_KEY_SCRIPT}" ]]; then
    echo "ERROR: missing make-key script: ${MAKE_KEY_SCRIPT}" >&2
    exit 1
  fi
}

liftover_to_hg19() {
  if [[ "${INPUT_BUILD_NORM}" = "hg19" ]]; then
    echo "==> Step 1: input is already hg19; standardize autosomal pfile..."
    plink2 \
      --pfile "${INPUT_PREFIX}" \
      --chr 1-22 \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${TMP_DIR}/topmed_lifted_hg19" \
      > "${LOG_DIR}/02_standardize_hg19_input.log" 2>&1
    return 0
  fi

  echo "==> Step 1: lift hg38 positions to hg19..."

  awk 'BEGIN{OFS="\t"} !/^#/ {
    c=$1; sub(/^chr/,"",c);
    print "chr" c, $2 - 1, $2, $3
  }' "${INPUT_PREFIX}.pvar" > "${TMP_DIR}/topmed_hg38.bed"

  "${LIFTOVER}" \
    "${TMP_DIR}/topmed_hg38.bed" \
    "${CHAIN_HG38_TO_HG19}" \
    "${TMP_DIR}/topmed_hg19.bed" \
    "${TMP_DIR}/topmed_hg19_unmapped.bed" \
    > "${LOG_DIR}/02_liftover.log" 2>&1

  awk 'BEGIN{OFS="\t"} {
    c=$1; sub(/^chr/,"",c);
    print $4, c, $3
  }' "${TMP_DIR}/topmed_hg19.bed" > "${TMP_DIR}/topmed_hg19.update_chr_pos"

  awk 'BEGIN{OFS="\t"} $2 ~ /^([1-9]|1[0-9]|2[0-2])$/ {print}' \
    "${TMP_DIR}/topmed_hg19.update_chr_pos" \
    > "${TMP_DIR}/topmed_hg19.update_chr_pos.autosome"

  awk 'BEGIN{OFS="\t"} {print $1}' \
    "${TMP_DIR}/topmed_hg19.update_chr_pos.autosome" \
    > "${TMP_DIR}/topmed_hg19.lifted_autosome.ids"

  plink2 \
    --pfile "${INPUT_PREFIX}" \
    --extract "${TMP_DIR}/topmed_hg19.lifted_autosome.ids" \
    --update-chr "${TMP_DIR}/topmed_hg19.update_chr_pos.autosome" 2 1 \
    --update-map "${TMP_DIR}/topmed_hg19.update_chr_pos.autosome" 3 1 \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${TMP_DIR}/topmed_lifted_hg19" \
    > "${LOG_DIR}/02_plink_liftover_update.log" 2>&1
}

build_manifest() {
  echo "==> Step 2: compare lifted BF2 alleles with 1KG/hg19 and build manifest..."

  "${RSCRIPT_CMD[@]}" - \
    "${TMP_DIR}/topmed_lifted_hg19.pvar" \
    "${REF_BIM_1KG}" \
    "${REPORT_DIR}" \
    "${RAW_REPORT_DIR}" \
    "${PAL_MAF_MAX}" \
    "${KEEP_STRAND_FLIP}" <<'RSCRIPT'
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
pvar_file <- args[[1]]
ref_bim <- args[[2]]
report_dir <- args[[3]]
raw_report_dir <- args[[4]]
pal_maf_max <- as.numeric(args[[5]])
keep_strand <- as.logical(as.integer(args[[6]]))

comp <- function(x) chartr("ACGT", "TGCA", toupper(x))
canon_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
allele_pair <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "/")

tgt <- fread(cmd = paste("grep -v '^##'", shQuote(pvar_file)))
setnames(tgt, c("#CHROM", "POS", "ID", "REF", "ALT"), c("CHR", "POS", "TOPMED_ID", "TOPMED_REF", "TOPMED_ALT"))
tgt[, TOPMED_MAF := suppressWarnings(as.numeric(sub(".*MAF=([^;]+).*", "\\1", INFO)))]
tgt <- tgt[, .(
  CHR = canon_chr(CHR),
  POS = as.integer(POS),
  TOPMED_ID = as.character(TOPMED_ID),
  TOPMED_REF = toupper(TOPMED_REF),
  TOPMED_ALT = toupper(TOPMED_ALT),
  TOPMED_MAF = TOPMED_MAF
)]
tgt[, target_key_hg19 := paste(CHR, POS, TOPMED_REF, TOPMED_ALT, sep = ":")]
tgt[, topmed_allele_pair := allele_pair(TOPMED_REF, TOPMED_ALT)]
tgt[, palindromic := topmed_allele_pair %in% c("A/T", "C/G")]
# Drop only when MAF is known and high; unknown MAF → keep (conservative for ADNI / odd INFO)
tgt[, palindromic_high_maf := palindromic & !is.na(TOPMED_MAF) & TOPMED_MAF > pal_maf_max]

ref <- fread(ref_bim, col.names = c("CHR", "RSID", "CM", "POS", "KG_A1", "KG_A2"))
ref <- ref[, .(
  CHR = canon_chr(CHR),
  POS = as.integer(POS),
  RSID = as.character(RSID),
  KG_A1 = toupper(KG_A1),
  KG_A2 = toupper(KG_A2)
)]
ref <- ref[RSID != "." & KG_A1 %in% c("A", "C", "G", "T") & KG_A2 %in% c("A", "C", "G", "T")]
ref[, kg_allele_pair := allele_pair(KG_A1, KG_A2)]
ref[, ref_pos_n := .N, by = .(CHR, POS)]

m <- merge(tgt, ref, by = c("CHR", "POS"), all.x = TRUE, allow.cartesian = TRUE)

m[, allele_status := fcase(
  is.na(RSID), "not_in_1kg",
  TOPMED_REF == KG_A1 & TOPMED_ALT == KG_A2, "same",
  TOPMED_REF == KG_A2 & TOPMED_ALT == KG_A1, "swapped",
  comp(TOPMED_REF) == KG_A1 & comp(TOPMED_ALT) == KG_A2, "strand_flip_same",
  comp(TOPMED_REF) == KG_A2 & comp(TOPMED_ALT) == KG_A1, "strand_flip_swapped",
  default = "mismatch"
)]

m[, duplicate_topmed_id := duplicated(TOPMED_ID) | duplicated(TOPMED_ID, fromLast = TRUE)]
m[, duplicate_rsid := !is.na(RSID) & (duplicated(RSID) | duplicated(RSID, fromLast = TRUE))]

safe_status <- c("same", "swapped")
if (keep_strand) safe_status <- c(safe_status, "strand_flip_same", "strand_flip_swapped")

m[, keep_for_prs_mr := allele_status %in% safe_status | allele_status == "not_in_1kg"]
m[palindromic_high_maf == TRUE, keep_for_prs_mr := FALSE]
m[duplicate_topmed_id == TRUE, keep_for_prs_mr := FALSE]
m[duplicate_rsid == TRUE, keep_for_prs_mr := FALSE]
m[ref_pos_n > 1, keep_for_prs_mr := FALSE]

m[, drop_reason := fcase(
  keep_for_prs_mr, "keep",
  is.na(RSID), "not_in_1kg",
  palindromic_high_maf, "palindromic_high_maf",
  duplicate_topmed_id, "duplicate_topmed_id",
  duplicate_rsid, "duplicate_rsid",
  ref_pos_n > 1, "multiple_1kg_records_at_position",
  default = allele_status
)]

# Final REF/ALT after strand-flip + REF-align for safe variants
aligned <- c("same", "swapped", "strand_flip_same", "strand_flip_swapped")
m[, new_REF := fifelse(allele_status %in% aligned, KG_A1, TOPMED_REF)]
m[, new_ALT := fifelse(allele_status %in% aligned, KG_A2, TOPMED_ALT)]

setcolorder(m, c(
  "TOPMED_ID", "RSID", "CHR", "POS",
  "TOPMED_REF", "TOPMED_ALT", "KG_A1", "KG_A2", "TOPMED_MAF",
  "new_REF", "new_ALT",
  "target_key_hg19", "topmed_allele_pair", "kg_allele_pair",
  "allele_status", "palindromic", "palindromic_high_maf",
  "duplicate_topmed_id", "duplicate_rsid",
  "ref_pos_n", "keep_for_prs_mr", "drop_reason"
))

fwrite(m, file.path(raw_report_dir, "02_variant_manifest.tsv"), sep = "\t")
fwrite(m[keep_for_prs_mr == TRUE, .(TOPMED_ID)], file.path(raw_report_dir, "02_keep_topmed_ids.txt"), col.names = FALSE)
fwrite(m[keep_for_prs_mr == TRUE & !is.na(RSID) & duplicate_rsid == FALSE, .(target_key_hg19, RSID)], file.path(raw_report_dir, "02_update_name_to_1kg_rsid.txt"), sep = "\t", col.names = FALSE)
fwrite(m[keep_for_prs_mr == FALSE], file.path(raw_report_dir, "02_no_safe_1kg_rsid_mapping.tsv"), sep = "\t")

# plink2 --update-alleles file: strand-flip BF2 letters to 1KG strand (both flip_same and flip_swapped)
flip_rows <- unique(
  m[keep_for_prs_mr == TRUE & allele_status %in% c("strand_flip_same", "strand_flip_swapped")],
  by = "target_key_hg19"
)
fwrite(
  flip_rows[, .(target_key_hg19,
                old_A1 = TOPMED_REF, old_A2 = TOPMED_ALT,
                new_A1 = comp(TOPMED_REF), new_A2 = comp(TOPMED_ALT))],
  file.path(raw_report_dir, "02_update_alleles.txt"), sep = "\t", col.names = FALSE
)

# plink2 --ref-allele file: set REF=KG_A1 for swapped + strand_flip_swapped
ref_rows <- unique(
  m[keep_for_prs_mr == TRUE & allele_status %in% c("swapped", "strand_flip_swapped")],
  by = "target_key_hg19"
)
fwrite(
  ref_rows[, .(target_key_hg19, new_ref = KG_A1)],
  file.path(raw_report_dir, "02_ref_allele.txt"), sep = "\t", col.names = FALSE
)

summary <- m[, .N, by = .(allele_status, drop_reason)]
setorder(summary, allele_status, drop_reason)
fwrite(summary, file.path(report_dir, "02_harmonization_summary.tsv"), sep = "\t")

message("Variants with safe 1KG rsID mapping: ", nrow(m[keep_for_prs_mr == TRUE]))
message("Variants without safe 1KG rsID mapping, retained with original IDs downstream: ", nrow(m[keep_for_prs_mr == FALSE]))
RSCRIPT

  # Fail loudly if R silently did nothing (e.g. conda run swallowed heredoc)
  assert_file "${RAW_REPORT_DIR}/02_variant_manifest.tsv" "build_manifest R script"
  assert_file "${RAW_REPORT_DIR}/02_update_name_to_1kg_rsid.txt" "build_manifest R script"
  # 02_ref_allele.txt / 02_update_alleles.txt may legitimately be empty if no
  # swapped / strand-flip variants exist, so we only check the manifest pair.
}

write_harmonized_target() {
  echo "==> Step 3a: rename to chr:hg19pos:ref:alt and dedupe..."

  awk 'BEGIN {OFS="\t"}
    !/^#/ {
      chr = $1; pos = $2; old_id = $3; ref = $4; alt = $5
      new_id = chr ":" pos ":" ref ":" alt
      print old_id, new_id
    }' "${TMP_DIR}/topmed_lifted_hg19.pvar" > "${TMP_DIR}/topmed_hg19_id_update.txt"

  plink2 \
    --pfile "${TMP_DIR}/topmed_lifted_hg19" \
    --update-name "${TMP_DIR}/topmed_hg19_id_update.txt" \
    --rm-dup force-first \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${TMP_DIR}/hg19_renamed" \
    > "${LOG_DIR}/02_rename_dedupe.log" 2>&1

  echo "==> Step 3b: strand-flip (--update-alleles) then REF-align (--ref-allele) in separate plink2 runs..."
  local AFTER_FLIP
  AFTER_FLIP="${TMP_DIR}/hg19_renamed"
  if [[ -s "${RAW_REPORT_DIR}/02_update_alleles.txt" ]]; then
    plink2 \
      --pfile "${TMP_DIR}/hg19_renamed" \
      --update-alleles "${RAW_REPORT_DIR}/02_update_alleles.txt" \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${TMP_DIR}/hg19_strand_flipped" \
      > "${LOG_DIR}/02_strand_flip_only.log" 2>&1
    AFTER_FLIP="${TMP_DIR}/hg19_strand_flipped"
  fi

  if [[ -s "${RAW_REPORT_DIR}/02_ref_allele.txt" ]]; then
    plink2 \
      --pfile "${AFTER_FLIP}" \
      --ref-allele force "${RAW_REPORT_DIR}/02_ref_allele.txt" 2 1 \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${RAW_HARMONIZED_PREFIX}" \
      > "${LOG_DIR}/02_ref_allele_align.log" 2>&1
  else
    plink2 \
      --pfile "${AFTER_FLIP}" \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${RAW_HARMONIZED_PREFIX}" \
      > "${LOG_DIR}/02_ref_allele_align.log" 2>&1
  fi
}

run_key_filter() {
  if [[ "${RUN_MAKE_KEY}" = "1" ]]; then
    echo "==> Step 4: build/refresh BF2 genetic key..."
    "${RSCRIPT_CMD[@]}" "${MAKE_KEY_SCRIPT}" \
      > "${LOG_DIR}/02_make_bf2_key.log" 2>&1
  fi

  if [[ "${RUN_KEY_FILTER}" != "1" ]]; then
    echo "==> Step 4: key filter disabled; copying raw harmonized pfile..."
    plink2 \
      --pfile "${RAW_HARMONIZED_PREFIX}" \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${KEYED_PREFIX}" \
      > "${LOG_DIR}/02_copy_raw_harmonized_to_key0.log" 2>&1
    return 0
  fi

  echo "==> Step 4: key filter BF2 samples..."
  # No fallback key: an old IID list would silently re-apply an outdated sample set.
  if [[ ! -s "${KEY_CSV}" ]]; then
    echo "ERROR: missing BF2 genetic key: ${KEY_CSV} (build it with ${MAKE_KEY_SCRIPT})" >&2
    exit 1
  fi
  # The key is sid,FID,IID: take the IID column by its header, not by position.
  awk -F',' '
    NR == 1 {
      for (i = 1; i <= NF; i++) { gsub(/\r/, "", $i); if ($i == "IID") col = i }
      if (!col) { print "ERROR: no IID column in key header" > "/dev/stderr"; exit 1 }
      next
    }
    {
      gsub(/\r/, "", $col)
      if ($col != "") {
        print $col
        x = $col
        sub(/^0+/, "", x)
        if (x != "") print x
      }
    }' "${KEY_CSV}" | sort -u > "${SAMPLE_QC_DIR}/bf2_key_iid.txt"
  echo "    key: ${KEY_CSV} ($(($(wc -l < "${KEY_CSV}") - 1)) rows)"

  # Genetic vs reported sex mismatches (05_array_sex_check.sh) leave here, so no
  # downstream file carries them. An empty list means none; a missing one is an error.
  # Samples whose sid is not in the phenotype table (03_make_genetic_key.R) leave
  # here too: they have no phenotype, and kept they only compete in het/KING.
  if [[ ! -f "${SEX_EXCLUDE}" ]]; then
    echo "ERROR: missing ${SEX_EXCLUDE} (run ${SCRIPT_DIR}/05_array_sex_check.sh)" >&2
    exit 1
  fi
  if [[ ! -f "${NO_PHENO_EXCLUDE}" ]]; then
    echo "ERROR: missing ${NO_PHENO_EXCLUDE} (run ${MAKE_KEY_SCRIPT})" >&2
    exit 1
  fi
  # FILENAME, not NR == FNR: with both lists empty, NR == FNR would read the key as
  # the drop list.
  awk 'FILENAME == ARGV[1] { x = $1; sub(/^0+/, "", x); if (x != "") drop[x] = 1; next }
       { y = $1; sub(/^0+/, "", y); if (!(y in drop)) print }' \
    <(cat "${SEX_EXCLUDE}" "${NO_PHENO_EXCLUDE}") "${SAMPLE_QC_DIR}/bf2_key_iid.txt" > "${SAMPLE_QC_DIR}/bf2_key_iid.nosex.txt"
  mv "${SAMPLE_QC_DIR}/bf2_key_iid.nosex.txt" "${SAMPLE_QC_DIR}/bf2_key_iid.txt"
  echo "    sex-check exclusions: ${SEX_EXCLUDE} ($(grep -c . "${SEX_EXCLUDE}" || true) IIDs)"
  echo "    no phenotype-table sid: ${NO_PHENO_EXCLUDE} ($(grep -c . "${NO_PHENO_EXCLUDE}" || true) IIDs)"

  awk '
    FILENAME == ARGV[1] {keep[$1] = 1; next}
    FNR == 1 {
      fid_col = 0
      iid_col = 0
      for (i = 1; i <= NF; i++) {
        if ($i == "#FID" || $i == "FID") fid_col = i
        if ($i == "#IID" || $i == "IID") iid_col = i
      }
      if (iid_col == 0) iid_col = (NF >= 2 ? 2 : 1)
      next
    }
    {
      iid = $iid_col
      iid_norm = iid
      sub(/^0+/, "", iid_norm)
      if (iid in keep || iid_norm in keep) {
        fid = (fid_col > 0 ? $fid_col : 0)
        print fid, iid
      }
    }' \
    "${SAMPLE_QC_DIR}/bf2_key_iid.txt" \
    "${RAW_HARMONIZED_PREFIX}.psam" \
    > "${SAMPLE_QC_DIR}/bf2_key_keep.txt"

  if [[ ! -s "${SAMPLE_QC_DIR}/bf2_key_keep.txt" ]]; then
    echo "ERROR: key filter matched zero samples. Check VCF/psam IID format vs ${KEY_CSV}." >&2
    exit 1
  fi
  # Every imputed sample comes from the raw array, so every one is in the key: only
  # the exclusion lists may take samples out here. More dropped than they hold means
  # IIDs that do not match the key (e.g. changed at TOPMed or on import).
  local n_psam n_kept n_excl
  n_psam=$(awk 'NR > 1 {n++} END {print n + 0}' "${RAW_HARMONIZED_PREFIX}.psam")
  n_kept=$(wc -l < "${SAMPLE_QC_DIR}/bf2_key_keep.txt" | tr -d ' ')
  n_excl=$(cat "${SEX_EXCLUDE}" "${NO_PHENO_EXCLUDE}" | grep -c . || true)
  echo "    imputed samples ${n_psam}, kept ${n_kept}, exclusion lists ${n_excl} IIDs"
  if (( n_psam - n_kept > n_excl )); then
    echo "ERROR: $((n_psam - n_kept - n_excl)) or more imputed samples match neither the key nor an exclusion list. Check IID format vs ${KEY_CSV}." >&2
    exit 1
  fi

  plink2 \
    --pfile "${RAW_HARMONIZED_PREFIX}" \
    --keep "${SAMPLE_QC_DIR}/bf2_key_keep.txt" \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${KEYED_PREFIX}" \
    > "${LOG_DIR}/02_sample_qc_step0_key_filter.log" 2>&1
}

run_sample_qc() {
  if [[ "${RUN_SAMPLE_QC}" != "1" ]]; then
    echo "==> Step 5: sample QC disabled; copying key-filtered pfile to final prefix..."
    plink2 \
      --pfile "${KEYED_PREFIX}" \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${OUT_PREFIX}" \
      > "${LOG_DIR}/02_copy_key0_to_final.log" 2>&1
    return 0
  fi

  echo "==> Step 5: run BF2 sample QC on harmonized 1KG/hg19 target..."

  plink2 \
    --pfile "${KEYED_PREFIX}" \
    --missing \
    --out "${SAMPLE_QC_DIR}/hg19_step0_miss" \
    > "${LOG_DIR}/02_sample_qc_step0_missing.log" 2>&1

  plink2 \
    --pfile "${KEYED_PREFIX}" \
    --geno "${GENO_MAX}" dosage \
    --mind "${MIND_MAX}" dosage \
    --make-pgen \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${SAMPLE_QC_DIR}/hg19_step1_missingness" \
    > "${LOG_DIR}/02_sample_qc_step1_geno_mind.log" 2>&1

  plink2 \
    --pfile "${SAMPLE_QC_DIR}/hg19_step1_missingness" \
    --autosome \
    --maf "${SAMPLE_QC_MAF_MIN}" \
    --hwe "${HWE_MIN}" \
    --make-pgen \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${SAMPLE_QC_DIR}/hg19_step2_hwe_maf" \
    > "${LOG_DIR}/02_sample_qc_step2_hwe_maf.log" 2>&1

  plink2 \
    --pfile "${SAMPLE_QC_DIR}/hg19_step2_hwe_maf" \
    --exclude range "${HIGH_LD_RANGE}" \
    --indep-pairwise 50 5 0.2 \
    --out "${SAMPLE_QC_DIR}/hg19_step2_prune" \
    > "${LOG_DIR}/02_sample_qc_step3_prune_for_het.log" 2>&1

  plink2 \
    --pfile "${SAMPLE_QC_DIR}/hg19_step2_hwe_maf" \
    --extract "${SAMPLE_QC_DIR}/hg19_step2_prune.prune.in" \
    --het \
    --out "${SAMPLE_QC_DIR}/hg19_step2" \
    > "${LOG_DIR}/02_sample_qc_step3_het.log" 2>&1

  : > "${SAMPLE_QC_DIR}/hg19_step3_het_remove_ind.txt"
  awk -v out="${SAMPLE_QC_DIR}/hg19_step3_het_remove_ind.txt" \
      -v sd_mult="${HET_SD}" '
    NR == 1 {
      fid_col = 0
      iid_col = 0
      f_col = 0
      for (i = 1; i <= NF; i++) {
        if ($i == "#FID" || $i == "FID") fid_col = i
        if ($i == "#IID" || $i == "IID") iid_col = i
        if ($i == "F") f_col = i
      }
      if (iid_col == 0) iid_col = (NF >= 2 ? 2 : 1)
      if (f_col == 0) {
        print "[het] ERROR: F column not found" > "/dev/stderr"
        exit 1
      }
      next
    }
    {
      n++
      fid[n] = (fid_col > 0 ? $fid_col : 0)
      iid[n] = $iid_col
      f[n] = $f_col
      sum += f[n]; sumsq += f[n] * f[n]
    }
    END {
      if (n < 2) {
        print "[het] too few samples for outlier detection; removed: 0" > "/dev/stderr"
        exit
      }
      mean = sum / n
      sd = sqrt((sumsq - sum * sum / n) / (n - 1))
      low = mean - sd_mult * sd
      high = mean + sd_mult * sd
      removed = 0
      for (i = 1; i <= n; i++) {
        if (f[i] < low || f[i] > high) {
          print fid[i], iid[i] >> out
          removed++
        }
      }
      print "[het] removed:", removed, "mean:", mean, "sd:", sd, "low:", low, "high:", high > "/dev/stderr"
    }
  ' "${SAMPLE_QC_DIR}/hg19_step2.het" \
    > "${LOG_DIR}/02_sample_qc_step3_het_outliers.log" 2>&1

  plink2 \
    --pfile "${SAMPLE_QC_DIR}/hg19_step2_hwe_maf" \
    --remove "${SAMPLE_QC_DIR}/hg19_step3_het_remove_ind.txt" \
    --make-pgen \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${SAMPLE_QC_DIR}/hg19_step3_het" \
    > "${LOG_DIR}/02_sample_qc_step3_remove_het.log" 2>&1

  plink2 \
    --pfile "${SAMPLE_QC_DIR}/hg19_step3_het" \
    --exclude range "${HIGH_LD_RANGE}" \
    --indep-pairwise 50 5 0.2 \
    --out "${SAMPLE_QC_DIR}/hg19_step3_prune" \
    > "${LOG_DIR}/02_sample_qc_step4_prune_for_king.log" 2>&1

  plink2 \
    --pfile "${SAMPLE_QC_DIR}/hg19_step3_het" \
    --extract "${SAMPLE_QC_DIR}/hg19_step3_prune.prune.in" \
    --king-cutoff "${KING_CUTOFF}" \
    --out "${SAMPLE_QC_DIR}/hg19_step3_king" \
    > "${LOG_DIR}/02_sample_qc_step4_king.log" 2>&1

  : > "${SAMPLE_QC_DIR}/hg19_remove_samples.txt"
  for remove_file in \
    "${SAMPLE_QC_DIR}/hg19_step3_het_remove_ind.txt" \
    "${SAMPLE_QC_DIR}/hg19_step3_king.king.cutoff.out.id"; do
    if [[ -s "${remove_file}" ]]; then
      awk '
        NR == 1 && ($1 == "#IID" || $1 == "IID" || $1 == "#FID" || $1 == "FID") { next }
        NF == 1 { print 0, $1; next }
        NF >= 2 { print $1, $2 }
      ' "${remove_file}" >> "${SAMPLE_QC_DIR}/hg19_remove_samples.txt"
    fi
  done
  sort -u "${SAMPLE_QC_DIR}/hg19_remove_samples.txt" \
    -o "${SAMPLE_QC_DIR}/hg19_remove_samples.txt"

  if [[ -s "${SAMPLE_QC_DIR}/hg19_remove_samples.txt" ]]; then
    plink2 \
      --pfile "${SAMPLE_QC_DIR}/hg19_step1_missingness" \
      --remove "${SAMPLE_QC_DIR}/hg19_remove_samples.txt" \
      --maf "${OUTPUT_MAF_MIN}" \
      --hwe "${HWE_MIN}" \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${OUT_PREFIX}" \
      > "${LOG_DIR}/02_sample_qc_step4_remove_related_final.log" 2>&1
  else
    plink2 \
      --pfile "${SAMPLE_QC_DIR}/hg19_step1_missingness" \
      --maf "${OUTPUT_MAF_MIN}" \
      --hwe "${HWE_MIN}" \
      --make-pgen \
      --sort-vars \
      --threads "${THREADS}" \
      --memory "${MEM_MB}" \
      --out "${OUT_PREFIX}" \
      > "${LOG_DIR}/02_sample_qc_step4_remove_related_final.log" 2>&1
  fi
}

run_pca() {
  if [[ "${RUN_PCA}" != "1" ]]; then
    return 0
  fi

  mkdir -p "$(dirname "${PCA_PRUNE_PREFIX}")"
  echo "==> Step 6: compute LD-pruned PCs..."
  plink2 \
    --pfile "${OUT_PREFIX}" \
    --indep-pairwise 200 50 0.2 \
    --exclude range "${HIGH_LD_RANGE}" \
    --out "${PCA_PRUNE_PREFIX}" \
    > "${LOG_DIR}/02_pca_prune.log" 2>&1

  plink2 \
    --pfile "${OUT_PREFIX}" \
    --extract "${PCA_PRUNE_PREFIX}.prune.in" \
    --pca "${PCA_COUNT}" \
    --out "${OUT_PREFIX}_pca" \
    > "${LOG_DIR}/02_pca.log" 2>&1
}

export_apoe() {
  if [[ "${RUN_APOE}" != "1" ]]; then
    return 0
  fi

  echo "==> Step 9: export APOE SNPs..."
  mkdir -p "${APOE_OUT}"
  plink2 \
    --pfile "${OUT_PREFIX}" \
    --snps rs429358,rs7412 \
    --export A ref-first \
    --freq \
    --out "${APOE_OUT}/hg19_apoe" \
    > "${LOG_DIR}/02_apoe.log" 2>&1
}

assign_rsid_final() {
  # build_manifest is required upstream, so the mapping file must exist by now.
  # Treat a missing file as a hard error rather than silently keeping chr:pos:ref:alt IDs.
  assert_file "${RAW_REPORT_DIR}/02_update_name_to_1kg_rsid.txt" "build_manifest R script (Step 2)"

  echo "==> Step 8: assign rsIDs from 1KG..."
  plink2 \
    --pfile "${OUT_PREFIX}" \
    --update-name "${RAW_REPORT_DIR}/02_update_name_to_1kg_rsid.txt" \
    --make-pgen \
    --sort-vars \
    --threads "${THREADS}" \
    --memory "${MEM_MB}" \
    --out "${TMP_DIR}/harmonized_rsid_tmp" \
    > "${LOG_DIR}/02_assign_rsid_final.log" 2>&1

  mv "${TMP_DIR}/harmonized_rsid_tmp.pgen" "${OUT_PREFIX}.pgen"
  mv "${TMP_DIR}/harmonized_rsid_tmp.pvar" "${OUT_PREFIX}.pvar"
  mv "${TMP_DIR}/harmonized_rsid_tmp.psam" "${OUT_PREFIX}.psam"
}

strip_iid_leading_zeros() {
  # pgen indexes samples by row position, so editing psam in place is safe.
  echo "==> Step 10: strip leading zeros from numeric IIDs..."
  strip_iid_zeros_in_file "${OUT_PREFIX}.psam"
  strip_iid_zeros_in_file "${OUT_PREFIX}_pca.eigenvec"
}

make_sample_report() {
  echo "==> Step 10: write sample and variant count reports..."

  awk 'NR > 1 {print $1, $2}' "${INPUT_PREFIX}.psam" | sort > "${TMP_DIR}/input.samples"
  awk 'NR > 1 {print $1, $2}' "${RAW_HARMONIZED_PREFIX}.psam" | sort > "${TMP_DIR}/raw_harmonized.samples"
  awk 'NR > 1 {print $1, $2}' "${KEYED_PREFIX}.psam" | sort > "${TMP_DIR}/keyed.samples"
  awk 'NR > 1 {print $1, $2}' "${OUT_PREFIX}.psam" | sort > "${TMP_DIR}/harmonized.samples"

  {
    echo -e "metric\tvalue"
    echo -e "input_samples\t$(wc -l < "${TMP_DIR}/input.samples" | tr -d ' ')"
    echo -e "raw_harmonized_samples\t$(wc -l < "${TMP_DIR}/raw_harmonized.samples" | tr -d ' ')"
    echo -e "key_filtered_samples\t$(wc -l < "${TMP_DIR}/keyed.samples" | tr -d ' ')"
    echo -e "final_samples\t$(wc -l < "${TMP_DIR}/harmonized.samples" | tr -d ' ')"
    echo -e "input_final_sample_overlap\t$(comm -12 "${TMP_DIR}/input.samples" "${TMP_DIR}/harmonized.samples" | wc -l | tr -d ' ')"
    echo -e "input_variants\t$(awk '!/^#/ {n++} END{print n + 0}' "${INPUT_PREFIX}.pvar")"
    echo -e "lifted_hg19_variants\t$(awk '!/^#/ {n++} END{print n + 0}' "${TMP_DIR}/topmed_lifted_hg19.pvar")"
    echo -e "raw_harmonized_variants\t$(awk '!/^#/ {n++} END{print n + 0}' "${RAW_HARMONIZED_PREFIX}.pvar")"
    echo -e "key_filtered_variants\t$(awk '!/^#/ {n++} END{print n + 0}' "${KEYED_PREFIX}.pvar")"
    echo -e "final_variants\t$(awk '!/^#/ {n++} END{print n + 0}' "${OUT_PREFIX}.pvar")"
    echo -e "run_make_key\t${RUN_MAKE_KEY}"
    echo -e "run_key_filter\t${RUN_KEY_FILTER}"
    echo -e "run_sample_qc\t${RUN_SAMPLE_QC}"
    echo -e "input_build\t${INPUT_BUILD_NORM}"
    echo -e "run_pca\t${RUN_PCA}"
    echo -e "pca_count\t${PCA_COUNT}"
    echo -e "run_apoe\t${RUN_APOE}"
    echo -e "geno_max\t${GENO_MAX}"
    echo -e "mind_max\t${MIND_MAX}"
    echo -e "sample_qc_maf_min\t${SAMPLE_QC_MAF_MIN}"
    echo -e "output_maf_min\t${OUTPUT_MAF_MIN}"
    echo -e "hwe_min\t${HWE_MIN}"
    echo -e "het_sd\t${HET_SD}"
    echo -e "king_cutoff\t${KING_CUTOFF}"
    echo -e "pal_maf_max\t${PAL_MAF_MAX}"
    echo -e "keep_strand_flip\t${KEEP_STRAND_FLIP}"
  } > "${REPORT_DIR}/02_target_counts.tsv"
}

main() {
  check_inputs
  liftover_to_hg19
  build_manifest
  write_harmonized_target
  run_key_filter
  run_sample_qc
  run_pca
  assign_rsid_final
  strip_iid_leading_zeros
  export_apoe
  make_sample_report

  echo "Done."
  echo "Harmonized pfile: ${OUT_PREFIX}"
  echo "Raw manifest: ${RAW_REPORT_DIR}/02_variant_manifest.tsv"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  main "$@"
fi
