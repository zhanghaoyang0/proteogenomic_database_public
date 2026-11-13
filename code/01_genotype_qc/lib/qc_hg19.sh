#!/usr/bin/env bash

run_preimpute_qc() {
  set -euo pipefail

  : "${DATASET:?Set DATASET in the dataset wrapper}"
  : "${INPUT_TYPE:?Set INPUT_TYPE to bfile or pfile in the dataset wrapper}"
  : "${INPUT_PREFIX:?Set INPUT_PREFIX in the dataset wrapper}"
  : "${OUT_DIR:?Set OUT_DIR in the dataset wrapper}"

  REF_BIM_1KG="${REF_BIM_1KG:-$HOME/db/shared_data/ref_geno/1kg_v3/EUR.bim}"
  HIGH_LD_RANGE="${HIGH_LD_RANGE:-$HOME/db/shared_data/ref_geno/1kg_v3/highld_hg19.range}"
  REF_FASTA_HG19="${REF_FASTA_HG19:-$HOME/db/shared_data/ref_geno/fasta/hg19.fa}"
  HET_SD="${HET_SD:-3}"
  MAX_ALLELE_MISMATCH_RATE="${MAX_ALLELE_MISMATCH_RATE:-0.05}"
  MAX_NO_REFPOS_RATE="${MAX_NO_REFPOS_RATE:-0.40}"
  PREIMPUTE_GENO="${PREIMPUTE_GENO:-0.02}"
  PREIMPUTE_MIND="${PREIMPUTE_MIND:-0.02}"

  PREFIX="${OUT_DIR}/${DATASET}"
  TMP_DEDUP="${OUT_DIR}/${DATASET}_tmp_dedup"
  TMP_QC0="${PREFIX}_qc0"
  TMP_QC1="${PREFIX}_qc1"
  TMP_QC2="${PREFIX}_qc2"
  FINAL="${PREFIX}"
  CHR_RENAME_TO_CHR="${OUT_DIR}/chr_rename_to_chr_hg19.txt"
  CHR_RENAME_TO_NOCHR="${OUT_DIR}/chr_rename_to_nochr_hg19.txt"
  LOG_DIR="${OUT_DIR}/log"
  VCF_DIR="${OUT_DIR}/vcf"

  mkdir -p "${OUT_DIR}" "${LOG_DIR}" "${VCF_DIR}"
  _require_file "${REF_BIM_1KG}" "1KG hg19 BIM"
  _require_file "${HIGH_LD_RANGE}" "high-LD range file"
  _require_file "${REF_FASTA_HG19}" "hg19 FASTA"
  if [ ! -s "${REF_FASTA_HG19}.fai" ]; then
    echo "==> Indexing hg19 FASTA..."
    samtools faidx "${REF_FASTA_HG19}"
  fi
  _write_chr_maps "${CHR_RENAME_TO_CHR}" "${CHR_RENAME_TO_NOCHR}"

  echo "==> [${DATASET}] Step 1: autosome ACGT SNPs + dedup..."
  if [ "${INPUT_TYPE}" = "bfile" ]; then
    plink2 --bfile "${INPUT_PREFIX}" \
           --chr 1-22 \
           --snps-only just-acgt \
           --set-all-var-ids '@:#:$r:$a' \
           --rm-dup force-first \
           --new-id-max-allele-len 500 \
           --make-pgen \
           --out "${TMP_DEDUP}"
  elif [ "${INPUT_TYPE}" = "pfile" ]; then
    plink2 --pfile "${INPUT_PREFIX}" \
           --chr 1-22 \
           --snps-only just-acgt \
           --set-all-var-ids '@:#:$r:$a' \
           --rm-dup force-first \
           --new-id-max-allele-len 500 \
           --make-pgen \
           --out "${TMP_DEDUP}"
  else
    echo "ERROR: INPUT_TYPE must be bfile or pfile, got ${INPUT_TYPE}" >&2
    exit 1
  fi
  cp "${TMP_DEDUP}.log" "${LOG_DIR}/step1_dedup.log"

  # plink2 applies --mind before --geno within one run, so on an array with many
  # failed probes (ADNI1 610-Quad, overall call rate 0.956) every sample fails
  # --mind. Drop high-missingness variants in their own pass first.
  echo "==> [${DATASET}] Step 2a: variant missingness (--geno) before sample missingness..."
  plink2 --pfile "${TMP_DEDUP}" \
         --geno "${PREIMPUTE_GENO}" \
         --make-pgen \
         --out "${TMP_DEDUP}_geno"
  cp "${TMP_DEDUP}_geno.log" "${LOG_DIR}/step2a_geno.log"

  echo "==> [${DATASET}] Step 2: missingness + HWE + low MAF filters..."
  plink2 --pfile "${TMP_DEDUP}_geno" \
         --geno "${PREIMPUTE_GENO}" --mind "${PREIMPUTE_MIND}" \
         --hwe 1e-6 \
         --maf 0.001 \
         --make-pgen \
         --out "${TMP_QC0}"
  cp "${TMP_QC0}.log" "${LOG_DIR}/step2_qc.log"

  echo "==> [${DATASET}] Step 3: hg19 build/strand/allele sanity check against 1KG EUR..."
  _allele_build_check "${TMP_QC0}.pvar" "${PREFIX}_allele_check"

  # combine allele mismatches + palindromic AT/CG SNPs into one exclusion list
  exclude_mismatch="${PREFIX}_allele_check.exclude_allele_mismatch.txt"
  exclude_ambiguous="${PREFIX}_allele_check.ambiguous_AT_CG.txt"
  exclude_combined="${PREFIX}_allele_check.exclude_combined.txt"
  cat "${exclude_mismatch}" "${exclude_ambiguous}" | sort -u | grep -v '^$' > "${exclude_combined}" || true

  if [ -s "${exclude_combined}" ]; then
    n_excl=$(wc -l < "${exclude_combined}")
    echo "==> [${DATASET}] Excluding ${n_excl} variants (allele mismatches + palindromic AT/CG)..."
    plink2 --pfile "${TMP_QC0}" \
           --exclude "${exclude_combined}" \
           --make-pgen \
           --out "${TMP_QC1}"
  else
    echo "==> [${DATASET}] No variants to exclude at step 3."
    plink2 --pfile "${TMP_QC0}" \
           --make-pgen \
           --out "${TMP_QC1}"
  fi
  cp "${TMP_QC1}.log" "${LOG_DIR}/step3_allele_check.log"

  echo "==> [${DATASET}] Step 4: heterozygosity QC on LD-pruned autosomal SNPs..."
  plink2 --pfile "${TMP_QC1}" \
         --exclude range "${HIGH_LD_RANGE}" \
         --indep-pairwise 50 5 0.2 \
         --out "${LOG_DIR}/step4_het_prune"

  plink2 --pfile "${TMP_QC1}" \
         --extract "${LOG_DIR}/step4_het_prune.prune.in" \
         --het \
         --out "${LOG_DIR}/step4_het"

  _write_het_outliers "${LOG_DIR}/step4_het.het" "${PREFIX}_het_remove_ind.txt"

  if [ -s "${PREFIX}_het_remove_ind.txt" ]; then
    echo "==> [${DATASET}] Removing heterozygosity outliers..."
    plink2 --pfile "${TMP_QC1}" \
           --remove "${PREFIX}_het_remove_ind.txt" \
           --make-pgen \
           --out "${TMP_QC2}"
  else
    echo "==> [${DATASET}] No heterozygosity outliers detected."
    plink2 --pfile "${TMP_QC1}" \
           --make-pgen \
           --out "${TMP_QC2}"
  fi
  cp "${TMP_QC2}.log" "${LOG_DIR}/step4_het_remove.log"

  echo "==> [${DATASET}] Step 5: write final pre-imputation pfile..."
  plink2 --pfile "${TMP_QC2}" \
         --make-pgen \
         --out "${FINAL}"
  cp "${FINAL}.log" "${LOG_DIR}/step5_final.log"

  echo "==> [${DATASET}] Step 6: export per-chromosome VCF + hg19 FASTA strand/ref fix..."
  fixref_summary="${PREFIX}_fixref.summary.tsv"
  echo -e "chr\tmetric\tvalue\tpct" > "${fixref_summary}"

  for chr in {1..22}; do
    raw_vcf="${VCF_DIR}/${DATASET}_chr${chr}.raw.vcf.gz"
    renamed_vcf="${VCF_DIR}/${DATASET}_chr${chr}.renamed.vcf.gz"
    fixed_vcf="${VCF_DIR}/${DATASET}_chr${chr}.fixed.vcf.gz"
    final_vcf="${VCF_DIR}/${DATASET}_chr${chr}.vcf.gz"
    fixref_log="${LOG_DIR}/step6_fixref_chr${chr}.log"

    plink2 --pfile "${FINAL}" \
           --chr "${chr}" \
           --export vcf bgz id-paste=iid \
           --out "${VCF_DIR}/${DATASET}_chr${chr}.raw"
    mv "${VCF_DIR}/${DATASET}_chr${chr}.raw.log" "${LOG_DIR}/step6_vcf_chr${chr}.log"

    # Always check and fix chromosome naming if needed. Use the header instead
    # of early-exiting a data stream; with pipefail, gunzip can report SIGPIPE.
    if bcftools view -h "${raw_vcf}" | grep -Eq '^##contig=<ID=([1-9]|1[0-9]|2[0-2])([,>])'; then
      echo "==> [${DATASET}] chr${chr}: renaming numeric to UCSC format..."
      bcftools annotate --rename-chrs <(echo "${chr} chr${chr}") "${raw_vcf}" -Oz -o "${renamed_vcf}"
      tabix -f -p vcf "${renamed_vcf}"
      vcf_for_fixref="${renamed_vcf}"
    else
      vcf_for_fixref="${raw_vcf}"
    fi

    bcftools +fixref "${vcf_for_fixref}" \
      -Oz -o "${fixed_vcf}" \
      -- -f "${REF_FASTA_HG19}" -m flip -d \
      > "${fixref_log}" 2>&1

    # fixref uses the local UCSC-style hg19 FASTA; TOPMed hg19/GRCh37
    # upload expects chromosome names without the chr prefix.
    bcftools annotate --rename-chrs <(echo "chr${chr} ${chr}") "${fixed_vcf}" -Oz -o "${final_vcf}"
    _append_fixref_summary "${chr}" "${fixref_log}" "${fixref_summary}"

    rm -f "${raw_vcf}" "${raw_vcf}.tbi"
    [ -n "${renamed_vcf}" ] && [ -f "${renamed_vcf}" ] && rm -f "${renamed_vcf}" "${renamed_vcf}.tbi"
    rm -f "${fixed_vcf}" "${fixed_vcf}.tbi"
  done

  rm -f "${TMP_DEDUP}".{pgen,pvar,psam,log}
  rm -f "${TMP_DEDUP}_geno".{pgen,pvar,psam,log}
  rm -f "${TMP_QC0}".{pgen,pvar,psam,log}
  rm -f "${TMP_QC1}".{pgen,pvar,psam,log}
  rm -f "${TMP_QC2}".{pgen,pvar,psam,log}
  rm -f "${CHR_RENAME_TO_CHR}" "${CHR_RENAME_TO_NOCHR}"

  echo "==> [${DATASET}] Cleaning up intermediate files..."
  rm -f "${PREFIX}_allele_check.tsv"
  rm -f "${PREFIX}_allele_check.exclude_allele_mismatch.txt"
  rm -f "${PREFIX}_allele_check.exclude_combined.txt"
  rm -f "${PREFIX}_allele_check.ambiguous_AT_CG.txt"
  rm -f "${PREFIX}_het_remove_ind.txt"

  echo ""
  echo "=== [${DATASET}] Done ==="
  echo "Final pfile : ${FINAL}.{pgen,pvar,psam}"
  echo "VCF for upload : ${VCF_DIR}/${DATASET}_chr{1..22}.vcf.gz"
  echo "Step logs   : ${LOG_DIR}/"
  echo "Allele check: ${PREFIX}_allele_check.summary.tsv"
  echo "Het QC      : ${PREFIX}_het_remove_ind.txt.summary.tsv"
  echo "Fixref      : ${PREFIX}_fixref.summary.tsv"
  echo "Upload ${VCF_DIR}/${DATASET}_chr{1..22}.vcf.gz to TOPMed; select input build hg19/GRCh37."
}

_require_file() {
  local path="$1"
  local label="$2"
  if [ ! -s "${path}" ]; then
    echo "ERROR: ${label} not found: ${path}" >&2
    exit 1
  fi
}

_write_chr_maps() {
  local to_chr="$1"
  local to_nochr="$2"
  : > "${to_chr}"
  : > "${to_nochr}"
  for chr in {1..22}; do
    echo -e "${chr}\tchr${chr}" >> "${to_chr}"
    echo -e "chr${chr}\tchr${chr}" >> "${to_chr}"
    echo -e "chr${chr}\t${chr}" >> "${to_nochr}"
    echo -e "${chr}\t${chr}" >> "${to_nochr}"
  done
}

_allele_build_check() {
  local pvar="$1"
  local out_prefix="$2"
  awk -v out_prefix="${out_prefix}" \
      -v max_mismatch="${MAX_ALLELE_MISMATCH_RATE}" \
      -v max_noref="${MAX_NO_REFPOS_RATE}" '
  BEGIN {
    OFS = "\t"
    check_tsv = out_prefix ".tsv"
    summary_tsv = out_prefix ".summary.tsv"
    mismatch_txt = out_prefix ".exclude_allele_mismatch.txt"
    ambiguous_txt = out_prefix ".ambiguous_AT_CG.txt"
    # truncate output files so stale runs do not leave old IDs
    printf "" > mismatch_txt; close(mismatch_txt)
    printf "" > ambiguous_txt; close(ambiguous_txt)
  }
  function comp(a, out, i, b) {
    out = ""
    for (i = 1; i <= length(a); i++) {
      b = substr(a, i, 1)
      out = out (b == "A" ? "T" : b == "T" ? "A" : b == "C" ? "G" : b == "G" ? "C" : b)
    }
    return out
  }
  function unordered_match(a, b, c, d) {
    return ((a == c && b == d) || (a == d && b == c))
  }
  NR == FNR {
    if ($1 ~ /^([1-9]|1[0-9]|2[0-2])$/) {
      key = $1 ":" $4
      a1[key] = toupper($5)
      a2[key] = toupper($6)
    }
    next
  }
  /^#/ { next }
  {
    total++
    chr = $1
    sub(/^chr/, "", chr)
    pos = $2
    id = $3
    ref = toupper($4)
    alt = toupper($5)
    key = chr ":" pos
    if (ref alt ~ /^(AT|TA|CG|GC)$/) {
      ambiguous++
      print id > ambiguous_txt
    }

    if (!(key in a1)) {
      no_refpos++
      status = "no_refpos"
    } else if (unordered_match(ref, alt, a1[key], a2[key])) {
      allele_match++
      status = "allele_match"
    } else if (unordered_match(comp(ref), comp(alt), a1[key], a2[key])) {
      strand_match++
      status = "strand_match"
    } else {
      allele_mismatch++
      status = "allele_mismatch"
      print id > mismatch_txt
    }

    print chr, pos, id, ref, alt, (key in a1 ? a1[key] : "NA"), (key in a1 ? a2[key] : "NA"), status > check_tsv
  }
  END {
    mismatch_rate = (total ? allele_mismatch / total : 1)
    noref_rate = (total ? no_refpos / total : 1)
    print "metric", "value" > summary_tsv
    print "total", total >> summary_tsv
    print "allele_match", allele_match + 0 >> summary_tsv
    print "strand_match", strand_match + 0 >> summary_tsv
    print "allele_mismatch", allele_mismatch + 0 >> summary_tsv
    print "no_refpos", no_refpos + 0 >> summary_tsv
    print "ambiguous_AT_CG", ambiguous + 0 >> summary_tsv
    print "allele_mismatch_rate", mismatch_rate >> summary_tsv
    print "no_refpos_rate", noref_rate >> summary_tsv
    if (mismatch_rate > max_mismatch) {
      print "ERROR: allele mismatch rate " mismatch_rate " > " max_mismatch > "/dev/stderr"
      exit 2
    }
    if (noref_rate > max_noref) {
      print "ERROR: no-ref-position rate " noref_rate " > " max_noref "; possible build mismatch" > "/dev/stderr"
      exit 3
    }
  }
  ' "${REF_BIM_1KG}" "${pvar}"
}

_write_het_outliers() {
  local het_file="$1"
  local out="$2"
  awk -v sd_cutoff="${HET_SD}" -v out="${out}" '
  BEGIN {
    OFS = "\t"
    summary_tsv = out ".summary.tsv"
  }
  # Columns by header: plink2 leaves out FID when every FID is 0, which moves F
  # from column 6 to 5 (a fixed $6 then removed nobody, until 2026-11-04).
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      if ($i == "#FID" || $i == "FID") fid_col = i
      if ($i == "#IID" || $i == "IID") iid_col = i
      if ($i == "F") f_col = i
    }
    if (!iid_col || !f_col) {
      print "ERROR: no IID or F column in " FILENAME > "/dev/stderr"
      bad = 1
      exit 1
    }
    next
  }
  {
    n++
    fid[n] = (fid_col ? $fid_col : 0)
    iid[n] = $iid_col
    f[n] = $f_col
    sum += f[n]
    sumsq += f[n] * f[n]
  }
  END {
    if (bad) exit 1
    if (n < 2) {
      print "WARN: fewer than 2 samples for het QC; no het outliers removed" > "/dev/stderr"
      print "metric", "value" > summary_tsv
      print "n_samples", n >> summary_tsv
      print "removed", 0 >> summary_tsv
      exit 0
    }
    mean = sum / n
    sd = sqrt((sumsq - (sum * sum / n)) / (n - 1))
    low = mean - sd_cutoff * sd
    high = mean + sd_cutoff * sd
    removed = 0
    for (i = 1; i <= n; i++) {
      if (f[i] < low || f[i] > high) {
        print fid[i], iid[i] > out
        removed++
      }
    }
    print "metric", "value" > summary_tsv
    print "n_samples", n >> summary_tsv
    print "mean_F", mean >> summary_tsv
    print "sd_F", sd >> summary_tsv
    print "low_cutoff", low >> summary_tsv
    print "high_cutoff", high >> summary_tsv
    print "removed", removed >> summary_tsv
  }
  ' "${het_file}"
}

_append_fixref_summary() {
  local chr="$1"
  local fixref_log="$2"
  local summary="$3"
  awk -v chr="${chr}" -F "\t" '
    /^NS/ {
      metric = $2
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", metric)
      pct = ($4 == "" ? "" : $4)
      print chr "\t" metric "\t" $3 "\t" pct
    }
  ' "${fixref_log}" >> "${summary}"
}
