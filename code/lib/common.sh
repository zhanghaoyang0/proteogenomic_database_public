# Shared bash helpers for code. Source with
#   source "<repo>/code/lib/common.sh"

# Same values as GWAS_BASE, CLUMP_SUBDIR and LD_REF in common.R.
GWAS_BASE="$HOME/db/shared_data/gwas"
LD_REF="$HOME/db/shared_data/ref_geno/1kg_v3/EUR"
CLUMP_SUBDIR="output_5e-8Pthres_1000kb_0.01r2"

# fm_clean_file <clean_dir> <trait>: the clean GWAS a fine-mapping trait reads. Mirrors
# the resolver in finemap/01_finemapping.R: flat <trait>.tsv.gz first, then the
# per-tissue <tissue>/<ENSG>.tsv.gz form of eQTL traits ("<tissue>_<ENSG...>").
fm_clean_file() {
  local dir=$1 trait=$2 f tissue
  f="$dir/$trait.tsv.gz"
  if [[ ! -e "$f" && "$trait" == *_ENSG* ]]; then
    tissue=${trait%%_ENSG*}
    f="$dir/$tissue/${trait#"${tissue}_"}.tsv.gz"
  fi
  printf '%s\n' "$f"
}

# fm_is_current <out> <clumped> <clean>: succeeds when <out> exists and is not older
# than either input, so the trait can be skipped. Skip-if-exists alone kept stale fits
# after a re-clump or a re-clean (the 2026-10-05 MetaBrain June fits). A missing clean
# file counts as stale, so the fine-mapping run is the one that fails loudly.
# FORCE=1 redoes every trait -- use it after a change to 01_finemapping.R itself.
# A Python copy of this rule is in finemap/submit_01_finemapping_gwwide.sh; change both.
# Whole seconds, not `-nt`: `-nt` compares to the nanosecond, and 9,929 GTEx .clumped
# files stamped in the same second as their tables but 45-60 ns later counted as stale
# (code review 2026-11-01, C92).
fm_is_current() {
  local out=$1 clumped=$2 clean=$3
  if [[ "${FORCE:-0}" == 1 || ! -e "$out" || ! -e "$clean" ]]; then return 1; fi
  local t; t=$(stat -c %Y "$out")
  [[ $(stat -c %Y "$clumped") -le $t && $(stat -c %Y "$clean") -le $t ]]
}
