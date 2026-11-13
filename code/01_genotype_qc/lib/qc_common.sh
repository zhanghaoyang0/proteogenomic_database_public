#!/usr/bin/env bash
# Common QC functions and utilities

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "ERROR: required command not found: $1" >&2
    exit 1
  }
}

setup_rscript() {
  # Resolve to an absolute Rscript binary. Never use `conda run` — it silently
  # swallows heredoc stdin and returns 0, masking total R failure.
  local candidates=(
    "${RSCRIPT_BIN}"
    "$HOME/mambaforge/envs/gwas/bin/Rscript"
    "$HOME/miniforge3/envs/gwas/bin/Rscript"
    "$(command -v Rscript 2>/dev/null || true)"
  )
  local c
  for c in "${candidates[@]}"; do
    if [[ -n "${c}" && -x "${c}" ]]; then
      RSCRIPT_CMD=("${c}")
      return 0
    fi
  done
  echo "ERROR: Rscript not found. Set RSCRIPT_BIN=/abs/path/to/Rscript explicitly." >&2
  exit 1
}

assert_file() {
  # assert_file <path> <human-readable message about what produced it>
  if [[ ! -s "$1" ]]; then
    echo "ERROR: expected non-empty file missing: $1" >&2
    echo "       producer: $2" >&2
    exit 1
  fi
}

ensure_dirs() {
  mkdir -p "$@"
}

strip_iid_zeros_in_file() {
  # In-place strip leading zeros from purely-numeric IIDs in a tab-delimited file
  # with a header containing #IID or IID. TOPMed VCFs zero-pad numeric IIDs (e.g.
  # BF2: "0NNNNNNN") which breaks joins to source key/phenotype tables ("NNNNNNN").
  # No-op for non-numeric IIDs.
  local file="$1"
  [[ -f "${file}" ]] || return 0
  local tmp="${file}.iid_strip.tmp"
  awk 'BEGIN{FS=OFS="\t"}
       NR==1 { for (i=1;i<=NF;i++) if ($i=="#IID"||$i=="IID") iid=i; print; next }
       { if ($iid ~ /^[0-9]+$/) { sub(/^0+/, "", $iid); if ($iid=="") $iid="0" } print }' \
       "${file}" > "${tmp}"
  mv "${tmp}" "${file}"
}
