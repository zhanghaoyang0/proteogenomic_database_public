#!/usr/bin/env bash
# Genetic vs reported sex for the BF1/BF2 arrays; writes the IIDs that
# 06_qc_harmonize_1kg.sh drops in its key filter.
#
# Genetic sex: raw array chrX (the imputed pfiles are autosomes only), non-PAR hg19,
# biallelic ACGT SNPs, MAF >= 0.05, missing <= 0.05; X het rate < 0.10 = M,
# > 0.20 = F, between = ambiguous. Same rule as process_wgs/03_bf2_qc/03_bf2_apoe.sh.
# Reported sex: PET export, Visit 0, gender_baseline_variable (0 = M, 1 = F),
# joined through bf{1,2}_sid_iid.csv -- rerun this after 03_make_genetic_key.R.
# Excluded: mismatch or ambiguous. No reported sex is not a reason to exclude.
#
# Output (secure): ~/db/secure_data/bf/array/sex_check/
#   bf{1,2}_array_sex_check.tsv   one row per raw array sample
#   bf{1,2}_sex_exclude.txt       IIDs to drop (may be empty, must exist)
# The screen gets counts only.
#   bash 05_array_sex_check.sh
set -euo pipefail

BF=${BF_ROOT:-$HOME/db/secure_data/bf}
PY=${PY:-python3}
PLINK2=${PLINK2:-$(command -v plink2 || echo /opt/homebrew/bin/plink2)}
PET=${PET:-$BF/pheno/raw/haoyang__20260227_112046_PET/haoyang__20260227_112046_x.csv}
OUT=$BF/array/sex_check
[[ -x "$PLINK2" ]] || { echo "ERROR: no plink2 ($PLINK2)" >&2; exit 1; }
for f in "$PET" "$BF/pheno/clean/bf1_sid_iid.csv" "$BF/pheno/clean/bf2_sid_iid.csv"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
w=$(mktemp -d)
trap 'rm -rf "$w"' EXIT
mkdir -p "$OUT"
printf 'X\t60000\t2699520\tPAR1\nX\t154931043\t155260560\tPAR2\n' > "$w/par_hg19.bed"

# SEX=NA so plink2 counts chrX as diploid for everyone (otherwise men's hets are zeroed).
count_x() {  # count_x <cohort> <plink2 input args...>
  local c=$1; shift
  "$PLINK2" "$@" --psam "$w/$c.psam" --allow-extra-chr --chr X --exclude bed0 "$w/par_hg19.bed" \
    --snps-only just-acgt --max-alleles 2 --maf 0.05 --geno 0.05 \
    --sample-counts --out "$w/$c" > "$w/$c.out" 2>&1 \
    || { grep -A3 '^Error' "$w/$c.out" | cut -c1-200 >&2; exit 1; }
  local n; n=$(grep -oE '[0-9]+ variants remaining' "$w/$c.out" | tail -1 | cut -d' ' -f1)
  echo "$c: chrX non-PAR common SNPs: ${n:-0}"
  [[ "${n:-0}" -ge 100 ]] || { echo "ERROR: $c: too few chrX SNPs; check chromosome codes" >&2; exit 1; }
}

# BF1 raw: the 1680-sample bfile that 01_preimpute_qc_bf1.sh reads (not the
# lab's QC'd 1546-sample raw/GRch37 pfile), so the lab's QC FAIL samples are checked too
R1=$BF/array/bf1/raw/BioFINDER1RAWDATA/GSP_2018-054-biofinder
{ printf '#FID\tIID\tSEX\n'; awk '{print $1 "\t" $2 "\tNA"}' "$R1.fam"; } > "$w/bf1.psam"
count_x bf1 --bed "$R1.bed" --bim "$R1.bim"

# BF2 raw: bfile, no sex in the fam
R2=$BF/array/bf2/raw/UB-2844_240229_PLINK_PCF_PLUS/1_BF-2_Genetic_Data
{ printf '#FID\tIID\tSEX\n'; awk '{print $1 "\t" $2 "\tNA"}' "$R2.fam"; } > "$w/bf2.psam"
count_x bf2 --bed "$R2.bed" --bim "$R2.bim"

"$PY" - "$w" "$OUT" "$BF" "$PET" <<'PY'
import sys, pandas as pd
w, out, bf, pet_path = sys.argv[1:5]
pet = pd.read_csv(pet_path, usecols=["sid", "Visit", "gender_baseline_variable"], dtype=str)
pet = pet[pd.to_numeric(pet["Visit"], errors="coerce") == 0]
pet["sex_reported"] = pet["gender_baseline_variable"].str.strip().map(
    {"0": "M", "0.0": "M", "1": "F", "1.0": "F"})
pet = pet.dropna(subset=["sex_reported"])[["sid", "sex_reported"]].drop_duplicates()
if pet["sid"].duplicated().any():
    sys.exit(f"{pet['sid'].duplicated().sum()} sids have two reported sexes at Visit 0")
zs = lambda s: s.astype(str).str.strip().str.replace(r"^0+(?=.)", "", regex=True)
for c in ["bf1", "bf2"]:
    sc = pd.read_csv(f"{w}/{c}.scount", sep="\t", dtype=str).rename(columns={"#IID": "IID", "#FID": "FID"})
    n = sum(sc[k].astype(int) for k in ["HOM_REF_CT", "HOM_ALT_SNP_CT", "HET_SNP_CT"])
    sc["x_het"] = sc["HET_SNP_CT"].astype(int) / n
    sc["sex_genetic"] = "ambiguous"
    sc.loc[sc["x_het"] < 0.10, "sex_genetic"] = "M"
    sc.loc[sc["x_het"] > 0.20, "sex_genetic"] = "F"
    key = pd.read_csv(f"{bf}/pheno/clean/{c}_sid_iid.csv", dtype=str)
    sc["iid_z"], key["iid_z"] = zs(sc["IID"]), zs(key["IID"])
    m = sc[["IID", "iid_z", "x_het", "sex_genetic"]].merge(key[["iid_z", "sid"]], on="iid_z", how="left")
    assert len(m) == len(sc), "IID repeated in the key"
    m = m.merge(pet, on="sid", how="left")
    assert len(m) == len(sc), "sid repeated in PET"
    has = m["sex_reported"].notna()
    m["status"] = "no_reported_sex"
    m.loc[has & (m["sex_genetic"] == m["sex_reported"]), "status"] = "match"
    m.loc[has & (m["sex_genetic"] == "ambiguous"), "status"] = "ambiguous"
    m.loc[has & m["sex_genetic"].isin(["M", "F"]) & (m["sex_genetic"] != m["sex_reported"]), "status"] = "mismatch"
    m.drop(columns="iid_z").to_csv(f"{out}/{c}_array_sex_check.tsv", sep="\t", index=False)
    bad = m[m["status"].isin(["mismatch", "ambiguous"])]
    bad[["IID"]].to_csv(f"{out}/{c}_sex_exclude.txt", index=False, header=False)
    print(f"\n== {c}: {len(m)} raw array samples")
    print("genetic sex:", m["sex_genetic"].value_counts().to_dict())
    print("status:", m["status"].value_counts().to_dict())
    print(f"excluded: {len(bad)} -> {out}/{c}_sex_exclude.txt")
    for _, r in bad.iterrows():
        print(f"  {r['status']:9s} genetic {r['sex_genetic']:9s} reported {r['sex_reported']}  x_het {r['x_het']:.3f}")
print("DONE")
PY
