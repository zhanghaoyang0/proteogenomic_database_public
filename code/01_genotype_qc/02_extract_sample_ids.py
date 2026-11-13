#!/usr/bin/env python3
"""Extract the identity columns of a BF phenotype export into a small table.

Writes ~/db/secure_data/bf/pheno/clean/<export_dir_name>_ids.csv: one row per
export row (sid x Visit), only the id, study, visit, age/sex/APOE and main date
columns (ID_COLS). ID and sample-identity checks (sid aliases, sample ids,
WGS/array keys) read this instead of the 1.5 GB export. Regenerable; rerun when
a new export arrives or the export is corrected.

The log carries counts only: rows, sids, filled cells per column, and for the
cross-person ids (aid, mid, subject_id) how many values are shared by more than
one sid -- the sign that a column links sids of one person.

Usage:
  python3 code/01_genotype_qc/02_extract_sample_ids.py > code/process_pheno/log/00_extract_ids.log 2>&1
  (--export-dir ~/db/secure_data/bf/pheno/raw/<other_export>)
"""
import argparse
import glob
import os
import sys
from os.path import expanduser

import pandas as pd

DEFAULT_EXPORT_DIR = expanduser("~/db/secure_data/bf/pheno/raw/haoyang__20260227_112046_PET")
OUT_DIR = expanduser("~/db/secure_data/bf/pheno/clean")
ID_COLS = [
    "uid", "Study", "MetaStudy", "sid", "Visit", "visit_date", "data_index",
    "aid", "mid", "subject_id", "OBARCODE", "blood_sample_id", "blood_sample_id_pc", "csf_sample_id",
    "age", "gender_baseline_variable", "apoe_genotype_baseline_variable", "excluded", "reason_missing",
    "blood_sample_date", "csf_date", "proc_mri_date", "proc_flutemetamol_pet_date",
]
# ids that may be shared by several sids of one person (checked, not assumed)
LINK_COLS = ["aid", "mid", "subject_id", "OBARCODE", "blood_sample_id", "csf_sample_id"]


def find_one(pattern, exclude=None):
    hits = [p for p in glob.glob(pattern) if exclude is None or exclude not in p]
    if len(hits) != 1:
        sys.exit(f"ERROR: expected exactly one match for {pattern}, got {len(hits)}")
    return hits[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", default=DEFAULT_EXPORT_DIR)
    parser.add_argument("--out-dir", default=OUT_DIR)
    args = parser.parse_args()

    export_dir = expanduser(args.export_dir.rstrip("/"))
    data_path = find_one(os.path.join(export_dir, "*_x.csv"))
    header = pd.read_csv(data_path, nrows=0).columns.tolist()
    missing = [c for c in ID_COLS if c not in header]
    if missing:
        sys.exit(f"ERROR: export lacks columns: {missing}")

    print(f"export: {data_path} ({len(header)} columns)", flush=True)
    df = pd.read_csv(data_path, usecols=ID_COLS, dtype=str)[ID_COLS]
    for c in df.columns:
        df[c] = df[c].str.strip()
        df[c] = df[c].mask(df[c].isin(["", "nan", "NA"]))
    # Rows without a sid exist in the export (614 in haoyang__20260227_112046,
    # 2026-10-13); kept and described, since an id column may still name them.
    no_sid = df["sid"].isna()
    if no_sid.any():
        print(f"rows without a sid: {int(no_sid.sum())}; their filled cells per column:",
              {c: int(df.loc[no_sid, c].notna().sum()) for c in ID_COLS if df.loc[no_sid, c].notna().any()})

    print(f"rows {len(df)}; sids {df['sid'].nunique()}; "
          f"sids per Study: {df.groupby('Study', dropna=False)['sid'].nunique().to_dict()}")
    print("filled cells per column:", {c: int(df[c].notna().sum()) for c in ID_COLS})
    print("values shared by more than one sid (a cross-person link if > 0):")
    for c in LINK_COLS:
        v = df[["sid", c]].dropna().drop_duplicates()
        per = v.groupby(c)["sid"].nunique()
        print(f"  {c:<16} distinct values {len(per):>6}; on >1 sid {int((per > 1).sum()):>5}; "
              f"max sids per value {int(per.max()) if len(per) else 0}")

    os.makedirs(args.out_dir, exist_ok=True)
    out_path = os.path.join(args.out_dir, f"{os.path.basename(export_dir)}_ids.csv")
    df.to_csv(out_path, index=False)
    print(f"written: {out_path} ({len(df)} rows x {df.shape[1]} columns, "
          f"{os.path.getsize(out_path) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
