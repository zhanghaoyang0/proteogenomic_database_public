#!/usr/bin/env python3
"""
Process covariates for PLINK: age, sex, PCs, mean_protein.
All proteins are baseline (Visit==0) only, so age/sex come from the PET export
(a superset of the main export, checked 2026-10-09) at Visit==0
and no per-protein measurement-visit column is emitted.
"""
import sys
import pandas as pd
from os.path import expanduser
import os

code_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(code_dir)
if code_dir not in sys.path:
    sys.path.insert(0, code_dir)

print("Starting 02_process_covar.py...")
print("Importing utility functions...")
from utils.process_data import (
    norm_id,
    format_covariates_for_plink,
)
print("Imports successful!")

# Same override as 01_process_protein.py, so both can be pointed at a sandbox
# copy for a dry run while a GWAS is reading the live covar file.
protein_out_dir = expanduser(os.environ.get(
    "PROTEIN_PLINK_DIR", "~/db/secure_data/bf/array/gwas_input/protein",
)) + "/"


def valid_id_map(df):
    df = df.copy()
    df["FID"] = norm_id(df["FID"])
    df["IID"] = norm_id(df["IID"])
    valid = (
        df["FID"].notna()
        & df["IID"].notna()
        & df["FID"].ne("Removed in QC")
        & df["IID"].ne("Removed in QC")
    )
    out = df.loc[valid].drop_duplicates(subset=["sid", "FID", "IID"]).copy()
    # plink2 psam files have no FID column → implicit FID=0. Match it.
    out["FID"] = "0"
    return out


def load_pca(path):
    pca = pd.read_csv(path, sep=r"\s+", header=None, comment="#")
    if pca.shape[1] == 12:
        pca.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, 11)]
    elif pca.shape[1] == 11:
        pca.columns = ['IID'] + [f'PC{i}' for i in range(1, 11)]
        pca.insert(0, 'FID', pca['IID'])
    else:
        raise ValueError(f"Unexpected PCA column count in {path}: {pca.shape[1]}")
    pca['IID'] = norm_id(pca['IID'])
    return pca.rename(columns={f'PC{i}': f'pc{i}' for i in range(1, 11)})

# ==============================================================
# Load covariates from df_main (all visits)
# ==============================================================
print("\n" + "="*60)
print("Loading covariates...")
print("="*60)

df_main = pd.read_csv(expanduser("~/db/secure_data/bf/pheno/raw/haoyang__20260227_112046_PET/haoyang__20260227_112046_x.csv"),
                      low_memory=False, usecols=['sid', 'Visit', 'age', 'gender_baseline_variable'])
print(f"Loaded df_main: {df_main.shape[0]} rows")

# Extract relevant columns without filtering Visit
df_covar_main = df_main[['sid', 'Visit', 'age', 'gender_baseline_variable']].copy()
df_covar_main['sid'] = norm_id(df_covar_main['sid'])
df_covar_main = df_covar_main.rename(columns={'gender_baseline_variable': 'sex'})

# Remove duplicates by (sid, Visit)
print(f"Before dedup: {df_covar_main.shape[0]} rows")
df_covar_main = df_covar_main.drop_duplicates(subset=['sid', 'Visit'], keep='first')
print(f"After dedup: {df_covar_main.shape[0]} rows")

# ==============================================================
# Load PCA results
# ==============================================================
print("\nLoading PCA results...")
bf1_pca = load_pca(expanduser("~/db/secure_data/bf/array/bf1/imputed_qc/bf1_hg19_pca.eigenvec"))
bf2_pca = load_pca(expanduser("~/db/secure_data/bf/array/bf2/imputed_qc/bf2_hg19_pca.eigenvec"))

print(f"bf1_pca: {bf1_pca.shape[0]} samples")
print(f"bf2_pca: {bf2_pca.shape[0]} samples")

# ==============================================================
# Load mean_protein from protein script
# ==============================================================
print("\nLoading mean_protein files from protein script...")
bf1_mean_protein = pd.read_csv(protein_out_dir + 'bf1_mean_protein.csv')
bf2_mean_protein = pd.read_csv(protein_out_dir + 'bf2_mean_protein.csv')

bf1_mean_protein['sid'] = norm_id(bf1_mean_protein['sid'])
bf2_mean_protein['sid'] = norm_id(bf2_mean_protein['sid'])

print(f"bf1_mean_protein: {bf1_mean_protein.shape[0]} samples, {bf1_mean_protein.shape[1]-1} columns")
print(f"bf2_mean_protein: {bf2_mean_protein.shape[0]} samples, {bf2_mean_protein.shape[1]-1} columns")

# ==============================================================
# Load sid_iid mappings and find an approximate subject-level protein visit
# ==============================================================
print("\nLoading sid_iid mappings...")

for cohort_label, pca_df, mean_protein_df in [
    ("bf1", bf1_pca, bf1_mean_protein),
    ("bf2", bf2_pca, bf2_mean_protein),
]:
    print(f"\n{cohort_label}:")

    # Load sid_iid mapping
    sid_iid = pd.read_csv(expanduser(f"~/db/secure_data/bf/pheno/clean/{cohort_label}_sid_iid.csv"), dtype=str)
    sid_iid['sid'] = norm_id(sid_iid['sid'])
    sid_iid = valid_id_map(sid_iid)
    print(f"  sid_iid mapping: {sid_iid.shape[0]} samples")

    baseline = df_covar_main[df_covar_main['Visit'].astype(float) == 0][['sid', 'age', 'sex']].copy()

    print(f"  Baseline rows (Visit==0): {baseline.shape[0]} sids")

    # Merge: mean_protein + baseline age/sex + FID/IID + PCA
    covar_df = mean_protein_df.merge(baseline, on='sid', how='left')
    covar_df = covar_df.merge(sid_iid[['sid', 'FID', 'IID']], on='sid', how='inner')

    print(f"  After merging with sid_iid: {covar_df.shape[0]} samples")

    # Create covariate dataframe: FID IID age sex pc1-pc10 mean_protein_*
    mp_cols = [c for c in covar_df.columns if c.startswith('mean_protein')]
    covar_formatted = format_covariates_for_plink(
        covar_df[['FID', 'IID'] + mp_cols],
        df_pca=pca_df,
        df_covar=covar_df[['IID', 'age', 'sex']]
    )

    out_path = protein_out_dir + f"{cohort_label}_covar.phe"
    covar_formatted.to_csv(out_path, sep="\t", index=False, na_rep="NA")

    cols_display = [c for c in covar_formatted.columns if not c.startswith('mean_protein')]
    print(f"  Saved {cohort_label}_covar.phe: {covar_formatted.shape[0]} samples, {covar_formatted.shape[1]-2} covariates")
    print(f"    Columns: {' '.join(cols_display)}")

print("\n" + "="*60)
print("Processing complete!")
print("="*60)
