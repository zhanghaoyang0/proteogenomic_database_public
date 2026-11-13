#!/usr/bin/env python3
"""
Process BioFinder protein data: QC, output protein matrix CSV + GWAS phenotype files.
Uses baseline (Visit==0) protein values only. Subjects without a Visit==0 measurement
for a given protein get NA for that protein.
No INT transformation. Platform QC matrices keep post-QC proteins; GWAS .phe files
only include proteins with >1000 non-missing values among GWAS-eligible samples.

Outputs:
  ~/db/secure_data/bf/protein/qc/{cohort}_{tag}_protein.csv  (蛋白矩阵：sid + proteins)
  ~/db/secure_data/bf/array/gwas_input/protein/{cohort}_{tag}.phe (GWAS用)
"""
import sys
import numpy as np
import pandas as pd
from os.path import expanduser
import os
import re

code_dir = os.path.abspath(os.path.dirname(__file__))
os.chdir(code_dir)
if code_dir not in sys.path:
    sys.path.insert(0, code_dir)

print("Starting 01_process_protein.py...")
print("Importing utility functions...")
from utils.process_data import (
    norm_id,
    qc_olink,
    qc_nulisa,
    qc_somalogic,
    add_nulisa2_lod,
    format_proteins_for_plink,
)
from utils.qc_report import qc_summary_row
print("Imports successful!")

# ==============================================================
# Load bf sample IDs
# ==============================================================
print("\n" + "="*60)
print("Loading bf sample IDs...")
print("="*60)
# BF1 sids come from the genetic key, not the array sheet: the sheet's Clinical_ID
# spells NOMAS2 in upper case and Parkinson without the dash, so those never matched
# the protein sids (code/process_pheno/docs/ids.md, "sid -> genotype IID").
file_bf1 = expanduser('~/db/secure_data/bf/pheno/clean/bf1_sid_iid.csv')
file_bf2 = expanduser('~/db/secure_data/bf/array/bf2/raw/GWAS_Full_Data_BF-2.xlsx')

print(f"Reading bf1 file: {file_bf1}")
bf1 = pd.read_csv(file_bf1, usecols=['sid'], dtype=str)
print(f"Reading bf2 file: {file_bf2}")
bf2 = pd.read_excel(file_bf2, usecols=['sid'],dtype=str)

bf1['sid'] = norm_id(bf1['sid'])
bf2['sid'] = norm_id(bf2['sid'])
bf1_id = set(bf1['sid'].dropna())
bf2_id = set(bf2['sid'])
if bf1_id & bf2_id:
    raise ValueError(f"BF1/BF2 sid overlap: {sorted(bf1_id & bf2_id)[:5]}")

print(f"Total unique bf1_id: {len(bf1_id)}")
print(f"Total unique bf2_id: {len(bf2_id)}")

# ==============================================================
# Load sid→FID/IID mappings
# ==============================================================
print("\n" + "="*60)
print("Loading sid→FID/IID mappings...")
print("="*60)

bf1_id_map = pd.read_csv(expanduser("~/db/secure_data/bf/pheno/clean/bf1_sid_iid.csv"), dtype=str)
bf2_id_map = pd.read_csv(expanduser("~/db/secure_data/bf/pheno/clean/bf2_sid_iid.csv"), dtype=str)

bf1_id_map['sid'] = norm_id(bf1_id_map['sid'])
bf1_id_map['IID'] = norm_id(bf1_id_map['IID'])
bf2_id_map['sid'] = norm_id(bf2_id_map['sid'])
bf2_id_map['IID'] = norm_id(bf2_id_map['IID'])

def _valid_id_map(df):
    df = df.copy()
    df['FID'] = norm_id(df['FID'])
    df['IID'] = norm_id(df['IID'])
    valid = (
        df['FID'].notna()
        & df['IID'].notna()
        & df['FID'].ne('Removed in QC')
        & df['IID'].ne('Removed in QC')
    )
    out = df.loc[valid].drop_duplicates(subset=['sid', 'FID', 'IID']).copy()
    # plink2 psam files in this project have no FID column → implicit FID=0.
    # Override to match so --pheno/--covar can find loaded sample IDs.
    out['FID'] = '0'
    return out

bf1_id_map = _valid_id_map(bf1_id_map)
bf2_id_map = _valid_id_map(bf2_id_map)

print(f"bf1_id_map: {bf1_id_map.shape[0]} samples")
print(f"bf2_id_map: {bf2_id_map.shape[0]} samples")

# ==============================================================
# Load protein feather files, baseline (Visit==0) only
# ==============================================================
print("\nLoading protein feather files...")

olink_dir = expanduser("~/db/secure_data/bf/pheno/raw/")
nulisa_dir = olink_dir + "haoyang__20260430_134828_nulisa/"

def _load_feather_long(path):
    """Same prep as _load_feather but keep all visits (no groupby)."""
    df = pd.read_feather(path)
    df['sid'] = norm_id(df['sid'])
    df = df[df['sid'].isin(bf1_id | bf2_id)].copy()
    df = df.dropna(axis=1, how='all')
    df = df.drop(columns=[c for c in df.columns if c.startswith('~RAW')])
    if 'Visit' in df.columns:
        df = df.sort_values(['sid', 'Visit'])
    return df.reset_index(drop=True)


def _load_olink_ht(path, kind):
    """Load Olink Explore HT (CSF / Plasma feather) and normalize to the layout
    expected by the existing Olink pipeline.

    Both fluids now come from the official database exports (olink6 = CSF,
    olink7 = plasma), which are feather and ship sid/Visit/Study/aid/... as real
    columns alongside uid. The earlier shorena_20260521 hand-delivery had the
    same values but a different envelope (plasma was parquet with only uid).

    - kind='csf'   : keep only rows labelled Explore HT (exclude Explore 3072
                     bridging-reference rows), then drop Bridged_o2csf /
                     BridgingSamples_o2csf and use __NPX / __LOD as-is.
    - kind='plasma': parse uid (Study__sid__Visit) → sid/Visit; replace
                     __NPX/__LOD with __PCNormalizedNPX/__PCNormalizedLOD
                     (plate/bridge-corrected; same proteins & samples as raw NPX).
                     sid is derived from uid rather than taken from the shipped
                     sid column so the sample set is provably unchanged against
                     the previous delivery.

    HT QC columns __AssayQC / __SampleQC are renamed to __Assay_Warning / __QC_Warning
    so qc_olink can be reused unchanged.
    """
    if kind not in ('csf', 'plasma'):
        raise ValueError(f"unknown kind: {kind}")
    df = pd.read_feather(path)

    if kind == 'plasma':
        if 'uid' not in df.columns:
            raise ValueError(f"plasma delivery has no uid column: {path}")
        parts = df['uid'].astype(str).str.split('__', n=2, expand=True)
        df['sid'] = norm_id(parts[1])
        df['Visit'] = pd.to_numeric(parts[2], errors='coerce')
    if 'sid' not in df.columns:
        raise ValueError(f"delivery has no sid column: {path}")

    # Delivery-level metadata shipped by the database export. Not markers, and
    # they collide across the CSF/plasma merges downstream, so drop them here.
    df = df.drop(columns=['uid', 'Project', 'Study', 'MetaStudy', 'aid',
                          'visit_date', 'data_index'], errors='ignore')

    df['sid'] = norm_id(df['sid'])
    df = df[df['sid'].isin(bf1_id | bf2_id)].copy()

    # plasma: swap in PC-normalized NPX/LOD, drop raw NPX/LOD + Adj_factor.
    if kind == 'plasma':
        df = df.drop(columns=[c for c in df.columns if c.endswith('__NPX') and 'PCNormalized' not in c], errors='ignore')
        df = df.drop(columns=[c for c in df.columns if c.endswith('__LOD') and 'PCNormalized' not in c], errors='ignore')
        df = df.drop(columns=[c for c in df.columns if c.endswith('__Adj_factor')], errors='ignore')
        rename = {}
        for c in df.columns:
            if c.endswith('__PCNormalizedNPX'):
                rename[c] = c.replace('__PCNormalizedNPX', '__NPX')
            elif c.endswith('__PCNormalizedLOD'):
                rename[c] = c.replace('__PCNormalizedLOD', '__LOD')
        df = df.rename(columns=rename)

    # CSF delivery also contains Explore 3072 reference measurements used for
    # bridging. They are processed separately from the original 3072 panels and
    # must not be duplicated in the ExploreHT panel.
    if kind == 'csf':
        platform_col = "ExploreHT_o6csf__Bridged_o2csf"
        if platform_col not in df.columns:
            raise ValueError(f"missing expected CSF platform column: {platform_col}")
        df = df.loc[df[platform_col].eq("Explore HT")].copy()
        df = df.drop(columns=[c for c in df.columns if 'Bridg' in c and not c.endswith('__BridgingRecommendation')], errors='ignore')

    # rename HT QC columns so qc_olink finds them under the legacy names.
    qc_rename = {}
    for c in df.columns:
        if c.endswith('__AssayQC'):
            qc_rename[c] = c.replace('__AssayQC', '__Assay_Warning')
        elif c.endswith('__SampleQC'):
            qc_rename[c] = c.replace('__SampleQC', '__QC_Warning')
    df = df.rename(columns=qc_rename)
    # drop BridgingRecommendation (unused by qc_olink)
    df = df.drop(columns=[c for c in df.columns if c.endswith('__BridgingRecommendation')], errors='ignore')

    df = df.dropna(axis=1, how='all')
    if 'Visit' in df.columns:
        df = df.sort_values(['sid', 'Visit'])
    return df.reset_index(drop=True)


def _collapse_baseline(df):
    """Keep one row per sid at Visit==0. Subjects without a V0 row are dropped."""
    if df.empty:
        return df
    out = df
    if 'Visit' in out.columns:
        v = pd.to_numeric(out['Visit'], errors='coerce')
        out = out.loc[v == 0]
    return out.drop_duplicates(subset=['sid'], keep='first').reset_index(drop=True)


def sanitize_phe_columns(columns):
    """Make protein column names PLINK/file-safe without changing biology."""
    protected = {"FID", "IID", "sid", "mean_protein"}
    sanitized = []
    for col in columns:
        if col in protected or str(col).startswith("mean_protein_"):
            sanitized.append(col)
            continue
        clean = (
            str(col)
            .replace(" ", "_")
            .replace("-", "_")
            .replace(".", "_")
            .replace("Î²", "b")
            .replace("β", "b")
        )
        sanitized.append(clean)
    return sanitized


olink_long = _load_feather_long(olink_dir + "haoyang__20260430_111026_olink/haoyang__20260430_111026.feather")
protein_olink = _collapse_baseline(olink_long)

# Olink Explore HT, official database exports (olink6 = CSF, olink7 = plasma).
# In the CSF delivery native Explore 3072 rows are bridged onto the HT scale; the
# pQTL panel keeps only rows labelled "Explore HT", while the network pipeline
# also uses the bridged 3072 rows. The CSF export ships data only for the
# bridging-assessed assays, so the assays it leaves empty are dropped by
# dropna(how='all') below -- see README.md for what that costs the panel.
ht_csf_long = _load_olink_ht(
    olink_dir + "haoyang__20260829_111908_olink6/haoyang__20260829_111908.feather",
    kind='csf',
)
ht_plasma_long = _load_olink_ht(
    olink_dir + "haoyang__20260829_120718_olink7/haoyang__20260829_120718.feather",
    kind='plasma',
)
ht_csf_b = _collapse_baseline(ht_csf_long)
ht_plasma_b = _collapse_baseline(ht_plasma_long)
print(
    f"  HT long rows — csf: {ht_csf_long.shape[0]} | plasma: {ht_plasma_long.shape[0]}"
    f" | baseline subjects — csf: {ht_csf_b.shape[0]} | plasma: {ht_plasma_b.shape[0]}"
)
# Merge HT panels onto the baseline Olink table (per-sid wide join). Drop HT
# helper non-marker columns before merge to avoid Visit/sid collisions.
for _ht in (ht_csf_b, ht_plasma_b):
    _ht.drop(columns=[c for c in _ht.columns if c == 'Visit'], errors='ignore', inplace=True)
protein_olink = protein_olink.merge(ht_csf_b, on='sid', how='outer')
protein_olink = protein_olink.merge(ht_plasma_b, on='sid', how='outer')
print(f"  After HT merge — protein_olink subjects: {protein_olink.shape[0]}, columns: {protein_olink.shape[1]}")

soma_long = _load_feather_long(olink_dir + "haoyang__20260429_234900_somalogic/haoyang__20260429_234900.feather")

nulisa_long = _load_feather_long(nulisa_dir + "haoyang__20260430_134828.feather")

print(
    f"  Long rows — olink: {olink_long.shape[0]} | soma: {soma_long.shape[0]} | nulisa: {nulisa_long.shape[0]}"
)
print(f"  After baseline (Visit==0) collapse — olink subjects: {protein_olink.shape[0]}")

# ==============================================================
# Olink: assign LOD_Flag, then QC per panel
# ==============================================================
# NOTE on LOD_Flag values "WARN"/"PASS":
#   These match the strings Olink ships in feather-provided LOD_Flag columns
#   and are used ONLY to compute protein-level detection rate (Step 1, ≥70%)
#   and to pick high-confidence proteins for mean_protein (Step 4, ≥90%).
#   They are NOT passed through bad_flags below — below-LOD individual
#   datapoints are KEPT per Olink/Alexa et al. recommendation. The
#   bad_flags={"WARN", ...} filter applies only to __Assay_Warning /
#   __QC_Warning columns inside qc_olink.
new_cols = {}
for col in protein_olink.columns:
    if col.endswith("__NPX"):
        prefix = col.rsplit("__", 1)[0]
        lod_col = prefix + "__LOD"
        lod_flag_col = prefix + "__LOD_Flag"
        if lod_col in protein_olink.columns and lod_flag_col not in protein_olink.columns:
            both_valid = protein_olink[col].notna() & protein_olink[lod_col].notna()
            flags = np.where(protein_olink[col] < protein_olink[lod_col], "WARN", "PASS").astype(object)
            flags[~both_valid] = None
            new_cols[lod_flag_col] = flags
n_assigned = len(new_cols)
protein_olink = pd.concat([protein_olink, pd.DataFrame(new_cols)], axis=1)
n_npx = sum(c.endswith('__NPX') for c in protein_olink.columns)
print(f"  Olink: {n_npx} proteins | LOD_Flag assigned for {n_assigned} proteins")

qc_dir = expanduser(os.environ.get(
    "PROTEIN_QC_DIR",
    "~/db/secure_data/bf/protein/qc",
)) + "/"
plink_out_dir = expanduser(os.environ.get(
    "PROTEIN_PLINK_DIR",
    "~/db/secure_data/bf/array/gwas_input/protein",
)) + "/"
info_dir = expanduser(os.environ.get(
    "PROTEIN_INFO_DIR",
    os.path.abspath(os.path.join(os.path.dirname(__file__), "../../info/bf_protein")),
))
os.makedirs(qc_dir, exist_ok=True)
os.makedirs(plink_out_dir, exist_ok=True)
os.makedirs(info_dir, exist_ok=True)

mean_proteins = {'bf1': {}, 'bf2': {}}
qc_report_rows = []
GWAS_MIN_PROTEIN_N = 1000

def _gwas_cols_above_min(df, protein_cols):
    """Keep GWAS protein columns with non-missing N > GWAS_MIN_PROTEIN_N."""
    keep = []
    for col in protein_cols:
        n_nonmissing = int(df[col].notna().sum())
        if n_nonmissing > GWAS_MIN_PROTEIN_N:
            keep.append(col)
    return keep

def _olink_panel(col):
    """Return raw panel string from PROTEIN__Neurology_o2csf__NPX."""
    parts = col.split('__')
    if len(parts) >= 3 and re.search(r'_o\d+(csf|plasma)$', parts[-2], re.IGNORECASE):
        return parts[-2]
    return None

combos = sorted(set(
    _olink_panel(c) for c in protein_olink.columns if c.endswith('__NPX')
) - {None})
print(f"  Olink panels: {combos}")

for tag in combos:
    st_cols = ['sid'] + [c for c in protein_olink.columns if _olink_panel(c) == tag]
    if 'Visit' in protein_olink.columns and 'Visit' not in st_cols:
        st_cols.append('Visit')
    if sum(c.endswith('__NPX') for c in st_cols) == 0:
        continue
    protein_olink_st = protein_olink[st_cols]
    mp_key = f'mean_protein_{tag}'

    for cohort_label, cohort_id, id_map in [
        ('bf1', bf1_id, bf1_id_map),
        ('bf2', bf2_id, bf2_id_map),
    ]:
        npx_cols_st = [c for c in protein_olink_st.columns if c.endswith('__NPX')]
        cohort_rows = protein_olink_st[protein_olink_st['sid'].isin(cohort_id)]
        if cohort_rows[npx_cols_st].notna().sum().sum() == 0:
            continue

        qc_df, steps = qc_olink(
            protein_olink_st,
            cohort_id,
            bad_flags={"WARN", "FAIL", "EXCLUDED", "MANUAL_WARN"},
        )
        n_qc = qc_df.shape[0]
        npx_count = sum(c.endswith('__NPX') for c in qc_df.columns)
        if n_qc == 0 or npx_count == 0:
            continue

        # Store mean_protein for covariate file
        mp_df = qc_df[['sid', 'mean_protein']].copy()
        mp_df = mp_df.rename(columns={'mean_protein': mp_key})
        mean_proteins[cohort_label][mp_key] = mp_df

        # File 1: protein matrix (all samples, sid-keyed)
        # Drop proteins with no non-missing values (e.g. assay-level QC FAIL/WARN
        # like APOE on ExploreHT plasma) so KNN and plink2 don't see empty cols.
        prot_cols_orig = [c for c in qc_df.columns if c.endswith('__NPX') and qc_df[c].notna().any()]
        prot_df_analysis = qc_df[['sid'] + prot_cols_orig].copy()
        # Drop the metric suffix from output column names: {assay}__{tag}__NPX -> {assay}__{tag}.
        prot_df_analysis = prot_df_analysis.rename(columns={c: c.removesuffix('__NPX') for c in prot_cols_orig})
        prot_cols_safe = sanitize_phe_columns(prot_df_analysis.columns)
        prot_df_analysis.columns = prot_cols_safe
        out_path_analysis = qc_dir + f'{cohort_label}_{tag}_protein.csv'
        prot_df_analysis.to_csv(out_path_analysis, index=False, na_rep='NA')

        # File 2: GWAS (genotyped samples only, FID/IID-keyed)
        # Subset to "kept" samples (≥1 non-missing in panel) so gwas ⊆ kept,
        # then dedupe in case id_map has multiple FID/IID per sid.
        panel_npx = [c for c in qc_df.columns if c.endswith("__NPX") and qc_df[c].notna().any()]
        qc_kept = qc_df.loc[qc_df[panel_npx].notna().any(axis=1)] if panel_npx else qc_df
        df_merged = qc_kept.merge(id_map, on='sid', how='inner').drop_duplicates(subset=['FID', 'IID'])
        gwas_npx = _gwas_cols_above_min(df_merged, panel_npx)
        if len(gwas_npx) < len(panel_npx):
            print(f"  {cohort_label}_{tag}: excluding {len(panel_npx) - len(gwas_npx)} proteins with GWAS non-missing N <= {GWAS_MIN_PROTEIN_N}")
        # Skip the GWAS .phe entirely when no protein passes the N filter (avoids
        # a confusing FID/IID-only file, e.g. small bf1 panels).
        if not gwas_npx:
            print(f"  {cohort_label}_{tag}: no protein passes GWAS N > {GWAS_MIN_PROTEIN_N}; skipping .phe")
        else:
            prot_df_gwas = format_proteins_for_plink(
                df_merged[['FID', 'IID', 'mean_protein'] + gwas_npx],
                gwas_npx,
                inverse_normal=False,
            )
            prot_df_gwas.columns = sanitize_phe_columns(prot_df_gwas.columns)
            prot_df_gwas = prot_df_gwas.drop(columns=['mean_protein'], errors='ignore')
            out_path_gwas = plink_out_dir + f'{cohort_label}_{tag}.phe'
            prot_df_gwas.to_csv(out_path_gwas, sep='\t', index=False, na_rep='NA')

        print(f"Saved {cohort_label}_{tag}: {prot_df_analysis.shape[0]} protein samples, {df_merged.shape[0]} GWAS samples")

        qc_report_rows.append(
            qc_summary_row(
                cohort_label, "Olink", tag, steps, qc_df, df_merged.shape[0],
                n_gwas_proteins=len(gwas_npx),
            )
        )

# ==============================================================
# NuLISA processing
# ==============================================================
lod_combined = pd.concat([
    pd.read_excel(nulisa_dir + 'BF1_LOD_nulisa2_CNS.xlsx'),
    pd.read_excel(nulisa_dir + 'BF2_LOD_nulisa2_CNS.xlsx'),
]).drop_duplicates(subset=['PlateID', 'TargetName'])
protein_nulisa = add_nulisa2_lod(nulisa_long.copy(), lod_combined)
protein_nulisa = _collapse_baseline(protein_nulisa)

npq_panels = sorted(set(
    c.split('__', 1)[1].removesuffix('_NPQ')
    for c in protein_nulisa.columns if c.endswith('_NPQ') and '__' in c
))

from collections import defaultdict
tag_to_panels = defaultdict(list)
for panel in npq_panels:
    tag_to_panels[panel].append(panel)

for tag, panels in sorted(tag_to_panels.items()):
    st_cols = ['sid'] + [c for c in protein_nulisa.columns
                         if any(f'_{p}_' in c for p in panels)]
    if 'Visit' in protein_nulisa.columns and 'Visit' not in st_cols:
        st_cols.append('Visit')
    if sum(c.endswith('_NPQ') and '__' in c for c in st_cols) == 0:
        continue
    protein_nulisa_st = protein_nulisa[st_cols]
    mp_key = f'mean_protein_{tag}'

    for cohort_label, cohort_id, id_map in [
        ('bf1', bf1_id, bf1_id_map),
        ('bf2', bf2_id, bf2_id_map),
    ]:
        npq_cols_st = [c for c in protein_nulisa_st.columns if c.endswith('_NPQ') and '__' in c]
        cohort_rows = protein_nulisa_st[protein_nulisa_st['sid'].isin(cohort_id)]
        if cohort_rows[npq_cols_st].notna().sum().sum() == 0:
            continue

        qc_df, steps = qc_nulisa(protein_nulisa_st, cohort_id)
        n_qc = qc_df.shape[0]
        # Drop proteins with no non-missing values.
        npq_cols = [c for c in qc_df.columns if c.endswith('_NPQ') and '__' in c and qc_df[c].notna().any()]
        if n_qc == 0 or len(npq_cols) == 0:
            continue

        # gwas ⊆ kept: filter to ≥1 non-missing in panel, then dedupe FID/IID.
        qc_kept = qc_df.loc[qc_df[npq_cols].notna().any(axis=1)] if npq_cols else qc_df
        df_merged = qc_kept.merge(id_map, on='sid', how='inner').drop_duplicates(subset=['FID', 'IID'])
        gwas_npq = _gwas_cols_above_min(df_merged, npq_cols)
        if len(gwas_npq) < len(npq_cols):
            print(f"  {cohort_label}_{tag}: excluding {len(npq_cols) - len(gwas_npq)} proteins with GWAS non-missing N <= {GWAS_MIN_PROTEIN_N}")

        mp_df = qc_df[['sid', 'mean_protein']].copy()
        mp_df = mp_df.rename(columns={'mean_protein': mp_key})
        mean_proteins[cohort_label][mp_key] = mp_df

        # File 1: protein matrix
        prot_cols_orig = [c for c in qc_df.columns if c.endswith('_NPQ') and '__' in c and qc_df[c].notna().any()]
        prot_df_analysis = qc_df[['sid'] + prot_cols_orig].copy()
        # Drop the metric suffix from output column names: {assay}__{tag}_NPQ -> {assay}__{tag}.
        rename_bare = {c: c.removesuffix('_NPQ') for c in prot_cols_orig}
        prot_df_analysis = prot_df_analysis.rename(columns=rename_bare)
        prot_cols_safe = sanitize_phe_columns(prot_df_analysis.columns)
        prot_df_analysis.columns = prot_cols_safe
        out_path_analysis = qc_dir + f'{cohort_label}_{tag}_protein.csv'
        prot_df_analysis.to_csv(out_path_analysis, index=False, na_rep='NA')

        # File 2: GWAS — skip entirely when no protein passes the N filter.
        if not gwas_npq:
            print(f"  {cohort_label}_{tag}: no protein passes GWAS N > {GWAS_MIN_PROTEIN_N}; skipping .phe")
        else:
            npx_rename = {c: c.removesuffix('_NPQ') + '__NPX' for c in gwas_npq}
            df_out = df_merged[['FID', 'IID', 'mean_protein'] + gwas_npq].rename(columns=npx_rename)
            df_out = format_proteins_for_plink(df_out, list(npx_rename.values()), inverse_normal=False)
            df_out.columns = sanitize_phe_columns(df_out.columns)
            df_out = df_out.drop(columns=['mean_protein'], errors='ignore')
            out_path_gwas = plink_out_dir + f'{cohort_label}_{tag}.phe'
            df_out.to_csv(out_path_gwas, sep='\t', index=False, na_rep='NA')

        print(f"Saved {cohort_label}_{tag}: {prot_df_analysis.shape[0]} protein samples, {df_merged.shape[0]} GWAS samples")

        qc_report_rows.append(
            qc_summary_row(
                cohort_label, "NuLISA", tag, steps, qc_df, df_merged.shape[0],
                n_gwas_proteins=len(gwas_npq),
            )
        )

# ==============================================================
# SomaLogic processing
# ==============================================================
soma_meta_path = olink_dir + 'haoyang__20260429_234900_somalogic/plasma_somalogic4.1_metadata_Jan2024.csv'
soma_meta = pd.read_csv(soma_meta_path)

soma_cols = ['sid'] + [c for c in soma_long.columns if c.endswith('__ANML')]
if 'Visit' in soma_long.columns and 'Visit' not in soma_cols:
    soma_cols.append('Visit')
if 'RowCheck' in soma_long.columns and 'RowCheck' not in soma_cols:
    soma_cols.append('RowCheck')
protein_soma_long = soma_long[soma_cols].copy()
protein_soma = _collapse_baseline(protein_soma_long)
if 'RowCheck' in protein_soma.columns:
    rowcheck = protein_soma['RowCheck'].astype('string').str.upper()
    n_rowcheck_flag = int(rowcheck.eq('FLAG').sum())
    if n_rowcheck_flag:
        print(f"  SomaLogic: dropping {n_rowcheck_flag} baseline rows with RowCheck=FLAG")
    protein_soma = protein_soma.loc[~rowcheck.eq('FLAG')].drop(columns=['RowCheck']).reset_index(drop=True)

print("  SomaLogic: no LOD information")
for cohort_label, cohort_id, id_map in [
    ('bf1', bf1_id, bf1_id_map),
    ('bf2', bf2_id, bf2_id_map),
]:
    qc_df, steps = qc_somalogic(protein_soma, cohort_id, soma_meta)
    n_qc = qc_df.shape[0]
    # Drop proteins with no non-missing values.
    anml_cols = [c for c in qc_df.columns if c.endswith('__ANML') and qc_df[c].notna().any()]
    if n_qc == 0 or len(anml_cols) == 0:
        continue

    # gwas ⊆ kept: filter to ≥1 non-missing in panel, then dedupe FID/IID.
    qc_kept = qc_df.loc[qc_df[anml_cols].notna().any(axis=1)] if anml_cols else qc_df
    df_merged = qc_kept.merge(id_map, on='sid', how='inner').drop_duplicates(subset=['FID', 'IID'])
    gwas_anml = _gwas_cols_above_min(df_merged, anml_cols)
    if len(gwas_anml) < len(anml_cols):
        print(f"  {cohort_label}_somalogic_plasma: excluding {len(anml_cols) - len(gwas_anml)} proteins with GWAS non-missing N <= {GWAS_MIN_PROTEIN_N}")

    mp_key = 'mean_protein_somalogic_plasma'
    mp_df = qc_df[['sid', 'mean_protein']].copy()
    mp_df = mp_df.rename(columns={'mean_protein': mp_key})
    mean_proteins[cohort_label][mp_key] = mp_df

    # File 1: protein matrix
    prot_df_analysis = qc_df[['sid'] + anml_cols].copy()
    # Drop the metric suffix from output column names: {anml-col} -> bare seq id (e.g. seq_9984_12).
    prot_df_analysis = prot_df_analysis.rename(columns={c: c.split('__')[1] for c in anml_cols})
    prot_cols_safe = sanitize_phe_columns(prot_df_analysis.columns)
    prot_df_analysis.columns = prot_cols_safe
    out_path_analysis = qc_dir + f'{cohort_label}_somalogic_plasma_protein.csv'
    prot_df_analysis.to_csv(out_path_analysis, index=False, na_rep='NA')

    # File 2: GWAS — skip entirely when no protein passes the N filter.
    if not gwas_anml:
        print(f"  {cohort_label}_somalogic_plasma: no protein passes GWAS N > {GWAS_MIN_PROTEIN_N}; skipping .phe")
    else:
        npx_rename = {c: c.split('__')[1] + '__NPX' for c in gwas_anml}
        df_out = df_merged[['FID', 'IID', 'mean_protein'] + gwas_anml].rename(columns=npx_rename)
        df_out = format_proteins_for_plink(df_out, list(npx_rename.values()), inverse_normal=False)
        df_out.columns = sanitize_phe_columns(df_out.columns)
        df_out = df_out.drop(columns=['mean_protein'], errors='ignore')
        out_path_gwas = plink_out_dir + f'{cohort_label}_somalogic_plasma.phe'
        df_out.to_csv(out_path_gwas, sep='\t', index=False, na_rep='NA')

    print(f"Saved {cohort_label}_somalogic_plasma: {prot_df_analysis.shape[0]} protein samples, {df_merged.shape[0]} GWAS samples")

    qc_report_rows.append(
        qc_summary_row(
            cohort_label, "SomaLogic", "somalogic_plasma", steps, qc_df, df_merged.shape[0],
            n_gwas_proteins=len(gwas_anml),
        )
    )

# ==============================================================
# QC summary table (all cohorts / panels)
# ==============================================================
if qc_report_rows:
    _cols = [
        "cohort", "platform", "panel",
        "start(sample x protein)", "remove_70LOD (sample x protein)",
        "kept (sample x protein)",
        "gwas (sample x protein)",
        "n_high_confident_protein",
        "n_aptamers_meta_pass",
    ]
    rep_df = pd.DataFrame(qc_report_rows)
    rep_df = rep_df[_cols]
    rep_path = os.path.join(info_dir, "qc_panel_summary.csv")
    rep_df.to_csv(rep_path, index=False)
    print(f"\nSaved QC panel summary: {rep_path} ({rep_df.shape[0]} rows)")

# ==============================================================
# Save mean_protein and protein_age files
# ==============================================================
print("\n" + "="*60)
print("Saving intermediate files...")
print("="*60)

for cohort_label in ['bf1', 'bf2']:
    # Save mean_protein
    mp_dict = mean_proteins[cohort_label]
    if mp_dict:
        mp_list = list(mp_dict.values())
        result = mp_list[0].copy()
        for df in mp_list[1:]:
            result = result.merge(df, on='sid', how='outer')
        out_path = plink_out_dir + f'{cohort_label}_mean_protein.csv'
        result.to_csv(out_path, index=False)
        print(f"  Saved {cohort_label}_mean_protein.csv: {result.shape[0]} samples, {result.shape[1]-1} columns")

print("\n" + "="*60)
print("Processing complete!")
print("="*60)
