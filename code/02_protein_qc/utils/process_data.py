#!/usr/bin/env python3
"""
Utility functions for Olink data processing
"""
import numpy as np
import pandas as pd
from statistics import NormalDist


_NORMAL = NormalDist()

def norm_id(s):
    """
    Normalize ID column: convert to string, remove .0 suffix, strip whitespace.
    
    Parameters:
    -----------
    s : Series or array-like
        ID column to normalize
    
    Returns:
    --------
    Series
        Normalized ID column
    """
    return (s.astype(str)
             .str.replace(r'\.0$', '', regex=True)
             .str.strip())


def inverse_normal_transform(series):
    """Rank-based inverse normal transform, preserving missing values."""
    x = pd.to_numeric(series, errors="coerce")
    out = pd.Series(np.nan, index=series.index, dtype=float)
    mask = x.notna()
    n = int(mask.sum())
    if n == 0:
        return out
    if n == 1:
        out.loc[mask] = 0.0
        return out
    ranks = x.loc[mask].rank(method="average")
    probs = (ranks - 0.5) / n
    out.loc[mask] = probs.map(_NORMAL.inv_cdf).astype(float)
    return out


def inverse_normal_transform_columns(df, columns):
    """Apply rank-based inverse normal transform to selected columns."""
    df = df.copy()
    for col in columns:
        df[col] = inverse_normal_transform(df[col])
    return df


def mean_inverse_normalized_columns(df, columns):
    """Mean of rank-based inverse normal transformed columns, preserving row missingness."""
    if not columns:
        return pd.Series(np.nan, index=df.index, dtype=float)
    return inverse_normal_transform_columns(df[columns], columns).mean(axis=1)


def qc_olink(
    protein_olink: pd.DataFrame,
    sid,
    bad_flags=None,
):
    """Run QC on Olink data: LOD≥70%, warnings→NA, then mean_protein.

    Step order (correct order per Alexa et al. Nature Neuroscience):
    1. Detectability: LOD detection ≥ 70% (protein-level filter, on original data)
    2. Warnings → NA (Assay_Warning/QC_Warning in bad_flags)
    3. mean_protein: LOD_Flag PASS ≥ 90%, INT-transformed

    No 0 → NA: exact 0.0 is a real value on the signed NPX log scale, so it is kept.
    0 → NA is applied only in qc_nulisa. See issues/resolved/protein_qc_0na/KNOWN_ISSUES.md #1.

    Note: bad_flags is applied ONLY to __Assay_Warning / __QC_Warning columns
    (Step 3). __LOD_Flag is never matched against bad_flags — below-LOD
    individual datapoints are kept (Olink/Alexa et al. recommendation); LOD_Flag
    is consumed only for the Step 1 detection-rate fallback and the Step 4
    high-confidence filter for mean_protein.

    Returns (qc_df, steps) where steps is a list of (label, n_samples, n_proteins, ds, dp).
    """
    _steps = []
    _prev = {}
    def _record(label, df, na_count=None):
        n_s = df.shape[0]
        n_p = sum(c.endswith("__NPX") for c in df.columns)
        ds = '' if not _prev else n_s - _prev['s']
        dp = '' if not _prev else n_p - _prev['p']
        _steps.append((label, n_s, n_p, ds, dp, na_count))
        _prev['s'], _prev['p'] = n_s, n_p

    if bad_flags is None:
        bad_flags = {"WARN", "FAIL", "EXCLUDED", "MANUAL_WARN"}

    sid_set = set(sid)
    protein_olink = protein_olink[protein_olink["sid"].isin(sid_set)].copy()
    npx_cols = [col for col in protein_olink.columns if col.endswith("__NPX")]
    if npx_cols:
        protein_olink = protein_olink.loc[protein_olink[npx_cols].notna().any(axis=1)].copy()
    _record("Raw data", protein_olink)

    # Step 1: LOD-based detection rate (≥ 70%) - PROTEIN-LEVEL FILTER on original data
    keep_proteins = []
    for col in protein_olink.columns:
        if col.endswith("__NPX"):
            prefix = col.rsplit("__", 1)[0]
            lod_col = prefix + "__LOD"
            flag_col = prefix + "__LOD_Flag"

            detect_rate = None

            # Prefer numeric LOD when available
            if lod_col in protein_olink.columns:
                detect_rate = (
                    protein_olink[col].ge(protein_olink[lod_col]).sum()
                    / len(protein_olink)
                )
            # Fallback: use LOD_Flag when no numeric LOD is present
            elif flag_col in protein_olink.columns:
                flag_series = protein_olink[flag_col]
                valid_mask = flag_series.notna()
                if valid_mask.any():
                    good_mask = flag_series[valid_mask].eq("PASS")
                    detect_rate = good_mask.mean()

            if detect_rate is not None and detect_rate >= 0.7:
                keep_proteins.append(prefix)
    keep_set = set(keep_proteins)
    protein_olink = protein_olink[
        [col for col in protein_olink.columns if col == "sid" or col.rsplit("__", 1)[0] in keep_set]
    ].copy()
    _record("LOD detectability ≥70%", protein_olink)

    # No 0 → NA for Olink: exact 0.0 is a real value on the signed NPX log scale.
    # (Applied only in qc_nulisa.) See issues/resolved/protein_qc_0na/KNOWN_ISSUES.md #1.

    # Step 2: set NPX to NA where Assay_Warning or QC_Warning is bad (sample-level)
    npx_cols_before = [col for col in protein_olink.columns if col.endswith("__NPX")]
    na_before = sum(protein_olink[col].isna().sum() for col in npx_cols_before)
    for col in protein_olink.columns:
        if col.endswith("__NPX"):
            prefix = col.rsplit("__", 1)[0]
            warn_col = prefix + "__Assay_Warning"
            if warn_col in protein_olink.columns:
                mask = protein_olink[warn_col].isin(bad_flags)
                protein_olink.loc[mask, col] = np.nan
            qc_col = prefix + "__QC_Warning"
            if qc_col in protein_olink.columns:
                mask_qc = protein_olink[qc_col].isin(bad_flags)
                protein_olink.loc[mask_qc, col] = np.nan
    na_after = sum(protein_olink[col].isna().sum() for col in npx_cols_before)
    _record("Warning/Failed → NA", protein_olink, na_count=int(na_after - na_before))

    # Step 3: high-confidence proteins (detection ≥ 90%) → INT-based mean_protein
    npx_cols = [col for col in protein_olink.columns if col.endswith("__NPX")]
    high_conf_prefix = []
    for col in npx_cols:
        prefix = col.rsplit("__", 1)[0]
        flag_col = prefix + "__LOD_Flag"
        if flag_col in protein_olink.columns:
            flag_series = protein_olink[flag_col]
            valid_mask = flag_series.notna()
            if valid_mask.any() and flag_series[valid_mask].eq("PASS").mean() >= 0.9:
                high_conf_prefix.append(prefix)
    high_conf_npx = [prefix + "__NPX" for prefix in high_conf_prefix]
    protein_olink["mean_protein"] = mean_inverse_normalized_columns(protein_olink, high_conf_npx)
    _steps.append(('mean_protein ≥90%', _prev.get('s', ''), len(high_conf_prefix), '', '', None))
    return protein_olink, _steps


def high_confidence_npx(df, min_pass=0.9):
    """Return the `__NPX` columns whose detection rate >= `min_pass`.

    "High confidence" matches the Step-4 mean_protein membership in qc_olink:
    LOD_Flag PASS rate >= 0.9. Detection is computed from the numeric `__LOD`
    when present (NPX >= LOD over rows where both are non-missing) -- this is
    exactly the rule Olink/GWAS use to derive the PASS/WARN flag
    (NPX < LOD -> WARN), so the selected set is identical. Falls back to a
    pre-computed `__LOD_Flag` column when no numeric LOD is available.

    Unlike qc_olink (which INT-transforms before averaging), this only selects
    columns; the caller decides whether to INT.
    """
    cols = []
    for col in df.columns:
        if not col.endswith("__NPX"):
            continue
        prefix = col.rsplit("__", 1)[0]
        lod_col = prefix + "__LOD"
        flag_col = prefix + "__LOD_Flag"
        if lod_col in df.columns:
            valid = df[col].notna() & df[lod_col].notna()
            if valid.any() and (df.loc[valid, col] >= df.loc[valid, lod_col]).mean() >= min_pass:
                cols.append(col)
        elif flag_col in df.columns:
            flag_series = df[flag_col]
            valid = flag_series.notna()
            if valid.any() and flag_series[valid].eq("PASS").mean() >= min_pass:
                cols.append(col)
    return cols


def format_proteins_for_plink(df_qc, protein_cols, inverse_normal=True):
    """
    Format dataframe for PLINK .phe format with proteins and mean_protein.
    df_qc should already have FID, IID and mean_protein columns.
    protein_cols: explicit list of protein columns in df_qc (each ending '__NPX').
    The '__NPX' suffix is stripped from the output header.
    Returns: FID IID mean_protein prot1 prot2 ...
    """
    df_plink = df_qc.copy()

    # Check required columns exist
    if 'FID' not in df_plink.columns or 'IID' not in df_plink.columns:
        raise ValueError("FID and IID columns are missing! Check that ID mapping was applied.")
    if 'mean_protein' not in df_plink.columns:
        raise ValueError("mean_protein column is missing! Check QC function output.")

    # Normalize FID and IID
    df_plink['FID'] = norm_id(df_plink['FID'])
    df_plink['IID'] = norm_id(df_plink['IID'])

    if inverse_normal:
        df_plink = inverse_normal_transform_columns(df_plink, protein_cols)

    # Create proteins dataframe: FID IID mean_protein prot1 prot2 ...
    proteins_df = df_plink[['FID', 'IID', 'mean_protein'] + protein_cols].copy()

    # Remove __NPX suffix from column names
    rename_dict = {col: col.removesuffix('__NPX') for col in protein_cols}
    proteins_df = proteins_df.rename(columns=rename_dict)

    return proteins_df


def _assert_unique_iid(df, what):
    n_dup = int(df.duplicated(subset=['IID']).sum())
    if n_dup:
        raise ValueError(f"{n_dup} duplicate IID(s) in {what}; fix the sid -> IID mapping upstream")


def _assert_rows(df, n_expected, what):
    if df.shape[0] != n_expected:
        raise ValueError(f"row count changed in {what}: {n_expected} -> {df.shape[0]}")


def format_covariates_for_plink(df_mean_protein, df_pca=None, df_covar=None):
    """
    Format dataframe for PLINK .phe format with covariates only.
    df_mean_protein: FID, IID, and one or more mean_protein* columns.
    Returns: FID IID age sex pc1..pc10 mean_protein_*
    """
    df_plink = df_mean_protein.copy()
    
    # Normalize IID
    df_plink['IID'] = norm_id(df_plink['IID'])
    # One row per IID on every side, checked here and not repaired: a duplicate means
    # the sid -> IID mapping upstream is wrong, and keep='first' (as before) picked one
    # covariate row without saying which. A left merge on a duplicated key adds rows.
    _assert_unique_iid(df_plink, 'mean_protein')
    n_rows = df_plink.shape[0]
    
    # Merge PCA data
    if df_pca is not None:
        # Normalize IID before merge
        df_pca = df_pca.copy()   # a caller's slice: no SettingWithCopyWarning, no side effect
        df_pca['IID'] = norm_id(df_pca['IID'])
        _assert_unique_iid(df_pca, 'PCA')
        # Merge on IID
        df_plink = df_plink.merge(
            df_pca[['IID'] + [f'pc{i}' for i in range(1, 11)]],
            on='IID',
            how='left'
        )
        _assert_rows(df_plink, n_rows, 'PCA merge')
        pc_matched = df_plink[[f'pc{i}' for i in range(1, 11)]].notna().all(axis=1).sum()
        if pc_matched < df_plink.shape[0] * 0.9:
            print(f"WARNING: Only {pc_matched}/{df_plink.shape[0]} samples have PCA data!")

    # Merge covariates (age, sex, visit)
    if df_covar is not None:
        df_covar = df_covar.copy()
        df_covar['IID'] = norm_id(df_covar['IID'])
        _assert_unique_iid(df_covar, 'covar')
        covar_cols = ['IID', 'age', 'sex']
        if 'visit' in df_covar.columns:
            covar_cols.append('visit')
        df_plink = df_plink.merge(df_covar[covar_cols], on='IID', how='left')
        _assert_rows(df_plink, n_rows, 'covar merge')
        covar_matched = df_plink[['age', 'sex']].notna().all(axis=1).sum()
        if covar_matched < df_plink.shape[0] * 0.9:
            print(f"WARNING: Only {covar_matched}/{df_plink.shape[0]} samples have both age and sex!")
    else:
        print("WARNING: df_covar is None, age and sex will be NA")
        df_plink['age'] = pd.NA
        df_plink['sex'] = pd.NA

    # Convert sex to integer
    if 'sex' in df_plink.columns:
        df_plink['sex'] = pd.to_numeric(df_plink['sex'], errors='coerce').astype('Int64')

    mean_protein_cols = sorted(c for c in df_plink.columns if c.startswith('mean_protein'))
    fixed_cols = [c for c in ['age', 'sex', 'visit'] + [f'pc{i}' for i in range(1, 11)] if c in df_plink.columns]
    return df_plink[['FID', 'IID'] + fixed_cols + mean_protein_cols].copy()


def map_covar_to_genotype_iid(df_covar, id_map, clinical_id_col, genotype_id_col='IID', clinical_id_in_df='sid'):
    """
    Map covariates dataframe from clinical ID to genotype sample ID.
    
    Parameters:
    -----------
    df_covar : DataFrame
        Covariates dataframe with clinical_id_in_df column
    id_map : DataFrame
        ID mapping dataframe with columns [clinical_id_col, genotype_id_col]
    clinical_id_col : str
        Column name in id_map for clinical ID
    genotype_id_col : str
        Column name in id_map for genotype sample ID
    clinical_id_in_df : str
        Column name in df_covar for clinical ID (default: 'sid')
    
    Returns:
    --------
    df_mapped : DataFrame
        Covariates dataframe with IID column (genotype sample ID) instead of clinical ID
    """
    df_mapped = df_covar.merge(
        id_map,
        left_on=clinical_id_in_df,
        right_on=clinical_id_col,
        how='inner'
    )
    # Extract FID, IID, age, sex
    df_mapped = df_mapped[['FID', genotype_id_col, 'age', 'sex']].copy()
    df_mapped = df_mapped.rename(columns={genotype_id_col: 'IID'})

    return df_mapped


def add_nulisa2_lod(protein_nulisa, lod_df):
    """Join per-plate LOD from xlsx into protein_nulisa for CNS_nulisa2 panel.

    lod_df must have columns: TargetName, PlateID, TargetLOD_NPQ.
    Uses PlateID__CNS_nulisa2_CSF column in protein_nulisa to match rows.
    Proteins absent from lod_df (e.g. pTau217, Ab38) are left without LOD.
    """
    plate_col = 'PlateID__CNS_nulisa2_CSF'
    if plate_col not in protein_nulisa.columns:
        return protein_nulisa
    lod_pivot = lod_df.pivot_table(index='PlateID', columns='TargetName', values='TargetLOD_NPQ')
    new_lod = {}
    for npq_col in [c for c in protein_nulisa.columns if c.endswith('__CNS_nulisa2_CSF_NPQ')]:
        protein = npq_col[:-len('__CNS_nulisa2_CSF_NPQ')]
        if protein in lod_pivot.columns:
            new_lod[protein + '__CNS_nulisa2_CSF_LOD'] = protein_nulisa[plate_col].map(lod_pivot[protein])
    if new_lod:
        protein_nulisa = pd.concat(
            [protein_nulisa, pd.DataFrame(new_lod, index=protein_nulisa.index)], axis=1
        )
    return protein_nulisa


def qc_nulisa(protein_nulisa, sid, bad_flags=None):
    """Run QC on NuLISA data: SampleQC→NA, 0→NA, LOD detection ≥70%, then mean_protein.

    Step order (no missing rate filtering):
    1. SampleQC → NA (CNS panels only; Infl has none)
    2. 0 → NA (exact 0.0 in _NPQ columns)
    3. LOD detection ≥ 70% (computed after masking; proteins without LOD are KEPT)
    4. mean_protein: LOD detection ≥ 90% (proteins without LOD are EXCLUDED from
       mean_protein, matching qc_olink; they remain as individual GWAS phenotypes)

    Key differences from Olink:
    - Metric suffix is _NPQ (not __NPX)
    - LOD col: try prefix + '_LOD'; if missing try single-underscore variant (Infl panel)
    - SampleQC col: prefix + '_SampleQC' (CNS panels only; Infl has none)

    Returns (qc_df, steps) where steps is a list of (label, n_s, n_p, ds, dp, na_count).
    """
    _steps = []
    _prev = {}
    def _record(label, df, na_count=None):
        n_s = df.shape[0]
        n_p = sum(c.endswith('_NPQ') and '__' in c for c in df.columns)
        ds = '' if not _prev else n_s - _prev['s']
        dp = '' if not _prev else n_p - _prev['p']
        _steps.append((label, n_s, n_p, ds, dp, na_count))
        _prev['s'], _prev['p'] = n_s, n_p

    if bad_flags is None:
        # CNS_v1 panel uses {"PASS","WARN"}; CNS_nulisa2_CSF uses {"passed","warning"}
        # (different case AND different word form). Cover both with a superset, then
        # match case-insensitively via bad_flags_lc.
        bad_flags = {"WARN", "FAIL", "WARNING", "FAILED"}
    bad_flags_lc = {f.lower() for f in bad_flags}

    protein_nulisa = protein_nulisa[protein_nulisa['sid'].isin(set(sid))].copy()
    npq_cols = [c for c in protein_nulisa.columns if c.endswith('_NPQ') and '__' in c]
    if npq_cols:
        protein_nulisa = protein_nulisa.loc[protein_nulisa[npq_cols].notna().any(axis=1)].copy()
    _record("Raw data", protein_nulisa, na_count=None)

    def _lod_col(npq_col, cols):
        prefix = npq_col[:-4]  # strip '_NPQ'
        c1 = prefix + '_LOD'
        if c1 in cols:
            return c1
        c2 = prefix.replace('__', '_', 1) + '_LOD'
        if c2 in cols:
            return c2
        return None

    def _sampleqc_col(npq_col, cols):
        prefix = npq_col[:-4]
        c = prefix + '_SampleQC'
        return c if c in cols else None

    # Step 1: set NPQ → NA where SampleQC is bad (CNS panels; Infl has no SampleQC)
    npq_cols = [c for c in protein_nulisa.columns if c.endswith('_NPQ') and '__' in c]
    n_sampleqc_na = 0
    for npq_col in npq_cols:
        sqc = _sampleqc_col(npq_col, protein_nulisa.columns)
        if sqc:
            sqc_vals = protein_nulisa[sqc].astype("string").str.lower()
            mask = sqc_vals.isin(bad_flags_lc)
            n_sampleqc_na += int(mask.sum())
            protein_nulisa.loc[mask, npq_col] = np.nan
    _record("SampleQC → NA", protein_nulisa, na_count=n_sampleqc_na)

    # Step 2: replace exact 0.0 with NA in _NPQ columns.
    # NuLISA is the ONLY platform doing this: NPQ is strictly positive with a 2.39%
    # spike at exactly 0.0 then a gap to ~0.0248 — a below-detection placeholder, not
    # a measurement. Olink (signed log, 0 is real) and SomaLogic (no zeros) skip it.
    # See issues/resolved/protein_qc_0na/KNOWN_ISSUES.md #1.
    npq_cols = [c for c in protein_nulisa.columns if c.endswith('_NPQ') and '__' in c]
    n_zeros = int(sum((protein_nulisa[c] == 0.0).sum() for c in npq_cols))
    for npq_col in npq_cols:
        mask_zero = protein_nulisa[npq_col] == 0.0
        protein_nulisa.loc[mask_zero, npq_col] = np.nan
    _record("0 → NA", protein_nulisa, na_count=n_zeros)

    # Step 3: LOD detection rate >= 70% (computed after masking)
    npq_cols = [c for c in protein_nulisa.columns if c.endswith('_NPQ') and '__' in c]
    keep_npq = []
    # LOD columns of the kept proteins. The Infl panel names its LOD columns with
    # a single underscore ({assay}_{tag}_LOD), so they do not start with the
    # '{assay}__{tag}' prefix and the subsetting below would drop them -- leaving
    # Step 4 unable to resolve any LOD and silently returning an all-NaN
    # mean_protein. Collect them here, while the LOD is already resolved.
    keep_lod = set()
    no_lod_count = 0
    for npq_col in npq_cols:
        lod = _lod_col(npq_col, protein_nulisa.columns)
        if lod:
            rate = protein_nulisa[npq_col].ge(protein_nulisa[lod]).sum() / len(protein_nulisa)
            if rate >= 0.7:
                keep_npq.append(npq_col)
                keep_lod.add(lod)
        else:
            keep_npq.append(npq_col)  # no LOD available, keep
            no_lod_count += 1
    if no_lod_count:
        _steps.append((f'note: {no_lod_count} proteins no LOD', '', '', '', '', None))
    keep_prefixes = {c[:-4] for c in keep_npq}
    protein_nulisa = protein_nulisa[
        [c for c in protein_nulisa.columns
         if c == 'sid' or any(c.startswith(p) for p in keep_prefixes)
         or c in keep_lod
         or c == 'PlateID__CNS_nulisa2_CSF']
    ].copy()
    _record("LOD detectability ≥70%", protein_nulisa, na_count=None)

    # Step 4: INT-based mean_protein from high-confidence proteins (LOD detection >= 90%)
    # Proteins without a LOD (e.g. Ab42, pTau217) are excluded from mean_protein,
    # matching qc_olink (which only admits proteins with a __LOD_Flag). They still
    # pass through as individual GWAS phenotypes via Step 3; mean_protein membership
    # is independent of that. See protein_qc notes.
    npq_cols = [c for c in protein_nulisa.columns if c.endswith('_NPQ') and '__' in c]
    high_conf = []
    for npq_col in npq_cols:
        lod = _lod_col(npq_col, protein_nulisa.columns)
        if lod is None:
            continue  # no LOD → not high-confidence; exclude from mean_protein
        rate = protein_nulisa[npq_col].ge(protein_nulisa[lod]).sum() / len(protein_nulisa)
        if rate >= 0.9:
            high_conf.append(npq_col)
    # Fail loudly: an all-NaN mean_protein is a valid covariate column as far as
    # every downstream check is concerned, so it propagates silently into
    # covar.phe and then drops every sample at --glm time. If any protein here
    # carries a LOD at all, some of them must clear 90%.
    if not high_conf and any(_lod_col(c, protein_nulisa.columns) for c in npq_cols):
        raise ValueError(
            "mean_protein would be all-NaN: LOD columns are present but no protein "
            "reached 90% detection. Check LOD column naming for this panel "
            "(see _lod_col) before trusting it."
        )
    protein_nulisa['mean_protein'] = mean_inverse_normalized_columns(protein_nulisa, high_conf)
    _steps.append(('mean_protein ≥90%', _prev.get('s', ''), len(high_conf), '', '', None))
    return protein_nulisa, _steps


def qc_somalogic(protein_soma, sid, meta_df):
    """QC for SomaLogic ANML data: aptamer filter, then mean_protein.

    Steps:
    1. Filter aptamers: Type == 'Protein' AND ColCheck PASS in all batches
    2. Subset to cohort samples
    3. INT-based mean_protein from aptamers with non-missing >= 90% (no LOD available)

    No 0 → NA (ANML has no zeros; 0 → NA is applied only in qc_nulisa — see
    issues/resolved/protein_qc_0na/KNOWN_ISSUES.md #1) and no missing-rate / LOD step.
    Returns (qc_df, steps) where steps is a list of (label, n_s, n_p, ds, dp, na_count).
    """
    _steps = []
    _prev = {}
    def _record(label, df, na_count=None):
        n_s = df.shape[0]
        n_p = sum(c.endswith('__ANML') for c in df.columns)
        ds = '' if not _prev else n_s - _prev['s']
        dp = '' if not _prev else n_p - _prev['p']
        _steps.append((label, n_s, n_p, ds, dp, na_count))
        _prev['s'], _prev['p'] = n_s, n_p

    # Step 0: filter aptamers by metadata
    colcheck_cols = [c for c in meta_df.columns if c.startswith('ColCheck')]
    approved_seqids = meta_df[
        (meta_df['Type'] == 'Protein') &
        (meta_df['Organism'] == 'Human') &
        meta_df[colcheck_cols].eq('PASS').all(axis=1)
    ]['SeqId']
    # SeqId format: '10000-28' -> column key: 'seq.10000.28'
    approved_keys = {'seq.' + s.replace('-', '.') for s in approved_seqids}
    approved_cols = [c for c in protein_soma.columns
                     if c.endswith('__ANML') and c.split('__')[1] in approved_keys]
    protein_soma = protein_soma[['sid'] + approved_cols].copy()
    _steps.append(('aptamer filter', '', len(approved_cols), '', '', None))

    # Subset to cohort
    protein_soma = protein_soma[protein_soma['sid'].isin(set(sid))].copy()
    anml_cols = [c for c in protein_soma.columns if c.endswith('__ANML')]
    if anml_cols:
        protein_soma = protein_soma.loc[protein_soma[anml_cols].notna().any(axis=1)].copy()
    _record("Raw data", protein_soma, na_count=None)

    # No 0 → NA for SomaLogic: ANML has no zeros, and 0 → NA is applied only in
    # qc_nulisa. See issues/resolved/protein_qc_0na/KNOWN_ISSUES.md #1.

    # Step 1: INT-based mean_protein from aptamers with non-missing rate >= 90%
    anml_cols = [c for c in protein_soma.columns if c.endswith('__ANML')]
    high_conf = [c for c in anml_cols if protein_soma[c].notna().mean() >= 0.9]
    protein_soma['mean_protein'] = mean_inverse_normalized_columns(protein_soma, high_conf)
    _steps.append(('mean_protein ≥90%', _prev.get('s', ''), len(high_conf), '', '', None))
    return protein_soma, _steps
