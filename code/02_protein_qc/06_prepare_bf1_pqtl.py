#!/usr/bin/env python3
"""Prepare selected BF1 post-QC protein panels for pQTL GWAS.

This reuses the existing post-QC, baseline protein matrices. It intersects them
with the final BF1 genotype sample set and retains proteins with at least the
requested number of model-complete samples (protein plus the covariates used by
01_gwas_protein.sh). No protein transformation is applied here; PLINK2 performs
the rank-based inverse-normal transformation at GWAS time.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import pandas as pd


BASE_COVARS = ["age", "sex"] + [f"pc{i}" for i in range(1, 11)]


def norm_id(series: pd.Series) -> pd.Series:
    out = series.astype("string").str.strip()
    return out.str.replace(r"\.0$", "", regex=True)


def read_manifest(path: Path) -> list[str]:
    manifest = pd.read_csv(path, sep="\t", dtype="string", comment="#")
    if "panel" not in manifest.columns:
        raise ValueError(f"Manifest must contain a panel column: {path}")
    if "cohort" in manifest.columns:
        invalid = manifest.loc[manifest["cohort"].notna() & manifest["cohort"].ne("bf1"), "cohort"]
        if not invalid.empty:
            raise ValueError(f"This preparer accepts BF1 only; found: {sorted(invalid.unique())}")
    panels = manifest["panel"].dropna().str.strip()
    panels = [panel for panel in panels if panel]
    if not panels:
        raise ValueError(f"No panels in manifest: {path}")
    if len(panels) != len(set(panels)):
        raise ValueError(f"Duplicate panels in manifest: {path}")
    return panels


def read_psam_iids(path: Path) -> set[str]:
    psam = pd.read_csv(path, sep=r"\s+", dtype="string")
    psam.columns = [column.lstrip("#") for column in psam.columns]
    if "IID" not in psam.columns:
        raise ValueError(f"Missing IID column: {path}")
    return set(norm_id(psam["IID"]).dropna())


def atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    os.close(handle)
    temp_path = Path(temp_name)
    try:
        frame.to_csv(temp_path, sep="\t", index=False, na_rep="NA")
        os.replace(temp_path, path)
    finally:
        temp_path.unlink(missing_ok=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel-manifest", type=Path, required=True)
    parser.add_argument("--min-model-n", type=int, default=500)
    parser.add_argument(
        "--qc-dir",
        type=Path,
        default=Path.home() / "db/secure_data/bf/protein/qc",
    )
    parser.add_argument(
        "--id-map",
        type=Path,
        default=Path.home() / "db/secure_data/bf/pheno/clean/bf1_sid_iid.csv",
    )
    parser.add_argument(
        "--psam",
        type=Path,
        default=Path.home() / "db/secure_data/bf/array/bf1/imputed_qc/bf1_hg19.psam",
    )
    parser.add_argument(
        "--covar",
        type=Path,
        default=Path.home() / "db/secure_data/bf/array/gwas_input/protein/bf1_covar.phe",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path.home() / "db/secure_data/bf/array/gwas_input/protein",
    )
    parser.add_argument("--summary-out", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.min_model_n <= 0:
        raise ValueError("--min-model-n must be positive")
    for path in (args.panel_manifest, args.qc_dir, args.id_map, args.psam, args.covar):
        if not path.exists():
            raise FileNotFoundError(path)

    panels = read_manifest(args.panel_manifest)
    genotype_iids = read_psam_iids(args.psam)

    id_map = pd.read_csv(args.id_map, dtype="string", usecols=["sid", "FID", "IID"])
    for column in ("sid", "FID", "IID"):
        id_map[column] = norm_id(id_map[column])
    valid = (
        id_map["sid"].notna()
        & id_map["IID"].notna()
        & id_map["FID"].notna()
        & id_map["IID"].ne("Removed in QC")
        & id_map["FID"].ne("Removed in QC")
        & id_map["IID"].isin(genotype_iids)
    )
    id_map = id_map.loc[valid, ["sid", "IID"]].drop_duplicates()
    if id_map["IID"].duplicated().any():
        raise ValueError("Multiple BF1 sid values map to the same final genotype IID")

    covar = pd.read_csv(args.covar, sep="\t", dtype="string", na_values=["NA", "nan", ""])
    covar["IID"] = norm_id(covar["IID"])
    missing_base = [column for column in BASE_COVARS if column not in covar.columns]
    if missing_base:
        raise ValueError(f"Missing required covariates: {missing_base}")
    covar = covar.drop_duplicates("IID")

    summary_rows: list[dict[str, object]] = []
    for panel in panels:
        matrix_path = args.qc_dir / f"bf1_{panel}_protein.csv"
        if not matrix_path.exists():
            raise FileNotFoundError(matrix_path)
        output_path = args.out_dir / f"bf1_{panel}.phe"
        if output_path.exists() and not args.overwrite:
            raise FileExistsError(f"Refusing to overwrite {output_path}; pass --overwrite to replace it")

        matrix = pd.read_csv(matrix_path, low_memory=False)
        if "sid" not in matrix.columns:
            raise ValueError(f"Missing sid column: {matrix_path}")
        matrix["sid"] = norm_id(matrix["sid"])
        protein_columns = [column for column in matrix.columns if column != "sid"]
        if not protein_columns:
            raise ValueError(f"No protein columns: {matrix_path}")

        analysis = matrix.merge(id_map, how="inner", on="sid").drop_duplicates("IID")
        analysis = analysis.merge(covar, how="left", on="IID", suffixes=("", "_covar"))
        panel_mean = f"mean_protein_{panel}"
        if panel_mean not in analysis.columns:
            raise ValueError(f"Missing panel mean covariate {panel_mean} in {args.covar}")
        model_covars = BASE_COVARS + [panel_mean]
        model_complete = analysis[model_covars].notna().all(axis=1)

        model_n = {
            protein: int((analysis[protein].notna() & model_complete).sum())
            for protein in protein_columns
        }
        selected = [protein for protein in protein_columns if model_n[protein] >= args.min_model_n]
        if not selected:
            raise RuntimeError(f"No proteins meet N >= {args.min_model_n} in {panel}")

        phe = analysis[["IID"] + selected].copy()
        phe.insert(0, "FID", "0")
        atomic_tsv(phe, output_path)

        selected_n = [model_n[protein] for protein in selected]
        summary_rows.append(
            {
                "cohort": "bf1",
                "panel": panel,
                "proteins_available": len(protein_columns),
                "proteins_written": len(selected),
                "N_model_min": min(selected_n),
                "N_model_median": float(pd.Series(selected_n).median()),
                "N_model_max": max(selected_n),
                "phe_path": str(output_path),
            }
        )
        print(
            f"{panel}: wrote {len(selected)}/{len(protein_columns)} proteins "
            f"(model N {min(selected_n)}-{max(selected_n)}) -> {output_path}"
        )

    summary = pd.DataFrame(summary_rows)
    if args.summary_out:
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.summary_out, sep="\t", index=False)
        print(f"Wrote aggregate summary: {args.summary_out}")


if __name__ == "__main__":
    main()
