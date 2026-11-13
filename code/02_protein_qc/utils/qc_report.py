"""Build one summary row per (cohort, platform, panel) for QC reporting.

- ``kept (sample x protein)``: rectangular dimensions of the post-QC matrix.
  Samples = those with ≥1 non-missing across the panel; proteins = post-LOD analytes.
  Individual cells inside the matrix can still be NA (warning→NA, SampleQC→NA, and
  0→NA for NuLISA only).
- ``gwas (sample x protein)``: genotyped merge count × GWAS-exported protein count
  after any .phe-only protein filters.
- ``n_aptamers_meta_pass``: Soma only — aptamers after meta filter (pre-cohort).

Dimension strings: start → after LOD protein drop (empty for Soma) → kept → gwas.
"""
from __future__ import annotations

from typing import Any

import pandas as pd


def _sx(n_s: int, n_p: int) -> str:
    return f"{n_s} x {n_p}"


def qc_summary_row(
    cohort: str,
    platform: str,
    panel: str,
    steps: list,
    qc_df,
    n_gwas: int,
    n_gwas_proteins: int | None = None,
) -> dict[str, Any]:
    """Align with Olink test table; NuLISA/Soma fill the same columns where applicable."""
    step_map = {s[0]: s for s in steps}

    row: dict[str, Any] = {
        "cohort": cohort,
        "platform": platform,
        "panel": panel,
        "start(sample x protein)": "",
        "remove_70LOD (sample x protein)": "",
        "kept (sample x protein)": "",
        "gwas (sample x protein)": "",
        "n_high_confident_protein": "",
        "n_aptamers_meta_pass": "",
    }

    if platform == "Olink":
        rd = step_map["Raw data"]
        lod = step_map["LOD detectability ≥70%"]
        row["start(sample x protein)"] = _sx(rd[1], rd[2])
        row["remove_70LOD (sample x protein)"] = _sx(lod[1], lod[2])
        row["n_high_confident_protein"] = int(step_map["mean_protein ≥90%"][2])
        pcols = [c for c in qc_df.columns if c.endswith("__NPX")]
    elif platform == "NuLISA":
        rd = step_map["Raw data"]
        lod = step_map["LOD detectability ≥70%"]
        row["start(sample x protein)"] = _sx(rd[1], rd[2])
        row["remove_70LOD (sample x protein)"] = _sx(lod[1], lod[2])
        row["n_high_confident_protein"] = int(step_map["mean_protein ≥90%"][2])
        pcols = [c for c in qc_df.columns if c.endswith("_NPQ") and "__" in c]
    elif platform == "SomaLogic":
        rd = step_map["Raw data"]
        row["start(sample x protein)"] = _sx(rd[1], rd[2])
        row["remove_70LOD (sample x protein)"] = ""
        apt = step_map.get("aptamer filter")
        if apt is not None and apt[2] != "":
            row["n_aptamers_meta_pass"] = int(apt[2])
        row["n_high_confident_protein"] = int(step_map["mean_protein ≥90%"][2])
        pcols = [c for c in qc_df.columns if c.endswith("__ANML")]
    else:
        raise ValueError(f"unknown platform {platform}")

    n_p = len(pcols)
    if n_gwas_proteins is None:
        n_gwas_proteins = n_p
    n_kept = int(qc_df[pcols].notna().any(axis=1).sum()) if pcols else 0
    row["kept (sample x protein)"] = _sx(n_kept, n_p)
    row["gwas (sample x protein)"] = _sx(n_gwas, n_gwas_proteins)

    return row
