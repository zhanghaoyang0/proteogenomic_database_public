#!/usr/bin/env python3
"""
Build the drug-target annotation for BioFINDER proteins from Open Targets
Platform 26.03 (info/druggable_genome/, downloaded from the URLs in url.txt).

Outputs:
    info/bf_protein/drug_targets.tsv.gz     one row per (ensembl_gene_id, drug_id)
    info/bf_protein/drug_targets_summary.tsv  coverage counts (metric, value)

Pairs come from drug_mechanism_of_action (ChEMBL mechanisms): a drug is listed
for a gene only when ChEMBL names that gene as its target. Genes are those in
info/gene_info.csv with n_bf_assay > 0.

drug_targets.tsv.gz columns:
    ensembl_gene_id, gene_symbol
    drug_id          ChEMBL id; the app links it to platform.opentargets.org/drug/
    drug_name, drug_type   from drug_molecule
    action_type      every ChEMBL action for the pair, sorted unique, "; "-joined
    mechanism        every ChEMBL mechanism for the pair, same aggregation
    max_stage        ordinal (0 Preclinical ... 9 Approved) of the highest
                     clinical_target.maxClinicalStage for THIS drug-target pair;
                     empty when Open Targets has no clinical record for the pair
    stage_label      display text for max_stage ("Not available" when empty)
    indications      clinical_target diseaseFromSource names for the pair,
                     sorted unique, "; "-joined

max_stage is deliberately never filled from drug_molecule.maximumClinicalStage:
that is the drug's stage across all its targets, so a drug approved through
target A would make target B look approved.
"""

from __future__ import annotations
import html
import re
import sys
from pathlib import Path

import pandas as pd
import pyarrow.dataset as ds

REPO = Path(__file__).resolve().parents[2]
OT_DIR = REPO / "info/druggable_genome"
GENE_INFO = REPO / "info/gene_info.csv"
OUT_TSV = REPO / "info/bf_protein/drug_targets.tsv.gz"
OUT_SUMMARY = REPO / "info/bf_protein/drug_targets_summary.tsv"

# Open Targets 26.03 codes maxClinicalStage as a string enum. Every code is
# mapped here on purpose: an unseen one stops the script instead of becoming NA.
STAGES = {
    "PRECLINICAL":   (0, "Preclinical"),
    "IND":           (1, "IND"),
    "EARLY_PHASE_1": (2, "Early Phase I"),
    "PHASE_1":       (3, "Phase I"),
    "PHASE_1_2":     (4, "Phase I/II"),
    "PHASE_2":       (5, "Phase II"),
    "PHASE_2_3":     (6, "Phase II/III"),
    "PHASE_3":       (7, "Phase III"),
    "PREAPPROVAL":   (8, "Pre-approval"),
    "APPROVAL":      (9, "Approved"),
    "UNKNOWN":       (None, "Not available"),
}
LABEL_OF = {rank: label for rank, label in STAGES.values() if rank is not None}
NA_LABEL = "Not available"


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def read(name: str, columns: list[str]) -> pd.DataFrame:
    path = OT_DIR / name
    if not path.exists():
        die(f"missing {path}")
    df = ds.dataset(path).to_table(columns=columns).to_pandas()
    if df.empty:
        die(f"{path} is empty")
    return df


def assert_unique(df: pd.DataFrame, key: list[str], what: str) -> None:
    n = int(df.duplicated(key).sum())
    if n:
        die(f"{what}: {n} duplicated {key} rows")


def join_all(values: pd.Series) -> str:
    return "; ".join(sorted({v.strip() for v in values.dropna() if v.strip()}))


# Trial registries spell one disease many ways ("alzheimer&#39;s disease",
# "alzheimer ’ s disease"); unescape and fold apostrophes and case so they collapse.
def clean_disease(x: str) -> str:
    x = html.unescape(x).lower()
    x = re.sub(r"\s*[’´`]\s*", "'", x)
    return re.sub(r"\s+", " ", x).strip()


# BioFINDER genes -------------------------------------------------------------
genes = pd.read_csv(GENE_INFO, usecols=["ensembl_gene_id", "gene_symbol", "n_bf_assay"])
genes = genes[genes["n_bf_assay"] > 0][["ensembl_gene_id", "gene_symbol"]]
assert_unique(genes, ["ensembl_gene_id"], "gene_info")
if genes.empty:
    die("no gene in gene_info.csv has n_bf_assay > 0")

# Mechanisms: (gene, drug) pairs ----------------------------------------------
moa = read("drug_mechanism_of_action", ["targets", "chemblIds", "actionType", "mechanismOfAction"])
moa = moa.explode("targets").explode("chemblIds")
moa = moa.rename(columns={"targets": "ensembl_gene_id", "chemblIds": "drug_id"})
moa = moa.dropna(subset=["ensembl_gene_id", "drug_id"])
moa = moa[moa["ensembl_gene_id"].isin(genes["ensembl_gene_id"])]
pairs = (moa.groupby(["ensembl_gene_id", "drug_id"], as_index=False)
            .agg(action_type=("actionType", join_all), mechanism=("mechanismOfAction", join_all)))
assert_unique(pairs, ["ensembl_gene_id", "drug_id"], "mechanism pairs")
if pairs.empty:
    die("no ChEMBL mechanism targets a BioFINDER gene")

# Clinical records per pair ---------------------------------------------------
ct = read("clinical_target", ["drugId", "targetId", "diseases", "maxClinicalStage"])
unmapped = set(ct["maxClinicalStage"].dropna()) - set(STAGES)
if unmapped:
    die(f"clinical_target has unmapped stage codes {sorted(unmapped)}; add them to STAGES")
ct["max_stage"] = ct["maxClinicalStage"].map(lambda s: STAGES[s][0] if s in STAGES else None)
ct["max_stage"] = pd.to_numeric(ct["max_stage"])
ct["disease"] = ct["diseases"].map(
    lambda d: [clean_disease(x["diseaseFromSource"]) for x in d if x and x.get("diseaseFromSource")] if d is not None else [])
ct = ct.rename(columns={"targetId": "ensembl_gene_id", "drugId": "drug_id"})
ct_stage = ct.groupby(["ensembl_gene_id", "drug_id"], as_index=False)["max_stage"].max()
ct_ind = (ct.explode("disease").groupby(["ensembl_gene_id", "drug_id"], as_index=False)
            .agg(indications=("disease", join_all)))
ct_pair = ct_stage.merge(ct_ind, on=["ensembl_gene_id", "drug_id"], how="inner", validate="1:1")
assert_unique(ct_pair, ["ensembl_gene_id", "drug_id"], "clinical_target pairs")

n = len(pairs)
out = pairs.merge(ct_pair, on=["ensembl_gene_id", "drug_id"], how="left", validate="1:1")
if len(out) != n:
    die(f"clinical join changed rows {n} -> {len(out)}")

# Drug names ------------------------------------------------------------------
mol = read("drug_molecule", ["id", "name", "drugType"])
mol = mol.rename(columns={"id": "drug_id", "name": "drug_name", "drugType": "drug_type"})
assert_unique(mol, ["drug_id"], "drug_molecule")
out = out.merge(mol, on="drug_id", how="left", validate="m:1")
if len(out) != n:
    die(f"drug_molecule join changed rows {n} -> {len(out)}")
out["drug_name"] = out["drug_name"].fillna(out["drug_id"])

out = out.merge(genes, on="ensembl_gene_id", how="left", validate="m:1")
if len(out) != n:
    die(f"gene join changed rows {n} -> {len(out)}")

out["max_stage"] = out["max_stage"].astype("Int64")
out["stage_label"] = out["max_stage"].map(lambda r: NA_LABEL if pd.isna(r) else LABEL_OF[int(r)])
out["indications"] = out["indications"].fillna("")
out = out[["ensembl_gene_id", "gene_symbol", "drug_id", "drug_name", "drug_type", "action_type",
           "mechanism", "max_stage", "stage_label", "indications"]]
out = out.sort_values(["gene_symbol", "max_stage", "drug_name"],
                      ascending=[True, False, True], na_position="last")
OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_TSV, sep="\t", index=False)

# Summary ---------------------------------------------------------------------
n_genes = len(genes)
n_with = out["ensembl_gene_id"].nunique()
rows = [
    ("n_biofinder_genes", n_genes),
    ("n_biofinder_genes_with_drugs", n_with),
    ("pct_biofinder_genes_with_drugs", round(100 * n_with / n_genes, 1)),
    ("n_unique_drugs", out["drug_id"].nunique()),
    ("n_drug_target_pairs", len(out)),
]
order = [label for _, label in sorted(STAGES.values(), key=lambda x: -1 if x[0] is None else x[0],
                                      reverse=True)]
by_stage = out["stage_label"].value_counts()
rows += [(f"n_drug_target_pairs_stage: {label}", int(by_stage.get(label, 0))) for label in order]
summary = pd.DataFrame([(m, str(v)) for m, v in rows], columns=["metric", "value"])
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

for path in (OUT_TSV, OUT_SUMMARY):
    if not path.exists() or path.stat().st_size == 0:
        die(f"{path} was not written")
print(summary.to_string(index=False))
print(f"wrote {OUT_TSV} ({len(out)} rows) and {OUT_SUMMARY}")
