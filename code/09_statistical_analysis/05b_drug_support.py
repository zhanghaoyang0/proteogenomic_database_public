#!/usr/bin/env python3
"""
Do protein-disease pairs with joint MR and colocalization support point at drug targets
already developed for that disease? (Open Targets Platform 26.03; Python because the
tables are parquet and the gwas R env has no arrow.)

Universe: every (gene, disease) pair that 05a could test, i.e. MR and colocalization both
run in at least one fluid (output/05a_disease_pair_evidence.csv.gz, written by 05a).
Supported: joint MR + colocalization in at least one fluid (05a's joint_sig), and per fluid.

Drug evidence for a pair: a ChEMBL mechanism names the gene as the drug's single-protein
target (protein-complex and family mechanisms would credit every member gene), and Open
Targets lists a clinical indication for that drug matching the disease. A FinnGen
endpoint is mapped to EFO/MONDO through the Open Targets FinnGen R12 study index
(FINNGEN_R12_<phenocode> -> diseaseIds); an indication matches only when it is the same
term. Direct parents matched too until 2026-10-30, when they credited IL6R with
ankylosing spondylitis through tocilizumab's approval for rheumatoid arthritis, a parent
term in EFO. Descendants matched until 2026-11-07: broad endpoints then took the drugs
of any sub-disease ("Autoimmune diseases" credited C5 with myasthenia gravis approvals
and CA12 with an MS trial; "Chronic diseases of tonsils and adenoids", mapped to
lymphatic system disease, credited FLT3 and IL6R with lymphoma trials). Stage = the
highest clinical stage over the matching drug-disease indications.

Outputs (output/):
  05b_drug_by_fluid_support.csv  matching drug evidence by mutually exclusive
                                  CSF-only, plasma-only and shared support classes
  05b_drug_supported.csv   supported pairs whose gene is a drug target, with the drugs
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

REPO = Path(__file__).resolve().parents[3]
PROJ = REPO
OUT = PROJ / "output"
OT_INFO = REPO / "info/druggable_genome"
OT_REF = Path.home() / "db/shared_data/ref_geno/opentargets_26.03"   # study, disease
PAIRS = OUT / "05a_disease_pair_evidence.csv.gz"
KEPT = OUT / "00_kept_assays.csv"

STAGE_RANK = {"PRECLINICAL": 0, "IND": 1, "EARLY_PHASE_1": 2, "PHASE_1": 3, "PHASE_1_2": 4,
              "PHASE_2": 5, "PHASE_2_3": 6, "PHASE_3": 7, "PREAPPROVAL": 8, "APPROVAL": 9,
              "UNKNOWN": -1}
CLINICAL = 3    # Phase I or later
APPROVED = 9


def die(msg: str) -> None:
    sys.exit(f"ERROR: {msg}")


def read(path: Path, columns: list[str]) -> pd.DataFrame:
    if not path.exists():
        die(f"missing {path}")
    return ds.dataset(path).to_table(columns=columns).to_pandas()


# ---- tested pairs ----------------------------------------------------------------
if not PAIRS.exists():
    die(f"missing {PAIRS}; run 05a_disease_data.R first")
ev = pd.read_csv(PAIRS, usecols=["protein_assay", "outcome_id", "disease", "disease_category",
                                 "fluid", "mr_ok", "coloc_tested", "joint_sig", "ivw_b"])
kept = pd.read_csv(KEPT, usecols=["protein_assay", "ensembl_gene_id", "gene"])
ev = ev.merge(kept, on="protein_assay", how="inner", validate="many_to_one")
ev = ev[ev["mr_ok"].astype(bool) & ev["coloc_tested"].astype(bool)]
# A complex assay ("ENSGa,ENSGb") is a test of each subunit gene.
ev["ensembl_gene_id"] = ev["ensembl_gene_id"].str.split(",")
ev = ev.explode("ensembl_gene_id")
ev["joint_sig"] = ev["joint_sig"].astype(bool)
for f in ["CSF", "Plasma"]:
    ev[f"tested_{f.lower()}"] = ev["fluid"] == f
    ev[f"supported_{f.lower()}"] = ev["joint_sig"] & (ev["fluid"] == f)
key = ["ensembl_gene_id", "outcome_id"]
pairs = (ev.groupby(key + ["gene", "disease", "disease_category"], as_index=False)
           .agg(supported=("joint_sig", "any"),
                supported_csf=("supported_csf", "any"), supported_plasma=("supported_plasma", "any"),
                tested_csf=("tested_csf", "any"), tested_plasma=("tested_plasma", "any")))
if pairs.duplicated(key).any():
    die("tested pairs not unique by gene x outcome")
print(f"tested gene x disease pairs: {len(pairs):,}; supported: {pairs.supported.sum():,}")

# ---- FinnGen endpoint -> disease ontology terms ---------------------------------
study = read(OT_REF / "study", ["studyId", "diseaseIds"])
study = study[study["studyId"].str.startswith("FINNGEN_R12_")]
study["phenocode"] = study["studyId"].str.replace("FINNGEN_R12_", "", regex=False)
outc = pd.DataFrame({"outcome_id": pairs["outcome_id"].unique()})
outc["phenocode"] = (outc["outcome_id"].str.replace(r"^finngen_R12_", "", regex=True)
                                       .str.replace(r"_meta_.*$", "", regex=True))
outc = outc.merge(study[["phenocode", "diseaseIds"]], on="phenocode", how="left")
n_map = outc["diseaseIds"].map(lambda x: not isinstance(x, float) and x is not None and len(x) > 0).sum()
print(f"outcomes mapped to an ontology term: {n_map} of {len(outc)}")
if n_map < 0.8 * len(outc):
    die("fewer than 80% of outcomes map to Open Targets FinnGen R12 studies")
# The endpoint's own terms only, no descendants (see the docstring).
outc_terms = (outc.explode("diseaseIds").dropna(subset=["diseaseIds"])
                  .rename(columns={"diseaseIds": "term"})[["outcome_id", "term"]].drop_duplicates())

# ---- drug -> gene (single protein) and drug -> indication -----------------------
moa = read(OT_INFO / "drug_mechanism_of_action", ["targets", "chemblIds", "targetType", "actionType"])
moa = moa[moa["targetType"] == "single protein"]
moa = moa.explode("targets").explode("chemblIds").rename(
    columns={"targets": "ensembl_gene_id", "chemblIds": "drug_id", "actionType": "action"})
moa = moa.dropna(subset=["ensembl_gene_id", "drug_id"])
moa = moa.groupby(["ensembl_gene_id", "drug_id"], as_index=False)["action"].agg(
    lambda s: "; ".join(sorted(set(s.dropna()))))
ind = read(OT_INFO / "clinical_indication", ["drugId", "diseaseId", "maxClinicalStage"])
bad = set(ind["maxClinicalStage"]) - set(STAGE_RANK)
if bad:
    die(f"unknown clinical stage codes: {bad}")
ind["stage"] = ind["maxClinicalStage"].map(STAGE_RANK)
ind = ind.rename(columns={"drugId": "drug_id", "diseaseId": "term"})
names = read(OT_INFO / "drug_molecule", ["id", "name"]).rename(columns={"id": "drug_id", "name": "drug_name"})

# gene x outcome x drug, with the best stage over matching indications
gd = (pairs[key].merge(moa, on="ensembl_gene_id")
                .merge(outc_terms, on="outcome_id")
                .merge(ind[["drug_id", "term", "stage"]], on=["drug_id", "term"]))
gd = gd.groupby(key + ["drug_id", "action"], as_index=False)["stage"].max()
best = gd.groupby(key, as_index=False)["stage"].max().rename(columns={"stage": "drug_stage"})
pairs = pairs.merge(best, on=key, how="left")
pairs["drug_stage"] = pairs["drug_stage"].fillna(-1).astype(int)
pairs["is_target"] = pairs["ensembl_gene_id"].isin(moa["ensembl_gene_id"])

# ---- support classes for the figure ---------------------------------------------
support_class = np.select(
    [
        pairs["supported_csf"] & ~pairs["supported_plasma"],
        ~pairs["supported_csf"] & pairs["supported_plasma"],
        pairs["supported_csf"] & pairs["supported_plasma"],
    ],
    ["CSF-supported only", "Plasma-supported only", "Supported in both"],
    default="",
)
by_fluid = pairs.assign(support_class=support_class)
by_fluid = by_fluid[by_fluid["support_class"] != ""].copy()
if len(by_fluid) != int(pairs["supported"].sum()):
    die("fluid-support classes do not partition supported pairs")

fluid_rows = []
for support in ["CSF-supported only", "Plasma-supported only", "Supported in both"]:
    sub = by_fluid[by_fluid["support_class"] == support]
    if sub.empty:
        die(f"no supported pairs in {support}")
    for label, cut in [("Clinical-stage drug", CLINICAL), ("Approved drug", APPROVED)]:
        n_with_drug = int((sub["drug_stage"] >= cut).sum())
        fluid_rows.append(dict(
            support_class=support,
            drug_evidence=label,
            n_pairs=len(sub),
            n_with_drug=n_with_drug,
            percent=100 * n_with_drug / len(sub),
        ))
by_fluid_out = pd.DataFrame(fluid_rows)
by_fluid_out.to_csv(OUT / "05b_drug_by_fluid_support.csv", index=False)

# ---- supported pairs on drug targets ----------------------------------------------
# For the table: drugs whose indication matches the disease, and, for the rest, whether
# the gene has an approved drug at all (a repurposing lead rather than a known use).
direction = (ev[ev["joint_sig"]].groupby(key)["ivw_b"]
               .agg(lambda b: "+" if (b > 0).all() else "-" if (b < 0).all() else "mixed")
               .rename("mr_direction").reset_index())
appr_any = (moa.merge(ind[ind["stage"] == APPROVED][["drug_id"]].drop_duplicates(), on="drug_id")
               .merge(names, on="drug_id").groupby("ensembl_gene_id")["drug_name"]
               .agg(lambda s: "; ".join(sorted(set(s))[:5])).rename("approved_drugs_any_indication"))
STAGE_NAME = {v: k for k, v in STAGE_RANK.items()}
# Stage, then name: sorting on stage alone left ties in merge order, so the listed drugs
# (and the first one, which the Results quote) changed whenever an upstream table did.
gdn = gd.merge(names, on="drug_id").sort_values(["stage", "drug_name"], ascending=[False, True],
                                                 kind="mergesort")
gdn["label"] = gdn["drug_name"] + " (" + gdn["action"].str.lower() + ", " + gdn["stage"].map(STAGE_NAME) + ")"
match = (gdn.groupby(key)["label"].agg(lambda s: "; ".join(list(s)[:5]))
            .rename("matching_drugs").reset_index())
sup = (pairs[pairs["supported"] & pairs["is_target"]]
         .merge(direction, on=key, how="left").merge(match, on=key, how="left")
         .merge(appr_any, on="ensembl_gene_id", how="left"))
sup = sup.sort_values(["drug_stage", "gene"], ascending=[False, True])
sup.to_csv(OUT / "05b_drug_supported.csv", index=False)
print(f"supported pairs on a drug target: {len(sup)}; with a matching indication: "
      f"{(sup.drug_stage >= CLINICAL).sum()} (approved {(sup.drug_stage == APPROVED).sum()})")
print("05b_drug_support.py done")
