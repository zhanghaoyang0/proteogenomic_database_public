"""Build assay -> UniProt mapping from panel target lists.

Input:
  - Olink/NULISA panel files in ~/db/secure_data/bf/pheno/raw/protein_info/
  - SomaLogic metadata in ~/db/secure_data/bf/pheno/raw/haoyang__20260429_234900_somalogic/
  - info/bf_protein/uniprot/rename.csv (manual Olink display name -> UniProt)

Output:
  - info/bf_protein/uniprot/assay_uniprot.csv (cols: tag, gene, uniprot_id, source)
  - Stdout: coverage report (tags covered vs missing).

Run:
  python3 code/02_protein_qc/04_annot_uniprot.py
  python3 code/02_protein_qc/04_annot_uniprot.py --merge-only   # assay_annot.csv only, no raw lists needed
"""
from pathlib import Path
import re
import sys

import pandas as pd

PROJ = Path(__file__).resolve().parents[2]
RAW = Path.home() / "db/secure_data/bf/pheno/raw"
INFO = PROJ / "info/bf_protein"
ENSG_DIR = INFO / "ensg"
UNIPROT_DIR = INFO / "uniprot"
PANEL_DIR = RAW / "protein_info"
SOMA_META = RAW / "haoyang__20260429_234900_somalogic/plasma_somalogic4.1_metadata_Jan2024.csv"
FEATHER = RAW / "haoyang__20260430_111026_olink/haoyang__20260430_111026.feather"
RENAME_CSV = UNIPROT_DIR / "rename.csv"
OUT = UNIPROT_DIR / "assay_uniprot.csv"
UNMATCHED_OUT = UNIPROT_DIR / "unmatched.csv"
ENSG_IN = ENSG_DIR / "assay_ensg.csv"
ANNOT_OUT = INFO / "assay_annot.csv"

def norm_name(name) -> str:
    """Collapse a protein display name to a joinable key.

    Panel target lists use display names ("BMP-7", "IL-10", "LDL_receptor")
    while 03_annot_ensg.R emits HGNC symbols ("BMP7", "IL10", "LDLR"), so the
    (tag, gene) join below misses every hyphenated NULISA target unless both
    sides are normalized first.
    """
    return re.sub(r"[^A-Z0-9]", "", str(name).upper())


def load_rename_map():
    """uniprot/rename.csv as {normalized name: (gene, uniprot_id)}.

    Hand-curated aliases that no panel list resolves. Used both for the custom
    Olink panels (feather column names) and as the last fallback of the
    ENSG x UniProt merge, so one file covers every alias in the project.
    """
    if not RENAME_CSV.exists():
        return {}
    return {norm_name(r.display_name): (r.gene, r.uniprot_id)
            for r in pd.read_csv(RENAME_CSV).itertuples()}


# panel_file -> list of tags it covers. CardiovascularIII_o4csf is a separate
# plate batch but a subset of the Target 96 CVD III 92-plex, so we reuse the
# same panel CSV for the alias tag.
PANEL_TAGS = {
    "Olink Target 96 Inflammation.csv":      ["INFLAMMATION_v3004_o1csf", "INFLAMMATION_v3004_o1plasma"],
    "Olink Target 96 Neurology.csv":         ["NEUROLOGY_v8001_o1csf", "NEUROLOGY_v8001_o1plasma"],
    "Olink Target 96 Cardiovascular III.csv":["CARDIOVASCULAR_III_v6101_o1csf", "CARDIOVASCULAR_III_v6101_o1plasma", "CardiovascularIII_o4csf"],
    "Olink Target 96 Neuro Exploratory.csv": ["NEURO_EXPLORATORY_v3911_o1csf", "NEURO_EXPLORATORY_v3911_o1plasma"],
    "Olink Explore Inflammation.csv":        ["Inflammation_o2csf", "Inflammation_o3plasma"],
    "Olink Explore Inflammation II.csv":     ["Inflammation_II_o2csf"],
    "Olink Explore Cardiometabolic.csv":     ["Cardiometabolic_o2csf", "Cardiometabolic_o3plasma"],
    "Olink Explore Cardiometabolic II.csv":  ["Cardiometabolic_II_o2csf"],
    "Olink Explore Neurology.csv":           ["Neurology_o2csf", "Neurology_o3plasma"],
    "Olink Explore Neurology II.csv":        ["Neurology_II_o2csf"],
    "Olink Explore Oncology.csv":            ["Oncology_o2csf", "Oncology_o3plasma"],
    "Olink Explore Oncology II.csv":         ["Oncology_II_o2csf"],
    "olink-explore-ht-assay-list.xlsx":      ["ExploreHT_o6csf", "ExploreHT_o7plasma"],
}

# nulisa.csv has multi-panel rows (Panel col may list several panels comma-joined);
# we explode by panel name then route to BF tags.
NULISA_PANEL_TAGS = {
    "CNS Panel 120":    ["CNS_v1_nulisa_CSF", "CNS_v1_nulisa_Plasma", "CNS_nulisa2_CSF"],
    "Inflammation 250": ["Infl_nulisa_CSF"],
}

# Custom Olink panels (not in public catalog): genes extracted from feather
# column headers; UniProt resolved via gene-match against the standard Olink
# panels + uniprot/rename.csv fallback for Olink display-name variants.
CUSTOM_PANEL_TAGS = [
    "Microglia_o5csf",
    "Biogen_Alz_Panel_v100010101_o1csf",
    "Biogen_Alz_Single_v100010101_o1csf",
]


def parse_olink_csv(path: Path) -> pd.DataFrame:
    # semicolon, UTF-8 BOM, cols: UniProt ID;Gene;Protein name
    df = pd.read_csv(path, sep=";", encoding="utf-8-sig")
    df = df.rename(columns={"UniProt ID": "uniprot_id", "Gene": "gene"})
    return df[["gene", "uniprot_id"]].dropna()


def parse_explore_ht(path: Path) -> pd.DataFrame:
    # header at row 11; cols: OlinkID, UniProt, ProteinName, Gene
    df = pd.read_excel(path, skiprows=10, header=None,
                       names=["olink_id", "uniprot_id", "protein_name", "gene"])
    return df[["gene", "uniprot_id"]].dropna()


def parse_nulisa_csv(path: Path) -> pd.DataFrame:
    # cols: Panel, Uniprot ID, Target List, Hallmarks..., Function
    # Panel cell may contain comma-joined list (e.g. "CNS Panel 120, Neuro 220 Panel").
    # Target List uses NULISA common name (e.g. "Abeta38 (ab38)" or "ACHE").
    df = pd.read_csv(path)
    df = df.rename(columns={"Uniprot ID": "uniprot_id", "Target List": "gene", "Panel": "panel"})
    # Strip "(alias)" suffix like "Abeta38 (ab38)" -> "Abeta38"
    df["gene"] = df["gene"].str.replace(r"\s*\(.*\)\s*$", "", regex=True)
    return df[["panel", "gene", "uniprot_id"]].dropna(subset=["gene", "uniprot_id"])


def parse_soma(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={"EntrezGeneSymbol": "gene", "UniProt": "uniprot_id"})
    # Some aptamers have non-human or composite UniProt strings; keep as-is, downstream handles.
    return df[["gene", "uniprot_id"]].dropna()


def main():
    rows = []
    for fname, tags in PANEL_TAGS.items():
        path = PANEL_DIR / fname
        if not path.exists():
            print(f"[skip] missing file: {fname}")
            continue
        if fname.endswith(".csv"):
            df = parse_olink_csv(path)
        elif fname.endswith(".xlsx"):
            df = parse_explore_ht(path)
        else:
            continue
        for tag in tags:
            sub = df.assign(tag=tag, source=fname)
            rows.append(sub[["tag", "gene", "uniprot_id", "source"]])
        print(f"[ok ] {fname}: {len(df)} entries -> {tags}")

    # NULISA: one csv with multi-panel rows
    nulisa_path = PANEL_DIR / "nulisa.csv"
    if nulisa_path.exists():
        nul = parse_nulisa_csv(nulisa_path)
        for panel_name, tags in NULISA_PANEL_TAGS.items():
            sub_panel = nul[nul["panel"].str.contains(panel_name, regex=False, na=False)]
            for tag in tags:
                sub = sub_panel.assign(tag=tag, source=f"{nulisa_path.name}:{panel_name}")
                rows.append(sub[["tag", "gene", "uniprot_id", "source"]])
            print(f"[ok ] nulisa.csv / {panel_name}: {len(sub_panel)} entries -> {tags}")
    else:
        print(f"[skip] missing nulisa.csv")

    # Custom Olink panels: pull gene list from feather column headers
    # ("{gene}__{tag}__NPX"), then resolve each gene to UniProt via the
    # standard-panel master + uniprot/rename.csv fallback for display-name
    # variants (e.g. "LDL_receptor" -> "LDLR / P01130").
    if FEATHER.exists() and RENAME_CSV.exists():
        olink_rows = [r for r in rows if r["source"].str.contains("Olink", na=False).any()]
        g2u = dict(zip(pd.concat(olink_rows)["gene"], pd.concat(olink_rows)["uniprot_id"]))
        rename_map = load_rename_map()
        feather_cols = pd.read_feather(FEATHER, columns=None).columns

        for tag in CUSTOM_PANEL_TAGS:
            suffix = f"__{tag}__NPX"
            genes = sorted({c.split("__")[0] for c in feather_cols if c.endswith(suffix)})
            resolved, unresolved = [], []
            for g in genes:
                if g in g2u:
                    resolved.append({"gene": g, "uniprot_id": g2u[g]})
                elif norm_name(g) in rename_map:
                    gene, up = rename_map[norm_name(g)]
                    resolved.append({"gene": gene, "uniprot_id": up})
                else:
                    unresolved.append(g)
            if unresolved:
                raise SystemExit(f"[FAIL] {tag}: unresolved genes {unresolved} -- "
                                 f"add to {RENAME_CSV.name}")
            df = pd.DataFrame(resolved).assign(tag=tag, source=f"feather:{tag}")
            rows.append(df[["tag", "gene", "uniprot_id", "source"]])
            print(f"[ok ] {tag}: {len(resolved)}/{len(genes)} genes resolved")
    else:
        print(f"[skip] feather or rename CSV missing; custom panels not added")

    if SOMA_META.exists():
        df = parse_soma(SOMA_META)
        df = df.assign(tag="somalogic_plasma", source=SOMA_META.name)
        rows.append(df[["tag", "gene", "uniprot_id", "source"]])
        print(f"[ok ] {SOMA_META.name}: {len(df)} entries -> ['somalogic_plasma']")
    else:
        print(f"[skip] missing SomaLogic metadata: {SOMA_META}")

    out = pd.concat(rows, ignore_index=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)

    # uniprot/unmatched.csv = rows where the panel source did not provide a
    # UniProt ID. In practice this is empty for Olink/NULISA; SomaLogic
    # aptamers occasionally have a missing UniProt field.
    missing_up = out["uniprot_id"].isna() | (out["uniprot_id"].astype(str).str.len() == 0)
    out[missing_up][["tag", "gene", "source"]].drop_duplicates().to_csv(UNMATCHED_OUT, index=False)
    out = out[~missing_up].drop_duplicates(subset=["tag", "gene", "uniprot_id"])

    # This used to also duplicate every ExploreHT row under the tags the
    # proteogenomic_network bridged matrices used (olink_explore3072_csf,
    # olink_exploreHT_csf, olink_exploreHT_plasma), so the merge below would
    # reach them. 03_annot_ensg.R no longer emits those tags, so the aliases
    # joined against nothing -- 16,239 rows that could never match.
    out_for_merge = out.copy()
    out_for_merge["protein_assay"] = (
        out_for_merge["gene"].astype(str) + "__" + out_for_merge["tag"]
    )

    out_for_merge.drop(columns=["protein_assay"]).drop_duplicates().to_csv(
        OUT, index=False
    )
    print(f"\nUnmatched (no uniprot_id): {int(missing_up.sum()):,} -> {UNMATCHED_OUT}")
    print(f"Wrote {len(out_for_merge):,} rows -> {OUT}")

    # Coverage report against the post-QC panel set (one matrix per cohort x panel).
    # protein/qc/ is the broadest panel list; .phe drops panels with no proteins
    # passing the N>1000 GWAS filter.
    qc_dir = Path.home() / "db/secure_data/bf/protein/qc"
    qc_tags = sorted({p.stem.split("_", 1)[1].removesuffix("_protein")
                      for p in qc_dir.glob("bf*_*_protein.csv")})
    covered = set(out_for_merge["tag"].unique())
    print("\n=== Coverage vs post-QC panels ===")
    for t in qc_tags:
        mark = "OK " if t in covered else "MISS"
        n_up = (out_for_merge["tag"] == t).sum() if t in covered else 0
        print(f"  [{mark}] {t:45s} uniprot_rows={n_up:>5d}")
    print(f"\nCovered: {len(covered & set(qc_tags))}/{len(qc_tags)}")

    merge_annot(out_for_merge)


# One UniProt accession, or several joined by "_" as the panel lists write complexes.
UNIPROT_ACC = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
UNIPROT_RE = re.compile(rf"{UNIPROT_ACC}(?:_{UNIPROT_ACC})*")


def merge_annot(out_for_merge: pd.DataFrame) -> None:
    """Merge ENSG + UniProt into the master assay annotation table.

    Each step fills only the rows still missing a UniProt ID:
      1. rename.csv by assay    -- hand-curated, so it beats the panel lists: names
                                   a list spells differently (MICB_MICA = MICA_MICB)
                                   and accessions a list gets wrong (IFNA1_IFNA13
                                   ships the retired P01562)
      2. (tag, assay name)      -- the assay's own name as the panel list writes it,
                                   punctuation ignored ("H1_1" = "H1-1"); only
                                   well-formed accessions, and only names that
                                   map to one accession within the tag
      3. (tag, gene)            -- the assay's representative HGNC gene
      4. (tag, protein_assay)   -- HGNC normalization changed the gene name
                                   (e.g. Olink JCHAIN -> HGNC IGJ)
      5. (tag, normalized gene) -- panel display name differs only in
                                   punctuation/case (NULISA "BMP-7" vs "BMP7")
      6. rename.csv by gene     -- e.g. FIGF -> VEGFD, MAPT -> tTau/pTau*
    The assay name goes first because a complex's representative gene can also
    name a single-subunit assay on the same panel: Explore HT has both
    IL12A_IL12B and IL12B, so a (tag, gene) join gave the complex IL12B's ID.
    Multiple UniProts are collapsed by first occurrence to keep one row per assay.
    """
    if not ENSG_IN.exists():
        print(f"[skip] {ENSG_IN} not found -- run 03_annot_ensg.R first")
        return
    ensg = pd.read_csv(ENSG_IN)
    merged = ensg.assign(uniprot_id=pd.NA,
                         assay_norm=ensg["protein_assay"].str.split("__").str[0].map(norm_name))

    valid = out_for_merge[out_for_merge["uniprot_id"].astype(str).str.fullmatch(UNIPROT_RE)]
    by_name = valid.assign(assay_norm=valid["gene"].map(norm_name))[
        ["tag", "assay_norm", "uniprot_id"]].drop_duplicates()
    ambiguous = by_name.duplicated(subset=["tag", "assay_norm"], keep=False)
    by_name = by_name[~ambiguous].rename(columns={"uniprot_id": "uniprot_id_by_name"})

    def fill_from(frame, keys, source_col):
        """Fill still-missing uniprot_id from `frame`, joined on `keys`."""
        nonlocal merged
        merged = merged.merge(frame, on=keys, how="left")
        gap = merged["uniprot_id"].isna()
        merged.loc[gap, "uniprot_id"] = merged.loc[gap, source_col]
        merged = merged.drop(columns=[source_col])
        return int(gap.sum() - merged["uniprot_id"].isna().sum())

    rename_map = load_rename_map()
    gap = merged["uniprot_id"].isna()
    merged.loc[gap, "uniprot_id"] = merged.loc[gap, "assay_norm"].map(
        lambda a: rename_map[a][1] if a in rename_map else None
    )
    n_by_rename_assay = int(gap.sum() - merged["uniprot_id"].isna().sum())

    n_by_name = fill_from(by_name, ["tag", "assay_norm"], "uniprot_id_by_name")
    merged = merged.drop(columns=["assay_norm"])

    up_lookup = out_for_merge.drop_duplicates(subset=["tag", "gene"], keep="first")[
        ["tag", "gene", "uniprot_id"]
    ].rename(columns={"uniprot_id": "uniprot_id_by_gene"})
    n_by_gene = fill_from(up_lookup, ["tag", "gene"], "uniprot_id_by_gene")

    by_assay = out_for_merge.drop_duplicates(
        subset=["tag", "protein_assay"], keep="first"
    )[["tag", "protein_assay", "uniprot_id"]].rename(
        columns={"uniprot_id": "uniprot_id_by_assay"}
    )
    n_by_assay = fill_from(by_assay, ["tag", "protein_assay"], "uniprot_id_by_assay")

    norm_lookup = out_for_merge.assign(gene_norm=out_for_merge["gene"].map(norm_name))
    norm_lookup = norm_lookup.drop_duplicates(subset=["tag", "gene_norm"], keep="first")[
        ["tag", "gene_norm", "uniprot_id"]
    ].rename(columns={"uniprot_id": "uniprot_id_by_norm"})
    merged["gene_norm"] = merged["gene"].map(norm_name)
    n_by_norm = fill_from(norm_lookup, ["tag", "gene_norm"], "uniprot_id_by_norm")

    gap = merged["uniprot_id"].isna()
    merged.loc[gap, "uniprot_id"] = merged.loc[gap, "gene_norm"].map(
        lambda g: rename_map[g][1] if g in rename_map else None
    )
    n_by_rename = int(gap.sum() - merged["uniprot_id"].isna().sum())
    merged = merged.drop(columns=["gene_norm"])
    print(f"  uniprot filled by: assay_name={n_by_name} "
          f"(ambiguous names skipped: {int(ambiguous.sum())}), "
          f"rename.csv by assay={n_by_rename_assay}, gene={n_by_gene}, "
          f"protein_assay={n_by_assay}, normalized_gene={n_by_norm}, "
          f"rename.csv by gene={n_by_rename}")
    if len(merged) != len(ensg):
        raise SystemExit(f"[FAIL] merge changed the row count: {len(ensg)} -> {len(merged)}")
    n_total = len(merged)
    n_with_up = merged["uniprot_id"].notna().sum()
    merged.to_csv(ANNOT_OUT, index=False)
    print(f"\nWrote {n_total:,} rows -> {ANNOT_OUT}")
    print(f"  with uniprot_id: {n_with_up:,} ({100*n_with_up/n_total:.1f}%)")
    miss = merged[merged["uniprot_id"].isna()]
    if len(miss):
        print(f"  missing uniprot by tag (top 10):")
        print(miss.groupby("tag").size().sort_values(ascending=False).head(10).to_string())


if __name__ == "__main__":
    # --merge-only rebuilds assay_annot.csv from the committed assay_uniprot.csv,
    # for machines without the raw panel lists (the cluster).
    if "--merge-only" in sys.argv[1:]:
        lists = pd.read_csv(OUT)
        merge_annot(lists.assign(protein_assay=lists["gene"].astype(str) + "__" + lists["tag"]))
    else:
        main()
