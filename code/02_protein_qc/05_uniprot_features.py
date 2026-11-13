#!/usr/bin/env python3
"""
Parse the local SwissProt human proteome flat file and extract per-UniProt features
for the accessions present in info/bf_protein/assay_annot.csv.

Outputs:
    info/bf_protein/uniprot/features.csv     one row per uniprot_id (wide)
    info/bf_protein/uniprot/xrefs.csv.gz     one row per cross-reference (long)

Source is the UniProt *reference proteome* (UP000005640_9606.dat), which keeps
one representative per gene. That is the right table for annotating assay
targets, but it means a handful of real accessions are absent because a sibling
gene product represents the gene (e.g. P58400, Neurexin-1-beta). Those are
reported as missing rather than silently filled.

features.csv (per uniprot_id):
    matched_accession       the accession whose record supplied these values
    is_isoform              1 when uniprot_id carries an isoform suffix and the
                            values therefore come from the canonical sequence
    length                  sequence length (aa)
    mass_da                 molecular weight (Da), from the record's SQ line
    signal_peptide          1 if a SIGNAL feature exists, else 0
    n_transmembrane         count of TRANSMEM features
    n_domain                count of DOMAIN features
    n_disulfide             count of DISULFID features
    n_glycosylation         count of CARBOHYD features
    n_phosphorylation       count of MOD_RES with "Phospho" in note
    n_modres                count of all MOD_RES features (any PTM)
    subcellular_location    every SUBCELLULAR LOCATION comment, evidence codes
                            and Note= tails removed, blocks joined by " | "
    is_secreted             1 if "Secreted" appears in subcellular location or keywords
    is_membrane             1 if "Membrane" appears in subcellular location or keywords
    n_go_bp/n_go_mf/n_go_cc GO term counts per aspect
    pfam / n_pfam           Pfam accessions (semicolon-joined) and count
    interpro / n_interpro   InterPro accessions (semicolon-joined) and count
    n_reactome              Reactome pathway count
    n_drugbank              DrugBank entries -- nonzero means a known drug target
    keywords                semicolon-joined keyword list

xrefs.csv.gz (per uniprot_id x cross-reference), for GO, Pfam, InterPro,
Reactome and DrugBank: uniprot_id, db, xref_id, name, detail. The wide table
carries counts; this is where the identities live.
"""

from __future__ import annotations
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SwissProt

REPO = Path(__file__).resolve().parents[2]
ASSAY_ANNOT = REPO / "info/bf_protein/assay_annot.csv"
OUT_CSV     = REPO / "info/bf_protein/uniprot/features.csv"
OUT_XREF    = REPO / "info/bf_protein/uniprot/xrefs.csv.gz"
DAT_FILE    = Path.home() / "db/shared_data/uniport/UP000005640_9606.dat"

XREF_DBS = ("GO", "Pfam", "InterPro", "Reactome", "DrugBank")

EVIDENCE_RE = re.compile(r"\{[^{}]*\}")
NOTE_RE = re.compile(r"\bNote=.*$", re.DOTALL)


def clean_location(comment: str) -> str:
    """One SUBCELLULAR LOCATION comment -> its locations, evidence removed.

    The previous version kept only the first comment and cut it at the first
    evidence brace, which truncated 3,419 of 6,249 proteins and dropped the
    second location of every multi-compartment protein. Here the whole comment
    is kept; only the ECO evidence tags and the free-text Note= tail go.
    """
    text = comment[len("SUBCELLULAR LOCATION:"):]
    text = EVIDENCE_RE.sub("", text)
    text = NOTE_RE.sub("", text)
    text = re.sub(r"\s+", " ", text)
    return text.strip(" .;")


def extract_record(rec) -> dict:
    feats = rec.features
    types = [ft.type for ft in feats]

    n_phospho = sum(
        1 for ft in feats
        if ft.type == "MOD_RES" and "phospho" in str(ft.qualifiers.get("note", "")).lower()
    )

    subloc_comments = [c for c in rec.comments if c.startswith("SUBCELLULAR LOCATION:")]
    locations = [clean_location(c) for c in subloc_comments]
    subloc = " | ".join(x for x in locations if x)

    # The stored field drops the free-text Note= tail, but the flags are read
    # off the full comment: UniProt sometimes states a compartment only inside
    # the Note (FGF1 P05230, PLG P00747). Judging the flags from the trimmed
    # field alone loses 11 of these.
    kw_blob = "; ".join(rec.keywords)
    full_text = " ".join(subloc_comments).lower()
    is_secreted = int("secreted" in full_text or "Secreted" in rec.keywords)
    is_membrane = int("membrane" in full_text or "Membrane" in rec.keywords)

    # Cross-references. Tuple shape is (db, id, *rest); rest is db-specific, so
    # index defensively.
    xrefs, go_aspect, per_db = [], Counter(), Counter()
    for x in rec.cross_references:
        db = x[0]
        if db not in XREF_DBS:
            continue
        xref_id = x[1] if len(x) > 1 else ""
        third = x[2] if len(x) > 2 else ""
        fourth = x[3] if len(x) > 3 else ""
        per_db[db] += 1
        if db == "GO":
            # third is "C:cell surface" / "P:..." / "F:..."; fourth is evidence.
            aspect, _, term = str(third).partition(":")
            go_aspect[aspect] += 1
            xrefs.append((db, xref_id, term, f"{aspect}|{fourth}"))
        else:
            # Pfam carries a copy count in the fourth field; the others do not.
            xrefs.append((db, xref_id, str(third), str(fourth)))

    pfam = [x[1] for x in xrefs if x[0] == "Pfam"]
    interpro = [x[1] for x in xrefs if x[0] == "InterPro"]

    # Mass from the record itself (the SQ line) rather than recomputed from the
    # sequence: it is what UniProt states, and it cannot fail on a sequence
    # containing an ambiguous residue.
    try:
        mass_da = float(rec.seqinfo[1])
    except (TypeError, IndexError, ValueError):
        mass_da = None

    row = {
        "length":               rec.sequence_length,
        "mass_da":              mass_da,
        "signal_peptide":       int("SIGNAL" in types),
        "n_transmembrane":      types.count("TRANSMEM"),
        "n_domain":             types.count("DOMAIN"),
        "n_disulfide":          types.count("DISULFID"),
        "n_glycosylation":      types.count("CARBOHYD"),
        "n_phosphorylation":    n_phospho,
        "n_modres":             types.count("MOD_RES"),
        "subcellular_location": subloc,
        "is_secreted":          is_secreted,
        "is_membrane":          is_membrane,
        "n_go_bp":              go_aspect.get("P", 0),
        "n_go_mf":              go_aspect.get("F", 0),
        "n_go_cc":              go_aspect.get("C", 0),
        "pfam":                 ";".join(pfam),
        "n_pfam":               len(pfam),
        "interpro":             ";".join(interpro),
        "n_interpro":           len(interpro),
        "n_reactome":           per_db.get("Reactome", 0),
        "n_drugbank":           per_db.get("DrugBank", 0),
        "keywords":             kw_blob,
    }
    return row, xrefs


def main() -> int:
    for p in (ASSAY_ANNOT, DAT_FILE):
        if not p.exists():
            sys.exit(f"ERROR: missing input {p}")

    annot = pd.read_csv(ASSAY_ANNOT, low_memory=False)
    wanted = set(annot["uniprot_id"].dropna().astype(str).str.strip())
    wanted.discard("")
    # The reference proteome holds canonical sequences only, so an isoform id
    # ("P12345-2") can only be served by its canonical record. Do that, but say
    # so in the output instead of letting the substitution pass unnoticed.
    base_to_orig: dict[str, set[str]] = {}
    for acc in wanted:
        base_to_orig.setdefault(acc.split("-", 1)[0], set()).add(acc)
    print(f"unique uniprot_id to look up: {len(wanted):,} "
          f"({len(base_to_orig):,} base accessions)", flush=True)

    found: dict[str, dict] = {}
    xref_rows: list[tuple] = []
    n_scanned = 0
    with open(DAT_FILE) as fh:
        for rec in SwissProt.parse(fh):
            n_scanned += 1
            hit_bases = set(base_to_orig).intersection(rec.accessions)
            if not hit_bases:
                continue
            row, xrefs = extract_record(rec)
            for base in hit_bases:
                for acc in base_to_orig[base]:
                    found[acc] = {
                        "matched_accession": base,
                        "is_isoform": int(acc != base),
                        **row,
                    }
                    xref_rows += [(acc, *x) for x in xrefs]
            if len(found) == len(wanted):
                break
    print(f"scanned {n_scanned:,} SwissProt records; matched {len(found):,}/{len(wanted):,}",
          flush=True)

    missing = sorted(wanted - found.keys())
    if missing:
        composite = [m for m in missing if "_" in m]
        single = [m for m in missing if "_" not in m]
        print(f"WARN: {len(missing)} uniprot_id not in the reference proteome.", flush=True)
        print(f"  {len(composite)} composite accession string(s) -- an assay targeting a "
              f"complex, which has no single record: {composite[:5]}", flush=True)
        print(f"  {len(single)} single accession(s) absent because a sibling gene product "
              f"represents the gene: {single[:10]}", flush=True)

    n_iso = sum(r["is_isoform"] for r in found.values())
    if n_iso:
        print(f"NOTE: {n_iso} isoform id(s) carry their canonical sequence's features "
              f"(is_isoform=1).", flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame.from_dict(found, orient="index").rename_axis("uniprot_id").reset_index()
    out = out.sort_values("uniprot_id").reset_index(drop=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"wrote {OUT_CSV} ({len(out):,} rows, {len(out.columns)} cols)")

    xr = pd.DataFrame(xref_rows, columns=["uniprot_id", "db", "xref_id", "name", "detail"])
    xr = xr.sort_values(["uniprot_id", "db", "xref_id"]).reset_index(drop=True)
    xr.to_csv(OUT_XREF, index=False, compression="gzip")
    print(f"wrote {OUT_XREF} ({len(xr):,} rows)")
    print("  per db:", dict(xr["db"].value_counts()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
