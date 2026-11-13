#!/usr/bin/env bash
# Genetic ancestry of the BF1/BF2 GWAS samples, by projection onto 1000 Genomes.
#
# Reference: 1000 Genomes phase 3 (hg19, ~/db/shared_data/ref_geno/1kg_v3), the five
# superpopulation files AFR AMR EAS EUR SAS (2,504 samples). SNPs: common in EUR
# (MAF >= 0.05), biallelic ACGT, not A/T or C/G, autosomal, outside the MHC and the
# high-LD ranges, LD-pruned in EUR (200 / 50 / r2 0.1), and present in all five
# reference files and in both BF GWAS pfiles. PCA (10 PCs) is fitted on the reference;
# reference and BF samples are then projected with the same allele weights, so both
# sit on one scale. Each BF sample takes the majority superpopulation of its 20 nearest
# reference samples on PC1-6; it is also flagged when it lies more than 6 SD from the
# EUR reference mean on any of PC1-4.
#
# Output (secure): ~/db/secure_data/bf/array/ancestry/
#   bf{1,2}_ancestry.tsv   one row per GWAS sample: PC1-10, nearest superpopulation,
#                          share of EUR neighbours, outside_eur_6sd
#   ancestry_pca.png       PC1 vs PC2, reference coloured, BF samples on top
# The screen gets counts only.
# Stage 7 of 02_postimpute_qc: run after Stage 2, whenever the GWAS pfiles change.
#   bash 07_ancestry_pca.sh
set -euo pipefail

BF=${BF_ROOT:-$HOME/db/secure_data/bf}
REF=${REF:-$HOME/db/shared_data/ref_geno/1kg_v3}
PY=${PY:-python3}
PLINK=${PLINK:-$(command -v plink || echo /opt/homebrew/bin/plink)}
PLINK2=${PLINK2:-$(command -v plink2 || echo /opt/homebrew/bin/plink2)}
OUT=$BF/array/ancestry
PF_BF1=$BF/array/bf1/imputed_qc/bf1_hg19
PF_BF2=$BF/array/bf2/imputed_qc/bf2_hg19
POPS="AFR AMR EAS EUR SAS"
for x in "$PLINK" "$PLINK2" "$PY"; do [[ -x "$x" ]] || { echo "ERROR: not executable: $x" >&2; exit 1; }; done
for p in $POPS; do [[ -s "$REF/$p.bed" ]] || { echo "ERROR: missing $REF/$p.bed" >&2; exit 1; }; done
for f in "$PF_BF1.pvar" "$PF_BF2.pvar" "$REF/highld_hg19.range" "$REF/mhc_hg19.range"; do
  [[ -s "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done
w=$(mktemp -d)
trap 'rm -rf "$w"' EXIT
mkdir -p "$OUT"
run() {  # run <log name> <command...>: quiet, but show plink's error and stop on failure
  local n=$1; shift
  "$@" > "$w/$n.out" 2>&1 || { grep -A3 -i '^error' "$w/$n.out" | cut -c1-200 >&2; tail -3 "$w/$n.out" >&2; exit 1; }
}

# 1. Common, unambiguous, pruned SNPs in EUR
cat "$REF/highld_hg19.range" "$REF/mhc_hg19.range" | awk 'NF >= 3 {print $1 "\t" $2 "\t" $3 "\t" (NF >= 4 ? $4 : "x")}' \
  | sed 's/^chr//' > "$w/exclude.range"
run prune "$PLINK2" --bfile "$REF/EUR" --autosome --snps-only just-acgt --max-alleles 2 --maf 0.05 \
  --exclude range "$w/exclude.range" --indep-pairwise 200 50 0.1 --out "$w/eur"
awk '($5 $6) !~ /^(AT|TA|CG|GC)$/ {print $2}' "$REF/EUR.bim" > "$w/nonpal.txt"
awk 'NR == FNR {k[$1]; next} $1 in k' "$w/nonpal.txt" "$w/eur.prune.in" > "$w/cand.txt"
echo "EUR common pruned unambiguous SNPs: $(wc -l < "$w/cand.txt")"

# 2. Keep SNPs present everywhere
cp "$w/cand.txt" "$w/keep.txt"
for p in AFR AMR EAS SAS; do
  awk 'NR == FNR {k[$1]; next} $2 in k {print $2}' "$w/keep.txt" "$REF/$p.bim" > "$w/k.tmp"; mv "$w/k.tmp" "$w/keep.txt"
done
for pv in "$PF_BF1.pvar" "$PF_BF2.pvar"; do
  awk 'NR == FNR {k[$1]; next} !/^#/ && ($3 in k) {print $3}' "$w/keep.txt" "$pv" > "$w/k.tmp"; mv "$w/k.tmp" "$w/keep.txt"
done
n_keep=$(wc -l < "$w/keep.txt")
echo "SNPs in all five reference files and both BF pfiles: $n_keep"
[[ "$n_keep" -ge 20000 ]] || { echo "ERROR: too few shared SNPs ($n_keep); check variant IDs" >&2; exit 1; }

# 3. Merge the reference (plink 1.9 merge; drop SNPs it cannot reconcile and retry once)
for p in $POPS; do run "ext_$p" "$PLINK" --bfile "$REF/$p" --extract "$w/keep.txt" --make-bed --out "$w/$p"; done
printf '%s\n' "$w/AMR" "$w/EAS" "$w/EUR" "$w/SAS" > "$w/merge.txt"
if ! "$PLINK" --bfile "$w/AFR" --merge-list "$w/merge.txt" --make-bed --out "$w/ref" > "$w/merge.out" 2>&1; then
  [[ -s "$w/ref-merge.missnp" ]] || { tail -5 "$w/merge.out" >&2; exit 1; }
  echo "reference merge: dropping $(wc -l < "$w/ref-merge.missnp") SNPs with allele conflicts"
  for p in $POPS; do run "ex2_$p" "$PLINK" --bfile "$w/$p" --exclude "$w/ref-merge.missnp" --make-bed --out "$w/${p}2"; done
  printf '%s\n' "$w/AMR2" "$w/EAS2" "$w/EUR2" "$w/SAS2" > "$w/merge.txt"
  run merge2 "$PLINK" --bfile "$w/AFR2" --merge-list "$w/merge.txt" --make-bed --out "$w/ref"
fi
for p in $POPS; do awk -v p="$p" '{print $2 "\t" p}' "$REF/$p.fam"; done > "$w/ref_pop.tsv"

# 4. PCA on the reference, then project reference and BF samples with the same weights
run pca "$PLINK2" --bfile "$w/ref" --freq counts --pca 10 allele-wts vcols=chrom,ref,alt --out "$w/refpca"
score() {  # score <name> <plink2 input args...>
  local n=$1; shift
  run "score_$n" "$PLINK2" "$@" --extract "$w/keep.txt" --read-freq "$w/refpca.acount" \
    --score "$w/refpca.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
    --score-col-nums 6-15 --out "$w/proj_$n"
}
score ref --bfile "$w/ref"
score bf1 --pfile "$PF_BF1"
score bf2 --pfile "$PF_BF2"

# 5. Nearest-neighbour superpopulation and EUR outliers; counts to the screen only
"$PY" - "$w" "$OUT" <<'PY'
import sys, numpy as np, pandas as pd
w, out = sys.argv[1:3]
def load(n):
    d = pd.read_csv(f"{w}/proj_{n}.sscore", sep="\t")
    d = d.rename(columns={"#IID": "IID"}) if "#IID" in d.columns else d.rename(columns={"#FID": "FID"})
    pcs = [c for c in d.columns if c.endswith("_AVG")]
    assert len(pcs) == 10, pcs
    d = d[["IID"] + pcs]; d.columns = ["IID"] + [f"PC{i}" for i in range(1, 11)]
    return d
ref = load("ref").merge(pd.read_csv(f"{w}/ref_pop.tsv", sep="\t", names=["IID", "pop"]), on="IID")
assert ref["pop"].notna().all() and len(ref) > 2000, len(ref)
X = ref[[f"PC{i}" for i in range(1, 7)]].to_numpy()
eur = ref[ref["pop"] == "EUR"]
mu, sd = eur[[f"PC{i}" for i in range(1, 5)]].mean(), eur[[f"PC{i}" for i in range(1, 5)]].std()
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(5, 4.5))
for p, g in ref.groupby("pop"):
    ax.scatter(g.PC1, g.PC2, s=4, alpha=.5, label=f"1000G {p}")
for n in ["bf1", "bf2"]:
    d = load(n)
    dist = ((d[[f"PC{i}" for i in range(1, 7)]].to_numpy()[:, None, :] - X[None, :, :]) ** 2).sum(2)
    nn = np.argsort(dist, axis=1)[:, :20]
    labs = ref["pop"].to_numpy()[nn]
    d["nearest_pop"] = pd.DataFrame(labs).mode(axis=1)[0].to_numpy()
    d["eur_share"] = (labs == "EUR").mean(1)
    z = (d[[f"PC{i}" for i in range(1, 5)]] - mu) / sd
    d["outside_eur_6sd"] = (z.abs() > 6).any(axis=1)
    d.to_csv(f"{out}/{n}_ancestry.tsv", sep="\t", index=False)
    print(f"\n{n.upper()}: {len(d)} GWAS samples")
    print("  nearest 1000G superpopulation:", d["nearest_pop"].value_counts().to_dict())
    print("  EUR share of 20 neighbours < 0.8:", int((d["eur_share"] < 0.8).sum()))
    print("  outside EUR mean +- 6 SD on PC1-4:", int(d["outside_eur_6sd"].sum()))
    ax.scatter(d.PC1, d.PC2, s=3, c="black", alpha=.4, label=n.upper())
ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.legend(fontsize=6, markerscale=2)
fig.tight_layout(); fig.savefig(f"{out}/ancestry_pca.png", dpi=200)
print("\nancestry_pca done:", out)
PY
