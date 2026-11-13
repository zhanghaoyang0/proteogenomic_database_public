#!/usr/bin/env Rscript

# Owns the FinnGen raw-category -> short disease-class taxonomy.
#
# Reads the FinnGen R12 manifests, enumerates every `category` value actually
# present, validates that the curated map below covers all of them (fail loud on
# drift -- a new FinnGen release must NOT silently fall through to "Drug/other"),
# then writes the lookup table consumed downstream:
#
#   OUT: info/finngen_category_map.csv   (finngen_category, class)
#        read by code/06_external_data/collect_trait_info.R to set trait_info category.
#
# Curation lives here (single source of truth). To change the taxonomy, edit the
# `m` vector and re-run this script. Semantic remaps fold small FinnGen
# endpoint-groups into their parent ICD chapter:
#   Diabetes endpoints            -> Endocrine & metabolic
#   Gastrointestinal endpoints    -> Digestive
#   Interstitial lung disease ep. -> Respiratory
# There is no count-based folding anywhere: a class keeps its name however few traits
# it has (the "< 5 traits -> Drug/other" rule in collect_trait_info.R was removed
# on 2026-10-22). Catch-all groups get explicit names (Medication, Comorbidity,
# Health status/contact, Other/unclassified) instead of one "Drug/other".
#
# Usage:  Rscript code/06_external_data/finngen_reclass.R

library(data.table)

args_all   <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- dirname(normalizePath(script_file))
repo_root  <- normalizePath(file.path(script_dir, "..", ".."))
setwd(repo_root)

FIN_INFO  <- path.expand("~/db/shared_data/gwas/finngen/info")
MANIFESTS <- file.path(FIN_INFO, c(
  "finngen_R12_manifest.tsv",                                  # category source used by collector
  "FinnGen_R12_MVP_UKBB_manifest.tsv",
  "finngen_R12_meta_analysis_mapping_with_definitions.tsv"
))
OUT_CSV <- "info/finngen_category_map.csv"

# ---- curated taxonomy: raw FinnGen `category` -> short disease class ----
m <- c(
  "I Certain infectious and parasitic diseases (AB1_)"                                                            = "Infectious",
  "II Neoplasms, from cancer register (ICD-O-3)"                                                                  = "Neoplasms",
  "II Neoplasms from hospital discharges (CD2_)"                                                                  = "Neoplasms",
  "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)" = "Blood & immune",
  "IV Endocrine, nutritional and metabolic diseases (E4_)"                                                        = "Endocrine & metabolic",
  "V Mental and behavioural disorders (F5_)"                                                                      = "Psychiatric",
  "VI Diseases of the nervous system (G6_)"                                                                       = "Neurological",
  "VII Diseases of the eye and adnexa (H7_)"                                                                      = "Eye",
  "VIII Diseases of the ear and mastoid process (H8_)"                                                            = "Ear",
  "IX Diseases of the circulatory system (I9_)"                                                                   = "Circulatory",
  "X Diseases of the respiratory system (J10_)"                                                                   = "Respiratory",
  "XI Diseases of the digestive system (K11_)"                                                                    = "Digestive",
  "XII Diseases of the skin and subcutaneous tissue (L12_)"                                                       = "Skin",
  "XIII Diseases of the musculoskeletal system and connective tissue (M13_)"                                      = "Musculoskeletal",
  "XIV Diseases of the genitourinary system (N14_)"                                                               = "Genitourinary",
  "XV Pregnancy, childbirth and the puerperium (O15_)"                                                            = "Pregnancy",
  "XVI Certain conditions originating in the perinatal period (P16_)"                                             = "Perinatal",
  "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)"                              = "Congenital",
  "XVIII Symptoms, signs and abnormal clinical and laboratory findings, not elsewhere classified (R18_)"          = "Symptoms",
  "XIX Injury, poisoning and certain other consequences of external causes (ST19_)"                              = "Injury",
  "XXI Factors influencing health status and contact with health services (Z21_)"                                = "Health status/contact",
  "XXII Codes for special purposes (U22_)"                                                                        = "Other/unclassified",
  "Alcohol related diseases"                                                                                      = "Alcohol-related",
  "Asthma and related endpoints"                                                                                  = "Asthma",
  "Cardiometabolic endpoints"                                                                                     = "Cardiometabolic",
  "Comorbidities of Asthma"                                                                                       = "Comorbidity",
  "Comorbidities of COPD"                                                                                         = "Comorbidity",
  "Comorbidities of Diabetes"                                                                                     = "Comorbidity",
  "Comorbidities of Gastrointestinal endpoints"                                                                   = "Comorbidity",
  "Comorbidities of Interstitial lung disease endpoints"                                                         = "Comorbidity",
  "Comorbidities of Neurological endpoints"                                                                       = "Comorbidity",
  "COPD and related endpoints"                                                                                    = "COPD",
  "Dental endpoints"                                                                                              = "Digestive",
  "Diabetes endpoints"                                                                                            = "Endocrine & metabolic",
  "Diseases marked as autimmune origin"                                                                           = "Autoimmune",
  "Drug purchase endpoints"                                                                                       = "Medication",
  "Gastrointestinal endpoints"                                                                                    = "Digestive",
  "Interstitial lung disease endpoints"                                                                           = "Respiratory",
  "Miscellaneous, not yet classified endpoints"                                                                   = "Other/unclassified",
  "Neurological endpoints"                                                                                        = "Neurological",
  "Other, not yet classified endpoints (same as #MISC)"                                                          = "Other/unclassified",
  "Psychiatric endpoints from Katri Räikkönen"                                                          = "Psychiatric",
  "Quantitative endpoints"                                                                                        = "Quantitative",
  "Rheuma endpoints"                                                                                              = "Rheumatological"
)
cmap <- data.table(finngen_category = names(m), class = unname(m))
if (anyDuplicated(cmap$finngen_category)) stop("raw category mapped twice in `m`")
if (anyNA(cmap$class) || any(!nzchar(cmap$class))) stop("empty class in `m`")

# ---- enumerate raw categories actually present across the manifests ----
present <- unique(unlist(lapply(MANIFESTS, function(p) {
  # Stop, not skip: the drift check below is only as good as the manifests it reads.
  if (!file.exists(p)) stop("missing FinnGen manifest: ", p)
  d <- fread(p, select = "category", quote = "", showProgress = FALSE)
  trimws(as.character(d$category))
})))
present <- present[!is.na(present) & nzchar(present)]
message("raw categories present across manifests: ", length(present))

# ---- fail loud on drift: any present category missing from the curated map ----
missing <- setdiff(present, cmap$finngen_category)
if (length(missing)) {
  message("UNCOVERED FinnGen categories (present in manifest, absent from curated map):")
  for (mc in missing) message("  - ", mc)
  stop(length(missing), " FinnGen category(ies) not classified -- add them to `m` in ",
       "code/06_external_data/finngen_reclass.R (do not let them silently fall to Drug/other)")
}

# ---- warn (non-fatal) on stale curated entries not in any current manifest ----
stale <- setdiff(cmap$finngen_category, present)
if (length(stale)) message("note: ", length(stale), " curated categories not seen in current manifests (kept): ",
                           paste(stale, collapse = " | "))

setorder(cmap, finngen_category)
# Write to a temp file and rename, so a failed run never leaves a partial map.
tmp_csv <- paste0(OUT_CSV, ".tmp")
fwrite(cmap, tmp_csv)
if (!file.rename(tmp_csv, OUT_CSV)) stop("could not move ", tmp_csv, " to ", OUT_CSV)
message("wrote ", nrow(cmap), " mappings -> ", OUT_CSV,
        "  (", uniqueN(cmap$class), " short classes)")
