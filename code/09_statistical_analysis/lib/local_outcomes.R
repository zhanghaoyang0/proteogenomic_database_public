# Individual-level BF2 outcomes shared by code/09_statistical_analysis/local/11a, 11c and 12, so the three use
# the same MRI and clinical definitions (Lu et al.'s rs5848 manuscript, BF2 part).
# Mac only; source after lib/local_data.R (PET_X, norm_id, parse_date). Nothing here is
# written to disk.

# MRI of CU participants, one row per scan, with the outcomes: FreeSurfer (longitudinal
# stream) volumes / ICV per hemisphere, regional thickness, bilateral HV and ILV / ICV, and
# the surface-area-weighted thickness composites; years from each participant's first scan.
# Returns list(mri, measures = the outcome table, e4 = APOE e4 carrier by sid).
load_mri_cu <- function() {
  hdr <- names(fread(PET_X, nrows = 0))
  SUB <- c("Lateral_Ventricle", "Inf_Lat_Vent", "Thalamus_Proper", "Caudate", "Putamen",
           "Pallidum", "Hippocampus", "Amygdala", "Accumbens_area", "VentralDC")
  SUB <- c(SUB, "Cerebellum_Cortex", "Cerebellum_White_Matter")
  sub_cols <- paste0("aseg_vol_", c(paste0("Left_", SUB), paste0("Right_", SUB),
                     "3rd_Ventricle", "4th_Ventricle", "Brain_Stem",
                     "CC_Anterior", "CC_Mid_Anterior", "CC_Central", "CC_Mid_Posterior", "CC_Posterior"))
  gv_cols <- grep("^aparc_grayvol_[a-z]+_[LR]$", hdr, value = TRUE)
  ct_cols <- grep("^aparc_ct_avg_[a-z]+_[LR]$", hdr, value = TRUE)
  # Bilateral surface-area-weighted thickness composites (the export's ct_composites block,
  # e.g. ct_adsign_surfwt_lr); Lu et al. map frontal, lateral parietal, ... composites.
  comp_cols <- grep("^ct_.*surfwt.*_lr$", hdr, value = TRUE)
  message("thickness composites: ", if (length(comp_cols)) paste(comp_cols, collapse = ", ") else "none found")
  need <- c("sid", "Study", "excluded", "cognitive_status_baseline_variable", "apoe_genotype_baseline_variable",
            "proc_mri_date", "icv_mm3", sub_cols, gv_cols, ct_cols, comp_cols)
  miss <- setdiff(need, hdr)
  if (length(miss)) stop("PET export is missing: ", paste(head(miss, 10), collapse = ", "))
  if (length(gv_cols) != 68 || length(ct_cols) != 68) stop("expected 68 Desikan regions per measure")

  mri <- fread(PET_X, select = need, colClasses = list(character = c("sid", "proc_mri_date")), showProgress = FALSE)
  mri[, sid := norm_id(sid)]
  mri <- mri[Study == "BF2" & !excluded %in% 1]
  cu_sid <- unique(mri[trimws(cognitive_status_baseline_variable) %in% c("Normal", "SCD"), sid])
  e4 <- unique(mri[, .(sid, code = suppressWarnings(as.integer(apoe_genotype_baseline_variable)))], by = "sid")
  e4[, e4_carrier := as.integer(grepl("4", code))][is.na(code), e4_carrier := NA_integer_]
  mri <- mri[sid %in% cu_sid & !is.na(icv_mm3)]
  # parse_date (lib/local_data.R): the export's date format is not documented.
  n_nonempty <- mri[!is.na(proc_mri_date) & trimws(proc_mri_date) != "", .N]
  mri[, mri_date := parse_date(proc_mri_date, "proc_mri_date")]
  if (mri[!is.na(mri_date), .N] < 0.99 * n_nonempty) stop("proc_mri_date does not parse for > 1% of non-empty values")
  mri <- mri[!is.na(mri_date)]
  # One row per scan: the export repeats a scan on every visit it is closest to.
  setorder(mri, sid, mri_date)
  mri <- unique(mri, by = c("sid", "mri_date"))
  num <- c("icv_mm3", sub_cols, gv_cols, ct_cols, comp_cols)
  mri[, (num) := lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols = num]
  mri[, years := as.numeric(mri_date - min(mri_date)) / 365.25, by = sid]

  # Outcomes: volumes / ICV per hemisphere, thickness as is, HV and ILV bilateral.
  vol <- c(sub_cols, gv_cols)
  mri[, (vol) := lapply(.SD, function(v) v / icv_mm3), .SDcols = vol]
  mri[, `:=`(HV  = (aseg_vol_Left_Hippocampus + aseg_vol_Right_Hippocampus) / 2,
             ILV = (aseg_vol_Left_Inf_Lat_Vent + aseg_vol_Right_Inf_Lat_Vent) / 2)]
  MEASURES <- rbind(
    data.table(col = c("HV", "ILV"), measure = "primary"),
    data.table(col = sub_cols, measure = "subcortical volume"),
    data.table(col = gv_cols, measure = "cortical volume"),
    data.table(col = ct_cols, measure = "cortical thickness"),
    data.table(col = comp_cols, measure = "thickness composite"))
  MEASURES[, region := sub("^(aseg_vol_|aparc_grayvol_|aparc_ct_avg_)", "", col)]
  MEASURES <- MEASURES[col %in% names(mri)]
  list(mri = mri, measures = MEASURES, e4 = e4)
}

# Clinical data: base (one row per BF2 sid at Visit 0, with baseline centiloid, CU and
# AD-focused flags), long (mPACC v3 visits, years from baseline) and cu (CU participants
# with time to AC-CI / AD-CI). Returns list(base, long, cu).
load_clinical <- function() {
  COLS <- c(sid = "sid", Study = "Study", excluded = "excluded", Visit = "Visit", visit_date = "visit_date",
            age = "age", mpacc = "mPACC_v3", cog = "cognitive_status_baseline_variable",
            apoe = "apoe_genotype_baseline_variable", edu = "education_level_years_baseline_variable",
            abeta = "Abnormal_CSF_Ab42_Ab40_Ratio", centiloid = "CL_bfnc_ber_com_composite",
            mci_dich = "converted_MCI_dich_baseline_variable", mci_date = "converted_MCI_date_baseline_variable",
            mci_last = "last_visit_date_MCI_nonconverted_baseline_variable",
            dem_dich = "converted_dementia_dich_baseline_variable", dem_date = "converted_dementia_date_baseline_variable",
            dem_last = "last_visit_date_dementia_nonconverted_baseline_variable")
  hdr <- names(fread(PET_X, nrows = 0))
  miss <- setdiff(COLS, hdr)
  if (length(miss)) stop("PET export is missing: ", paste(miss, collapse = ", "))
  cl <- fread(PET_X, select = unname(COLS), colClasses = "character", showProgress = FALSE)
  setnames(cl, unname(COLS), names(COLS))
  cl[, sid := norm_id(sid)]
  cl <- cl[Study == "BF2" & !excluded %in% "1"]
  num <- c("Visit", "age", "mpacc", "edu", "abeta", "centiloid", "mci_dich", "dem_dich")
  cl[, (num) := lapply(.SD, function(v) suppressWarnings(as.numeric(v))), .SDcols = num]
  for (v in c("visit_date", "mci_date", "mci_last", "dem_date", "dem_last"))
    set(cl, j = v, value = parse_date(cl[[v]], v))
  if (anyDuplicated(cl[, .(sid, Visit)])) stop("several rows per sid x Visit")

  base <- cl[Visit == 0, .(sid, base_date = visit_date, age, edu, abeta, cog = trimws(cog),
                           e4_carrier = as.integer(grepl("4", apoe)),
                           mci_dich, mci_date, mci_last, dem_dich, dem_date, dem_last)]
  base[is.na(suppressWarnings(as.integer(cl[Visit == 0, apoe]))), e4_carrier := NA_integer_]
  # Baseline centiloid: the first available Abeta-PET.
  cen0 <- cl[!is.na(centiloid)][order(Visit)][, .SD[1], by = sid][, .(sid, centiloid)]
  base <- merge(base, cen0, by = "sid", all.x = TRUE)
  base[, cu := cog %in% c("Normal", "SCD")]
  base[, ad_focused := cu | (cog %in% c("MCI", "Dementia") & abeta %in% 1)]

  long <- merge(cl[!is.na(mpacc) & !is.na(visit_date), .(sid, visit_date, mpacc)],
                base[, .(sid, base_date)], by = "sid")
  long[, years := as.numeric(visit_date - base_date) / 365.25]
  long <- long[years >= -0.25]   # a test dated just before the baseline visit counts as baseline
  long[years < 0, years := 0]

  # Progression in CU: first MCI or dementia; censored at the last non-converted visit.
  cu <- base[cu == TRUE]
  cu[, event_date := pmin(fifelse(mci_dich %in% 1, mci_date, as.IDate(NA)),
                          fifelse(dem_dich %in% 1, dem_date, as.IDate(NA)), na.rm = TRUE)]
  # A converter without a conversion date has no event time: excluded, not censored.
  no_date <- cu[(mci_dich %in% 1 | dem_dich %in% 1) & is.na(event_date), .N]
  cu <- cu[!((mci_dich %in% 1 | dem_dich %in% 1) & is.na(event_date))]
  cu[, ac_ci := as.integer(!is.na(event_date))]
  cu[, end_date := fifelse(ac_ci == 1, event_date, pmax(mci_last, dem_last, na.rm = TRUE))]
  cu[, t_years := as.numeric(end_date - base_date) / 365.25]
  n_bad <- cu[is.na(t_years) | t_years <= 0, .N]
  message(sprintf("CU at baseline: %d; AC-CI events %d; dropped: converted without a date %d, without positive follow-up %d",
                  nrow(cu), cu[ac_ci == 1, .N], no_date, n_bad))
  cu <- cu[!is.na(t_years) & t_years > 0]
  # AD-CI: an event only if Abeta+ at baseline; Abeta- converters are censored at conversion.
  # Without baseline Abeta the AD-CI status is unknown: NA, so the AD-CI models exclude the
  # participant (AC-CI keeps them).
  cu[, ad_ci := fifelse(is.na(abeta), NA_integer_, as.integer(ac_ci == 1 & abeta == 1))]
  message(sprintf("AD-CI: %d CU without baseline Abeta excluded; %d events", cu[is.na(ad_ci), .N],
                  cu[ad_ci == 1, .N]))
  list(base = base, long = long, cu = cu)
}
