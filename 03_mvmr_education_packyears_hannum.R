# ============================================================
# MVMR (Multivariable Mendelian Randomization)
# EDUCATION + PACK-YEARS + HANNUM  ->  ORAL CANCER
# X1: Education (ebi-a-GCST90029013)
# X2: Pack-years (ukb-b-10831)
# X3: Hannum Age (ebi-a-GCST90014301)
# Y : Oral cavity cancer (ieu-b-4961)
# AUTHOR: IRMA VALDEBENITO
# ============================================================

# --------------------------
# STEP 0 — PACKAGES
# --------------------------
pkgs <- c("TwoSampleMR","ieugwasr","dplyr","readr","ggplot2","openxlsx","beepr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install)
lapply(pkgs, library, character.only = TRUE)

# --------------------------
# STEP 1 — IDS + OUTPUT FOLDER
# --------------------------
X1_id <- "ebi-a-GCST90029013"  # Education
X2_id <- "ukb-b-10831"         # Pack-years
X3_id <- "ebi-a-GCST90014301"  # Hannum
Y_id  <- "ieu-b-4961"          # Oral cavity cancer

exposure_ids <- c(X1_id, X2_id, X3_id)

base_dir <- "YOUR_FOLDER"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(base_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(base_dir, paste0("runlog_", run_tag, ".txt"))
sink(log_file, split = TRUE)

cat("====================================\n")
cat("MVMR RUN TAG:", run_tag, "\n")
cat("Base dir:", base_dir, "\n")
cat("Exposures:", paste(exposure_ids, collapse = " | "), "\n")
cat("Outcome:", Y_id, "\n")
cat("====================================\n\n")

# --------------------------
# STEP 2 — EXTRACT + COMBINED CLUMPING (KEY STEP)
# --------------------------
# mv_extract_exposures does: extract instruments for each exposure, then combines + clumps
# You can tune pval_threshold / clump_r2 / clump_kb
cat(">> Extracting instruments + combined clumping...\n")
exp_mv <- TwoSampleMR::mv_extract_exposures(
  exposure_ids,
  pval_threshold = 5e-8,
  clump_r2 = 0.001,
  clump_kb = 10000
)

cat("Instruments after combined clumping:\n")
print(table(exp_mv$exposure))
cat("\nUnique SNPs:", length(unique(exp_mv$SNP)), "\n\n")

# Save exposure instruments
readr::write_csv(exp_mv, file.path(base_dir, "data", paste0("mv_exposures_", run_tag, ".csv")))

# --------------------------
# STEP 3 — EXTRACT OUTCOME ASSOCIATIONS (Y) FOR THE SAME SNPS
# --------------------------
cat(">> Extracting outcome associations...\n")
out_dat <- TwoSampleMR::extract_outcome_data(
  snps = unique(exp_mv$SNP),
  outcomes = Y_id
)

cat("Outcome rows:", nrow(out_dat), "\n\n")
readr::write_csv(out_dat, file.path(base_dir, "data", paste0("outcome_", run_tag, ".csv")))

# --------------------------
# STEP 4 — HARMONISE (MVMR FORMAT)  ✅ FIX
# --------------------------
cat(">> Harmonising multivariable data...\n")
mvdat <- TwoSampleMR::mv_harmonise_data(exposure_dat = exp_mv, outcome_dat = out_dat)

# Sanity checks (mvdat is a LIST, not a data.frame)
cat("Harmonised SNPs:", nrow(mvdat$exposure_beta), "\n")
cat("Exposures (columns):", paste(colnames(mvdat$exposure_beta), collapse = " | "), "\n")
cat("Outcome rows:", length(mvdat$outcome_beta), "\n\n")

# Save the harmonised object properly
saveRDS(mvdat, file.path(base_dir, "data", paste0("mv_harmonised_", run_tag, ".rds")))

# OPTIONAL: also export a 'long' data.frame version for inspection
mvdat_long <- data.frame(
  SNP = rownames(mvdat$exposure_beta),
  outcome_beta = mvdat$outcome_beta,
  outcome_se   = mvdat$outcome_se,
  outcome_pval = mvdat$outcome_pval
)

# add each exposure beta/se as columns
for(j in seq_len(ncol(mvdat$exposure_beta))){
  ex_name <- colnames(mvdat$exposure_beta)[j]
  mvdat_long[[paste0(ex_name, "_beta")]] <- mvdat$exposure_beta[, j]
  mvdat_long[[paste0(ex_name, "_se")]]   <- mvdat$exposure_se[, j]
}

readr::write_csv(mvdat_long, file.path(base_dir, "data", paste0("mv_harmonised_long_", run_tag, ".csv")))

# --------------------------
# STEP 5 — RUN MVMR (DIRECT EFFECTS)
# --------------------------
cat(">> Running mv_multiple (MVMR)...\n")
res_mvmr <- TwoSampleMR::mv_multiple(mvdat)

cat("\n===== MVMR RESULTS (Direct effects) =====\n")
print(res_mvmr)
cat("========================================\n\n")

readr::write_csv(res_mvmr$result,
                 file.path(base_dir, "results",
                           paste0("MVMR_results_", run_tag, ".csv")))


# Optional: also export to XLSX
openxlsx::write.xlsx(
  list(
    mvmr_results = res_mvmr,
    instruments = exp_mv,
    outcome = out_dat
  ),
  file = file.path(base_dir, "results", paste0("MVMR_bundle_", run_tag, ".xlsx")),
  overwrite = TRUE
)

# --------------------------
# ==========================
# STEP 6 — FOREST PLOT (OK)
# ==========================
cat(">> Plotting forest of direct effects...\n")

# IMPORTANT: use the DATA FRAME, not the list
res_df <- res_mvmr$result

plot_df <- res_df %>%
  dplyr::mutate(
    lo = b - 1.96 * se,
    hi = b + 1.96 * se
  )

p_forest <- ggplot(plot_df,
                   aes(x = b,
                       y = reorder(exposure, b))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2.6) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  labs(
    title = "MVMR direct effects: Education + Pack-years + Hannum → Oral cancer",
    x = "Direct effect (β) with 95% CI",
    y = ""
  ) +
  theme_bw()

ggsave(
  filename = file.path(base_dir, "figures",
                       paste0("MVMR_forest_", run_tag, ".png")),
  plot = p_forest,
  width = 9,
  height = 4.8,
  dpi = 300
)


# --------------------------
# STEP 7 — UNIVARIABLE MR (for comparison) ✅ FIX
# --------------------------
cat(">> Univariable MR for each exposure -> Y (same run)...\n")

uv_list <- lapply(exposure_ids, function(eid){
  
  # 1) Extract instruments (no clump args here, your version doesn't support them)
  exp_uv <- TwoSampleMR::extract_instruments(outcomes = eid, p1 = 5e-8)
  
  # 2) Clump separately (here you can set r2/kb)
  exp_uv <- TwoSampleMR::clump_data(exp_uv, clump_r2 = 0.001, clump_kb = 10000)
  
  # 3) Extract outcome associations
  out_uv <- TwoSampleMR::extract_outcome_data(snps = exp_uv$SNP, outcomes = Y_id)
  
  # 4) Harmonise + MR
  dat_uv <- TwoSampleMR::harmonise_data(exp_uv, out_uv)
  mr_uv  <- TwoSampleMR::mr(dat_uv)
  
  mr_uv$exposure_id <- eid
  mr_uv
})

res_uv <- dplyr::bind_rows(uv_list)

readr::write_csv(res_uv,
                 file.path(base_dir, "results",
                           paste0("Univariable_MR_", run_tag, ".csv")))

cat("\n===== UNIVARIABLE MR (for comparison) =====\n")
print(res_uv)
cat("===========================================\n\n")

# --------------------------
# STEP 8 — DONE
# --------------------------
cat("DONE. Outputs saved under:\n", base_dir, "\n")
sink()

beepr::beep(2)
