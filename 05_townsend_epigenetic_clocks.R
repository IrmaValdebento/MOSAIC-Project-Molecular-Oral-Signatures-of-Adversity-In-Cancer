rm(list = ls())

# ============================================================
# TWO-STEP MENDELIAN RANDOMIZATION (MEDIATION)
# TOWNSEND -> EPIGENETIC CLOCKS (PhenoAge, Hannum) -> ORAL CANCER
#
# X  : Townsend deprivation index            [ukb-b-10011]
# M1 : DNA methylation PhenoAge acceleration [ebi-a-GCST90014292]
# M2 : DNA methylation Hannum age accel      [ebi-a-GCST90014301]
# Y  : Oral cavity cancer                    [ieu-b-4961]
# AUTHOR: IRMA VALDEBENITO
# ============================================================

# --------------------------
# PACKAGES
# --------------------------
pkgs <- c("TwoSampleMR","ieugwasr","dplyr","readr","ggplot2","openxlsx","tibble","beepr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# --------------------------
# FOLDERS (GitHub)
# --------------------------
base_dir <- "YOUR_FOLDER"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(base_dir, paste0("X3_Townsend__M_EpiClocks(PhenoAge_Hannum)__Y_OralCancer__", run_tag))
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(run_dir, "01_step1_Townsend_to_Clock"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "02_step2_Clock_to_OralCancer"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "03_summary"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(run_dir, "logs", paste0("run_log_", run_tag, ".txt"))
sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

cat("Run started:", as.character(Sys.time()), "\n")
cat("Run folder :", run_dir, "\n")
cat("TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("ieugwasr   :", as.character(packageVersion("ieugwasr")), "\n\n")
try(beepr::beep(1), silent = TRUE)

# --------------------------
# AUTH (OpenGWAS JWT)  ✅ no token hardcoded
# --------------------------
cat("AUTH: OpenGWAS authentication...\n")
jwt <- trimws(Sys.getenv("YOUR_TOKEN"))

if (jwt == "") {
  cat("\nNo OPENGWAS_JWT found in environment.\nPaste your OpenGWAS JWT token and press ENTER:\n")
  jwt <- trimws(readline(prompt = "> "))
  if (jwt == "") stop("No token provided. Stopping.")
  Sys.setenv(OPENGWAS_JWT = jwt)
  cat("Token loaded into session (not saved to disk).\n")
} else {
  cat("Token found in environment.\n")
}

if (length(strsplit(jwt, "\\.")[[1]]) != 3) stop("Token does not look like a JWT (expected 3 dot-separated parts).")

cat("Testing API with ieugwasr::user()...\n")
u <- try(ieugwasr::user(), silent = TRUE)
if (inherits(u, "try-error")) stop("Token present but API call failed (401/403). Generate a NEW token and retry.")
cat("✅ Auth confirmed.\n\n")
try(beepr::beep(2), silent = TRUE)

# --------------------------
# IDS
# --------------------------
X_id <- "ukb-b-10011"          # Townsend deprivation index
Y_id <- "ieu-b-4961"           # Oral cavity cancer

clocks <- tibble::tibble(
  clock = c("PhenoAge", "Hannum"),
  M_id  = c("ebi-a-GCST90014292", "ebi-a-GCST90014301")
)

# quick checks (optional, but nice)
cat("Checking GWAS IDs...\n")
stopifnot(nrow(ieugwasr::gwasinfo(X_id)) > 0)
stopifnot(nrow(ieugwasr::gwasinfo(Y_id)) > 0)
for (j in seq_len(nrow(clocks))) stopifnot(nrow(ieugwasr::gwasinfo(clocks$M_id[j])) > 0)
cat("✅ IDs exist.\n\n")

# --------------------------
# HELPERS
# --------------------------
safe_write_csv <- function(df, path){
  tryCatch(readr::write_csv(df, path),
           error = function(e) message("write_csv error: ", e$message))
}

safe_save_plot <- function(p, path, w=8, h=6){
  tryCatch(ggsave(filename = path, plot = p, width = w, height = h, dpi = 300),
           error = function(e) message("ggsave error: ", e$message))
}

get_ivw <- function(mr_res){
  if (is.null(mr_res) || nrow(mr_res) == 0) return(NULL)
  mr_res %>% dplyr::filter(method == "Inverse variance weighted") %>% dplyr::slice(1)
}

run_mr_block <- function(exposure_id, outcome_id, label, outdir){
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "data_rds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "tables_csv"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
  
  cat("\n=============================\n")
  cat("RUN:", label, "\nExposure:", exposure_id, "\nOutcome:", outcome_id, "\n")
  cat("=============================\n")
  
  exp_dat <- tryCatch(
    TwoSampleMR::extract_instruments(outcomes = exposure_id, p1 = 5e-8, clump = TRUE),
    error = function(e){ cat("extract_instruments error:", e$message, "\n"); return(NULL) }
  )
  if (is.null(exp_dat) || nrow(exp_dat) < 3) {
    cat("Not enough instruments (<3). Skipping.\n")
    return(list(ok=FALSE, reason="too_few_instruments", mr=NULL, n_snps=NA_integer_))
  }
  
  out_dat <- tryCatch(
    TwoSampleMR::extract_outcome_data(
      snps = exp_dat$SNP, outcomes = outcome_id,
      proxies = TRUE, rsq = 0.8, align_alleles = 1,
      palindromes = 1, maf_threshold = 0.3
    ),
    error = function(e){ cat("extract_outcome_data error:", e$message, "\n"); return(NULL) }
  )
  if (is.null(out_dat) || nrow(out_dat) == 0) {
    cat("No outcome data returned. Skipping.\n")
    return(list(ok=FALSE, reason="no_outcome_data", mr=NULL, n_snps=NA_integer_))
  }
  
  dat <- tryCatch(
    TwoSampleMR::harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat),
    error = function(e){ cat("harmonise_data error:", e$message, "\n"); return(NULL) }
  )
  if (is.null(dat) || nrow(dat) < 3) {
    cat("Not enough harmonised SNPs (<3). Skipping.\n")
    return(list(ok=FALSE, reason="too_few_harmonised", mr=NULL, n_snps=NA_integer_))
  }
  
  mr_res <- TwoSampleMR::mr(dat)
  het    <- tryCatch(TwoSampleMR::mr_heterogeneity(dat), error=function(e) NULL)
  pleio  <- tryCatch(TwoSampleMR::mr_pleiotropy_test(dat), error=function(e) NULL)
  single <- tryCatch(TwoSampleMR::mr_singlesnp(dat), error=function(e) NULL)
  loo    <- tryCatch(TwoSampleMR::mr_leaveoneout(dat), error=function(e) NULL)
  
  saveRDS(dat, file.path(outdir, "data_rds", paste0(label, "_harmonised.rds")))
  safe_write_csv(dat,    file.path(outdir, "tables_csv", paste0(label, "_harmonised.csv")))
  safe_write_csv(mr_res, file.path(outdir, "tables_csv", paste0(label, "_mr.csv")))
  if (!is.null(het))   safe_write_csv(het,   file.path(outdir, "tables_csv", paste0(label, "_heterogeneity.csv")))
  if (!is.null(pleio)) safe_write_csv(pleio, file.path(outdir, "tables_csv", paste0(label, "_pleiotropy_egger.csv")))
  if (!is.null(single))safe_write_csv(single,file.path(outdir, "tables_csv", paste0(label, "_singleSNP.csv")))
  if (!is.null(loo))   safe_write_csv(loo,   file.path(outdir, "tables_csv", paste0(label, "_leaveoneout.csv")))
  
  # Plots
  p_scatter <- tryCatch(TwoSampleMR::mr_scatter_plot(mr_res, dat)[[1]], error=function(e) NULL)
  if (!is.null(p_scatter)) safe_save_plot(p_scatter, file.path(outdir, "figures", paste0(label, "_scatter.png")), 8, 6)
  
  p_forest <- tryCatch({ if (!is.null(single) && nrow(single) > 0) TwoSampleMR::mr_forest_plot(single)[[1]] else NULL }, error=function(e) NULL)
  if (!is.null(p_forest)) safe_save_plot(p_forest, file.path(outdir, "figures", paste0(label, "_forest.png")), 9, 7)
  
  p_funnel <- tryCatch({ if (!is.null(single) && nrow(single) > 0) TwoSampleMR::mr_funnel_plot(single)[[1]] else NULL }, error=function(e) NULL)
  if (!is.null(p_funnel)) safe_save_plot(p_funnel, file.path(outdir, "figures", paste0(label, "_funnel.png")), 8, 7)
  
  p_loo <- tryCatch({ if (!is.null(loo) && nrow(loo) > 0) TwoSampleMR::mr_leaveoneout_plot(loo)[[1]] else NULL }, error=function(e) NULL)
  if (!is.null(p_loo)) safe_save_plot(p_loo, file.path(outdir, "figures", paste0(label, "_leaveoneout.png")), 9, 7)
  
  return(list(ok=TRUE, reason="ok", mr=mr_res, n_snps=nrow(dat), het=het, pleio=pleio))
}

# --------------------------
# MAIN: run two-step per clock
# --------------------------
summary_rows <- list()

for (i in seq_len(nrow(clocks))) {
  clk  <- clocks$clock[i]
  M_id <- clocks$M_id[i]
  
  cat("\n#############################################\n")
  cat("CLOCK:", clk, "\nMediator ID:", M_id, "\n")
  cat("#############################################\n")
  
  # Step 1: Townsend -> Clock
  out1 <- file.path(run_dir, "01_step1_Townsend_to_Clock", clk)
  step1 <- run_mr_block(
    exposure_id = X_id, outcome_id = M_id,
    label = paste0("STEP1_Townsend_to_", clk),
    outdir = out1
  )
  
  # Step 2: Clock -> Cancer
  out2 <- file.path(run_dir, "02_step2_Clock_to_OralCancer", clk)
  step2 <- run_mr_block(
    exposure_id = M_id, outcome_id = Y_id,
    label = paste0("STEP2_", clk, "_to_OralCancer"),
    outdir = out2
  )
  
  s1_ivw <- get_ivw(step1$mr)
  s2_ivw <- get_ivw(step2$mr)
  
  summary_rows[[i]] <- tibble::tibble(
    clock = clk,
    clock_id = M_id,
    
    step1_ok = step1$ok,
    step1_nsnps = step1$n_snps,
    step1_ivw_b = if (!is.null(s1_ivw)) s1_ivw$b else NA_real_,
    step1_ivw_p = if (!is.null(s1_ivw)) s1_ivw$pval else NA_real_,
    
    step2_ok = step2$ok,
    step2_nsnps = step2$n_snps,
    step2_ivw_b = if (!is.null(s2_ivw)) s2_ivw$b else NA_real_,
    step2_ivw_p = if (!is.null(s2_ivw)) s2_ivw$pval else NA_real_
  )
  
  try(beepr::beep(2), silent = TRUE)
}

summary_df <- dplyr::bind_rows(summary_rows)

# Save summary
safe_write_csv(summary_df, file.path(run_dir, "03_summary", "two_step_summary_Townsend_PhenoAge_Hannum.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Summary")
openxlsx::writeData(wb, "Summary", summary_df)
openxlsx::saveWorkbook(wb, file.path(run_dir, "03_summary", "two_step_summary_Townsend_PhenoAge_Hannum.xlsx"), overwrite = TRUE)

cat("\nDONE ✅\nSaved in:\n", run_dir, "\n")
cat("Run finished:", as.character(Sys.time()), "\n")
try(beepr::beep(8), silent = TRUE)

# End
