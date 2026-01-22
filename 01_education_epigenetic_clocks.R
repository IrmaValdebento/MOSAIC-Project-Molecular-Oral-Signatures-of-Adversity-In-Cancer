# ============================================================
# TWO-STEP MENDELIAN RANDOMIZATION (MEDIATION)
# EDUCATION -> EPIGENETIC AGING (PhenoAge / Hannum) -> ORAL CANCER
#
##
# X  : Educational attainment (years)        [ebi-a-GCST90029013]
# M1 : PhenoAge acceleration                 [ebi-a-GCST90014292]
# M2 : Hannum age acceleration               [ebi-a-GCST90014301]
# Y  : Oral cavity cancer                    [ieu-b-4961]
#Author: Irma Valdebenito
# ============================================================

rm(list = ls())

# --------------------------
# PACKAGES
# --------------------------
pkgs <- c(
  "TwoSampleMR","ieugwasr","dplyr","readr",
  "ggplot2","openxlsx","tibble","stringr","beepr"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# --------------------------
# OUTPUT FOLDERS
# --------------------------
base_dir <- "your_folder"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

dirs <- c(
  "01_step1_Education_to_Clock",
  "02_step2_Clock_to_OralCancer",
  "03_summary",
  "logs"
)
invisible(lapply(file.path(base_dir, dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

log_file <- file.path(
  base_dir, "logs",
  paste0("run_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)
sink(log_file, split = TRUE)

cat("Run started:", as.character(Sys.time()), "\n")
cat("TwoSampleMR:", as.character(packageVersion("TwoSampleMR")), "\n")
cat("ieugwasr:", as.character(packageVersion("ieugwasr")), "\n\n")

# --------------------------
# AUTHENTICATION (OpenGWAS)
# --------------------------
cat("Authenticating OpenGWAS...\n")
jwt <- trimws(Sys.getenv("TOKEN"))

if (jwt == "") {
  cat("No OPENGWAS_JWT found. Paste token and press ENTER:\n")
  jwt <- trimws(readline(prompt = "> "))
  if (jwt == "") stop("No token provided.")
  Sys.setenv(OPENGWAS_JWT = jwt)
}

if (length(strsplit(jwt, "\\.")[[1]]) != 3) {
  stop("Token does not look like a valid JWT.")
}

u <- try(ieugwasr::user(), silent = TRUE)
if (inherits(u, "try-error")) stop("Token invalid or expired.")
cat("✅ OpenGWAS authentication confirmed.\n\n")
try(beepr::beep(2), silent = TRUE)

# --------------------------
# IDS
# --------------------------
X_id  <- "ebi-a-GCST90029013"   # Education (years)
Y_id  <- "ieu-b-4961"           # Oral cavity cancer

clocks <- tibble(
  clock   = c("PhenoAge", "Hannum"),
  M_id    = c("ebi-a-GCST90014292", "ebi-a-GCST90014301")
)

# --------------------------
# HELPERS
# --------------------------
safe_write_csv <- function(df, path){
  tryCatch(readr::write_csv(df, path),
           error = function(e) message("write_csv error: ", e$message))
}

safe_save_plot <- function(p, path, w=8, h=6){
  tryCatch(ggsave(path, plot = p, width = w, height = h, dpi = 300),
           error = function(e) message("ggsave error: ", e$message))
}

get_ivw <- function(mr_res){
  if (is.null(mr_res)) return(NULL)
  mr_res %>% filter(method == "Inverse variance weighted") %>% slice(1)
}

# --------------------------
# CORE MR FUNCTION
# --------------------------
run_mr <- function(exposure_id, outcome_id, label, outdir){
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  exp_dat <- tryCatch(
    extract_instruments(outcomes = exposure_id, p1 = 5e-8, clump = TRUE),
    error = function(e) NULL
  )
  if (is.null(exp_dat) || nrow(exp_dat) < 3) return(NULL)
  
  out_dat <- tryCatch(
    extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id),
    error = function(e) NULL
  )
  if (is.null(out_dat) || nrow(out_dat) == 0) return(NULL)
  
  dat <- harmonise_data(exp_dat, out_dat)
  if (nrow(dat) < 3) return(NULL)
  
  mr_res <- mr(dat)
  het    <- tryCatch(mr_heterogeneity(dat), error=function(e) NULL)
  pleio  <- tryCatch(mr_pleiotropy_test(dat), error=function(e) NULL)
  single <- tryCatch(mr_singlesnp(dat), error=function(e) NULL)
  loo    <- tryCatch(mr_leaveoneout(dat), error=function(e) NULL)
  
  safe_write_csv(dat,    file.path(outdir, paste0(label, "_harmonised.csv")))
  safe_write_csv(mr_res, file.path(outdir, paste0(label, "_mr.csv")))
  if (!is.null(het))   safe_write_csv(het,   file.path(outdir, paste0(label, "_heterogeneity.csv")))
  if (!is.null(pleio)) safe_write_csv(pleio, file.path(outdir, paste0(label, "_pleiotropy.csv")))
  if (!is.null(single))safe_write_csv(single,file.path(outdir, paste0(label, "_singleSNP.csv")))
  if (!is.null(loo))   safe_write_csv(loo,   file.path(outdir, paste0(label, "_leaveoneout.csv")))
  
  if (!is.null(single)) {
    safe_save_plot(mr_forest_plot(single)[[1]],
                   file.path(outdir, paste0(label, "_forest.png")), 9, 7)
    safe_save_plot(mr_funnel_plot(single)[[1]],
                   file.path(outdir, paste0(label, "_funnel.png")), 8, 7)
  }
  
  safe_save_plot(mr_scatter_plot(mr_res, dat)[[1]],
                 file.path(outdir, paste0(label, "_scatter.png")))
  
  if (!is.null(loo)) {
    safe_save_plot(mr_leaveoneout_plot(loo)[[1]],
                   file.path(outdir, paste0(label, "_leaveoneout.png")), 9, 7)
  }
  
  list(
    nsnps = nrow(dat),
    ivw   = get_ivw(mr_res)
  )
}

# --------------------------
# MAIN ANALYSIS
# --------------------------
summary_list <- list()

for (i in seq_len(nrow(clocks))) {
  
  clk <- clocks$clock[i]
  M_id <- clocks$M_id[i]
  
  cat("Running clock:", clk, "\n")
  
  s1 <- run_mr(
    X_id, M_id,
    paste0("STEP1_Education_to_", clk),
    file.path(base_dir, "01_step1_Education_to_Clock", clk)
  )
  
  s2 <- run_mr(
    M_id, Y_id,
    paste0("STEP2_", clk, "_to_OralCancer"),
    file.path(base_dir, "02_step2_Clock_to_OralCancer", clk)
  )
  
  summary_list[[i]] <- tibble(
    clock = clk,
    step1_nsnps = s1$nsnps,
    step1_ivw_b = s1$ivw$b,
    step1_ivw_p = s1$ivw$pval,
    step2_nsnps = s2$nsnps,
    step2_ivw_b = s2$ivw$b,
    step2_ivw_p = s2$ivw$pval
  )
  
  try(beepr::beep(2), silent = TRUE)
}

summary_df <- bind_rows(summary_list)

safe_write_csv(
  summary_df,
  file.path(base_dir, "03_summary", "two_step_summary_PhenoAge_Hannum.csv")
)

wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
saveWorkbook(
  wb,
  file.path(base_dir, "03_summary", "two_step_summary_PhenoAge_Hannum.xlsx"),
  overwrite = TRUE
)

cat("\nDONE ✅\nSaved in:\n", base_dir, "\n")
cat("Run finished:", as.character(Sys.time()), "\n")

sink()
# End
