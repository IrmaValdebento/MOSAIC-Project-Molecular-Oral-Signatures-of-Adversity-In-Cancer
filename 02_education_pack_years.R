rm(list = ls())
############################################################
# TWO-STEP MENDELIAN RANDOMIZATION
# Step 1: Education -> Pack-years (ukb-b-10831)
# Step 2: Pack-years -> Oral cavity cancer (ieu-b-4961)
# Author: Irma Valdebenito
############################################################

# --------------------------
# STEP 0 — PACKAGES
# --------------------------
pkgs <- c("TwoSampleMR","ieugwasr","dplyr","readr","ggplot2","beepr","openxlsx","tibble")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install)

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(readr)
library(ggplot2)
library(beepr)
library(openxlsx)
library(tibble)

# --------------------------
# STEP 1 — FOLDERS (GitHub test folder)
# --------------------------
base_dir <- "YOUR_FOLDER"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(base_dir, paste0("X1_Education__M_PackYears__Y_OralCancer__", run_tag))

step1_dir <- file.path(run_dir, "step1_X_to_M_education_to_packyears")
step2_dir <- file.path(run_dir, "step2_M_to_Y_packyears_to_oral_cancer")

dir.create(step1_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(step2_dir, recursive = TRUE, showWarnings = FALSE)

for (d in c("data_rds","tables_csv","figures","logs")) {
  dir.create(file.path(step1_dir, d), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(step2_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# logging
log_file <- file.path(run_dir, "run_log.txt")
sink(log_file, split = TRUE)
on.exit({ try(sink(), silent = TRUE) }, add = TRUE)

cat("Run started:", Sys.time(), "\n")
cat("Run folder :", run_dir, "\n\n")
beepr::beep(1)

# --------------------------
# STEP 2 — IDs
# --------------------------
X_id <- "ebi-a-GCST90029013"   # Education (years of education)
M_id <- "ukb-b-10831"          # Pack-years
Y_id <- "ieu-b-4961"           # Oral cavity cancer

# --------------------------
# STEP 3 — AUTH (OpenGWAS JWT)
# --------------------------
cat("AUTH: OpenGWAS authentication...\n")

jwt <- trimws(Sys.getenv("YOUR_TOKEN"))

if (jwt == "") {
  cat("\nNo OPENGWAS_JWT found.\nPaste your OpenGWAS JWT token and press ENTER:\n")
  jwt <- trimws(readline(prompt = "> "))
  if (jwt == "") stop("No token provided. Stopping.")
  Sys.setenv(OPENGWAS_JWT = jwt)
  cat("Token loaded into session (not saved to disk).\n")
} else {
  cat("Token found in environment.\n")
}

# Sanity check: JWT has 3 dot-separated parts
if (length(strsplit(jwt, "\\.")[[1]]) != 3) stop("Token does not look like a JWT (expected 3 dot-separated parts).")

# Live API check
cat("Testing API with ieugwasr::user()...\n")
u <- try(ieugwasr::user(), silent = TRUE)
if (inherits(u, "try-error")) stop("Token present but API call failed (401/403). Generate a NEW token and retry.")
cat("✅ Auth confirmed.\n\n")
beepr::beep(2)

# Quick ID checks
xi <- ieugwasr::gwasinfo(X_id); if (nrow(xi) == 0) stop("❌ X_id not found.")
mi <- ieugwasr::gwasinfo(M_id); if (nrow(mi) == 0) stop("❌ M_id not found.")
yi <- ieugwasr::gwasinfo(Y_id); if (nrow(yi) == 0) stop("❌ Y_id not found.")

cat("X trait:", xi$trait[1], "\n")
cat("M trait:", mi$trait[1], "\n")
cat("Y trait:", yi$trait[1], "\n\n")
beepr::beep(2)

# --------------------------
# STEP 4 — STEP 1 MR: X → M
# --------------------------
cat("===== STEP 1: Education -> Pack-years =====\n")

exp_X <- extract_instruments(X_id, clump = TRUE)
out_M <- extract_outcome_data(
  snps = exp_X$SNP, outcomes = M_id,
  proxies = TRUE, rsq = 0.8, align_alleles = 1,
  palindromes = 1, maf_threshold = 0.3
)

dat_XM <- harmonise_data(exp_X, out_M)
mr_XM  <- mr(dat_XM)

saveRDS(dat_XM, file.path(step1_dir, "data_rds", "harm_XM.rds"))
write_csv(as_tibble(mr_XM), file.path(step1_dir, "tables_csv", "mr_XM.csv"))

# Sensitivity
ss_XM  <- mr_singlesnp(dat_XM)
loo_XM <- mr_leaveoneout(dat_XM)

ggsave(file.path(step1_dir, "figures", "forest_X_to_M.png"),
       mr_forest_plot(ss_XM)[[1]], width=7, height=10, dpi=300)
ggsave(file.path(step1_dir, "figures", "scatter_X_to_M.png"),
       mr_scatter_plot(mr_XM, dat_XM)[[1]], width=7, height=6, dpi=300)
ggsave(file.path(step1_dir, "figures", "funnel_X_to_M.png"),
       mr_funnel_plot(ss_XM)[[1]], width=7, height=6, dpi=300)
ggsave(file.path(step1_dir, "figures", "leaveoneout_X_to_M.png"),
       mr_leaveoneout_plot(loo_XM)[[1]], width=7, height=10, dpi=300)

cat("STEP 1 complete.\n\n")
beepr::beep(3)

# --------------------------
# STEP 5 — STEP 2 MR: M → Y
# --------------------------
cat("===== STEP 2: Pack-years -> Oral cancer =====\n")

exp_M <- extract_instruments(M_id, clump = TRUE)
out_Y <- extract_outcome_data(
  snps = exp_M$SNP, outcomes = Y_id,
  proxies = TRUE, rsq = 0.8, align_alleles = 1,
  palindromes = 1, maf_threshold = 0.3
)

dat_MY <- harmonise_data(exp_M, out_Y)
mr_MY  <- mr(dat_MY)

saveRDS(dat_MY, file.path(step2_dir, "data_rds", "harm_MY.rds"))
write_csv(as_tibble(mr_MY), file.path(step2_dir, "tables_csv", "mr_MY.csv"))

# Sensitivity
ss_MY  <- mr_singlesnp(dat_MY)
loo_MY <- mr_leaveoneout(dat_MY)

ggsave(file.path(step2_dir, "figures", "forest_M_to_Y.png"),
       mr_forest_plot(ss_MY)[[1]], width=7, height=10, dpi=300)
ggsave(file.path(step2_dir, "figures", "scatter_M_to_Y.png"),
       mr_scatter_plot(mr_MY, dat_MY)[[1]], width=7, height=6, dpi=300)
ggsave(file.path(step2_dir, "figures", "funnel_M_to_Y.png"),
       mr_funnel_plot(ss_MY)[[1]], width=7, height=6, dpi=300)
ggsave(file.path(step2_dir, "figures", "leaveoneout_M_to_Y.png"),
       mr_leaveoneout_plot(loo_MY)[[1]], width=7, height=10, dpi=300)

cat("STEP 2 complete.\n\n")
beepr::beep(3)

# --------------------------
# DONE
# --------------------------
cat("✅ Two-step MR finished:", Sys.time(), "\n")
cat("Outputs saved in:", run_dir, "\n")
beepr::beep(8)

