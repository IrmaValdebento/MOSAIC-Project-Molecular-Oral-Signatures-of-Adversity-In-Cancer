# MOSAIC Project: Molecular Oral Signatures of Adversity In Cancer

**Grant Application:** MOSAIC - Letter of Intent (2026)

## Overview
This repository contains the analytic code for the preliminary data presented in the MOSAIC Letter of Intent. Our research identifies divergent biological pathways for social adversity in oral cancer, distinguishing between "biological weathering" (direct aging) and behavioral risk accumulation.

## Repository Structure & Analysis Workflow

The analysis is divided into 5 scripts covering Univariable and Multivariable Mendelian Randomization (MVMR):

### 1. Education Pathway (The "Biological Shield")
* **`01_education_epigenetic_clocks.R`**: Primary analysis. Demonstrates that higher educational attainment causally reduces epigenetic age acceleration (e.g., PhenoAge, Hannum), acting as a protective biological shield.
* **`02_education_pack_years.R`**: Confirming the behavioral link. Shows that education also reduces cumulative smoking exposure.
* **`03_mvmr_education_packyears_hannum.R`**: Multivariable MR analysis. Crucial step that disentangles the effects, suggesting that education's protective effect on biological aging persists even after adjusting for smoking behavior.

### 2. Townsend Pathway (The "Behavioral Driver")
* **`04_townsend_pack_years.R`**: Demonstrates that material deprivation (Townsend Index) drastically increases cumulative smoking risk (Behavioral pathway).
* **`05_townsend_epigenetic_clocks.R`**: Negative Control. Shows that material deprivation *does not* significantly accelerate epigenetic aging directly when behavioral mediators are not considered, highlighting the specificity of the biological pathway found for Education.

## Requirements
* R version 4.x
* Packages: `TwoSampleMR`, `ieugwasr`, `tidyverse`, `MendelianRandomization`

## Contact
For questions regarding the analysis or reproducibility, please contact: [Tu Correo]
