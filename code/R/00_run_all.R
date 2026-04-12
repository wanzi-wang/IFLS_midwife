# =====================================================================
# Master orchestrator — Indonesia Village Midwife Program
# Author: Wanzi Wang
#
# Runs the full data-construction and analysis pipeline end to end.
# Individual stage scripts to be written in the data-pipeline plan.
#
# Usage:  Rscript code/R/00_run_all.R
# Prereq: set IFLS_DATA_ROOT env var (default: ~/Academic/Thesis_Data_IFLS/old/IFLS)
# =====================================================================

set.seed(20260412)

# --- Environment ------------------------------------------------------
IFLS_DATA_ROOT <- Sys.getenv(
  "IFLS_DATA_ROOT",
  unset = "/Users/wwz/Academic/Thesis_Data_IFLS/old/IFLS"
)
stopifnot(dir.exists(IFLS_DATA_ROOT))

PROJECT_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."))

# --- Pipeline stages (to be implemented) ------------------------------
# source("code/R/01_load_ifls.R")
# source("code/R/02_construct_treatment.R")
# source("code/R/03_birth_roster.R")
# source("code/R/04_outcomes.R")
# source("code/R/05_merge_analysis_df.R")
# source("code/R/06_descriptives.R")
# source("code/R/07_main_regressions.R")
# source("code/R/08_event_study.R")
# source("code/R/09_robustness.R")
# source("code/R/10_heterogeneity.R")
# source("code/R/11_tables.R")
# source("code/R/12_figures.R")

message("Pipeline scaffolding loaded. Stage scripts are stubs — see data-pipeline plan.")
