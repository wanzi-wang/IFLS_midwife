# =====================================================================
# Master orchestrator — Indonesia Village Midwife Programme upgrade
# Author: Wanzi Wang
#
# Usage:  Rscript code/R/00_run_all.R
# Prereq: IFLS_DATA_ROOT env var (default: old/IFLS under project root)
#
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
})

# Sanity: the raw root must exist. `ifls_root()` hard-fails otherwise,
# but we also want to verify before sourcing stages so error messages
# don't mention a stage the user hasn't gotten to yet.
source(here::here("code", "R", "lib", "paths.R"))
invisible(ifls_root())

# ---------------------------------------------------------------------
# Pipeline stages. Commented stages are not yet implemented.
# ---------------------------------------------------------------------
t0 <- Sys.time()

source(here::here("code", "R", "01_load_ifls.R"))
source(here::here("code", "R", "02_birth_roster.R"))
source(here::here("code", "R", "03_midwife_rollout.R"))
source(here::here("code", "R", "04_outcomes.R"))
source(here::here("code", "R", "04b_mediators.R"))
source(here::here("code", "R", "05_analysis_frame.R"))

source(here::here("code", "R", "06_sample_weights_indices.R"))
source(here::here("code", "R", "07_balance_descriptives.R"))
source(here::here("code", "R", "08_main_regressions.R"))
source(here::here("code", "R", "09_event_study.R"))
source(here::here("code", "R", "10_sibling_fe.R"))
source(here::here("code", "R", "11_robustness.R"))
source(here::here("code", "R", "12_heterogeneity_finalize.R"))
source(here::here("code", "R", "13_multiple_inference.R"))
source(here::here("code", "R", "14_mechanisms.R"))
source(here::here("code", "R", "16_appendix_figures.R"))

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
message(sprintf("Pipeline finished in %s s.", elapsed))
