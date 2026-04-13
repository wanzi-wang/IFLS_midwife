# =====================================================================
# lib/paths.R — Path resolution for IFLS raw data and pipeline outputs
#
# Single source of truth for filesystem layout. Every stage script sources
# this file. No other file should hardcode IFLS subpaths.
# =====================================================================

#' Resolve the IFLS raw-data root directory.
#'
#' Reads the `IFLS_DATA_ROOT` env var, falling back to the default on the
#' author's machine. Fails fast if the directory does not exist.
#'
#' @return Absolute path to the IFLS raw-data root (contains W1..W5).
ifls_root <- function() {
  root <- Sys.getenv(
    "IFLS_DATA_ROOT",
    unset = "/Users/wwz/Academic/Thesis_Data_IFLS/old/IFLS"
  )
  if (!dir.exists(root)) {
    stop(sprintf("IFLS_DATA_ROOT '%s' does not exist. Set env var.", root))
  }
  normalizePath(root)
}

#' Build an absolute path under the IFLS raw-data root.
#'
#' @param ... Path components appended to `ifls_root()`.
#' @return Absolute path; does NOT check existence (callers decide).
ifls_path <- function(...) file.path(ifls_root(), ...)

#' Absolute path under the project's `data/intermediate/` cache.
#'
#' Creates the directory on first call.
int_path <- function(...) {
  d <- here::here("data", "intermediate")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, ...)
}

#' Per-stage cache directory under `data/intermediate/`.
#'
#' Stages that emit multiple raw-derived files (e.g. 01_load_ifls) write here.
int_subdir <- function(sub) {
  d <- here::here("data", "intermediate", sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Catalog of raw IFLS files used by stage 01.
#'
#' Each row: wave, module, logical name (stable across waves), and relative
#' path under `ifls_root()`. Verified against the actual raw data layout on
#' 2026-04-13. Add rows here — never sprinkle paths across stage scripts.
ifls_file_catalog <- function() {
  tibble::tribble(
    ~wave, ~module,       ~name,         ~relpath,
    # --- W1 (1993): individual roster + household chars + mother vars ---
    "W1",  "individual",  "roster_main", "W1/hh93dta/bukkar1.dta",
    "W1",  "individual",  "roster_alt",  "W1/hh93dta/bukkar2.dta",
    "W1",  "individual",  "height",      "W1/hh93dta/bukcca2.dta",
    "W1",  "individual",  "smoke",       "W1/hh93dta/buk3km1.dta",
    "W1",  "community",   "sar",         "W1/cf93dta/bkps_cf.dta",
    # --- W2 (1997): individual roster + tracking + community ---
    "W2",  "individual",  "roster_main", "W2/hh97dta/hh97bk/bk_ar1.dta",
    "W2",  "individual",  "roster_hh",   "W2/hh97dta/hh97bk/bk_ar0.dta",
    "W2",  "tracking",    "ptrack",      "W2/hh97dta/hh97bk/ptrack.dta",
    "W2",  "tracking",    "htrack",      "W2/hh97dta/hh97bk/htrack.dta",
    "W2",  "community",   "sar",         "W2/cf97dta/sar_v12.dta",
    "W2",  "community",   "pkk",         "W2/cf97dta/pkk_v12.dta",
    # --- W3 (2000) ---
    "W3",  "individual",  "roster_main", "W3/hh00_all_dta/bk_ar1.dta",
    "W3",  "individual",  "roster_hh",   "W3/hh00_all_dta/bk_ar0.dta",
    "W3",  "tracking",    "ptrack",      "W3/hh00_all_dta/ptrack.dta",
    "W3",  "tracking",    "htrack",      "W3/hh00_all_dta/htrack.dta",
    "W3",  "community",   "sar",         "W3/cf00_all_dta/sar.dta",
    "W3",  "community",   "pkk",         "W3/cf00_all_dta/pkk.dta",
    # --- W4 (2007): tracking + community only (outcomes come from W5) ---
    "W4",  "tracking",    "ptrack",      "W4/hh07_all_dta/ptrack.dta",
    "W4",  "tracking",    "htrack",      "W4/hh07_all_dta/htrack.dta",
    "W4",  "individual",  "roster_main", "W4/hh07_all_dta/bk_ar1.dta",
    "W4",  "individual",  "roster_hh",   "W4/hh07_all_dta/bk_ar0.dta",
    "W4",  "community",   "sar",         "W4/cf07_all_dta/sar.dta",
    "W4",  "community",   "pkk",         "W4/cf07_all_dta/pkk.dta",
    # --- W5 (2014): tracking + birth-history + outcomes + community ---
    "W5",  "tracking",    "ptrack",      "W5/hh14_all_dta/ptrack.dta",
    "W5",  "tracking",    "htrack",      "W5/hh14_all_dta/htrack.dta",
    "W5",  "individual",  "roster_main", "W5/hh14_all_dta/bk_ar1.dta",
    "W5",  "individual",  "roster_hh",   "W5/hh14_all_dta/bk_ar0.dta",
    "W5",  "individual",  "b3a_cov",     "W5/hh14_all_dta/b3a_cov.dta",
    "W5",  "community",   "sar",         "W5/cf14_all_dta/sar.dta",
    "W5",  "community",   "pkk",         "W5/cf14_all_dta/pkk.dta",
    # W5 outcomes (adult: respondents from book 3b) ---
    "W5",  "outcome",     "health",      "W5/hh14_all_dta/bus_us.dta",
    "W5",  "outcome",     "early_health","W5/hh14_all_dta/b3b_eh.dta",
    "W5",  "outcome",     "depression",  "W5/hh14_all_dta/b3b_kp.dta",
    "W5",  "outcome",     "personality", "W5/hh14_all_dta/b3b_psn.dta",
    # Adult cognition (book 3B, respondents ~15+): Raven-style + math + recall
    "W5",  "outcome",     "cog_adult1",  "W5/hh14_all_dta/b3b_co1.dta",
    "W5",  "outcome",     "cog_adultb",  "W5/hh14_all_dta/b3b_cob.dta",
    "W5",  "outcome",     "math_a",      "W5/hh14_all_dta/b3b_ma1.dta",
    "W5",  "outcome",     "math_b",      "W5/hh14_all_dta/b3b_ma2.dta",
    "W5",  "outcome",     "recall_1",    "W5/hh14_all_dta/b3b_kk1.dta",
    "W5",  "outcome",     "recall_2",    "W5/hh14_all_dta/b3b_kk2.dta",
    "W5",  "outcome",     "recall_3",    "W5/hh14_all_dta/b3b_kk3.dta",
    "W5",  "outcome",     "recall_4",    "W5/hh14_all_dta/b3b_kk4.dta",
    # Children cognition (book EK, respondents <15 in 2014): Raven for kids
    "W5",  "outcome",     "cog_child1",  "W5/hh14_all_dta/ek_ek1.dta",
    "W5",  "outcome",     "cog_child2",  "W5/hh14_all_dta/ek_ek2.dta"
  )
}
