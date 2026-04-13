# =====================================================================
# 01_load_ifls.R — Raw IFLS .dta → typed .rds cache
#
# Single place where raw IFLS files are touched. Ingests every wave and
# module enumerated in `ifls_file_catalog()`, zaps Stata value-labels,
# standardizes identifier types, and caches one .rds per logical file
# under `data/intermediate/raw/`.
#
# Downstream stages (02..05) read only from `data/intermediate/raw/`,
# never from `$IFLS_DATA_ROOT`.
#
# Pipeline contract:
#   - No sample restrictions here.
#   - No variable renaming here (except ID-type coercion).
#   - No joins here.
# =====================================================================

set.seed(20260412)

library(here)
library(haven)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
catalog <- ifls_file_catalog()

stage01_load_one <- function(wave, name, relpath) {
  out_key <- paste0(wave, "__", name)
  out_path <- file.path(int_subdir("raw"), paste0(out_key, ".rds"))
  if (file.exists(out_path)) {
    return(tibble::tibble(
      wave = wave, name = name, rows = NA_integer_,
      cols = NA_integer_, status = "cached"
    ))
  }
  df <- read_ifls_dta(relpath)
  saveRDS(df, out_path, compress = "xz")
  tibble::tibble(
    wave = wave, name = name,
    rows = nrow(df), cols = ncol(df), status = "loaded"
  )
}

log_stage("stage01", sprintf(
  "begin: %d files in catalog; cache -> %s",
  nrow(catalog), int_subdir("raw")
))

summary <- purrr::pmap_dfr(
  list(catalog$wave, catalog$name, catalog$relpath),
  stage01_load_one
)

write_intermediate(summary, "raw", "_manifest")

log_stage("stage01", sprintf(
  "done: %d loaded, %d cached, total %d files",
  sum(summary$status == "loaded"),
  sum(summary$status == "cached"),
  nrow(summary)
))

message("--- Stage 01 summary ---")
print(summary, n = Inf)
