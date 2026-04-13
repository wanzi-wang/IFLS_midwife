# =====================================================================
# lib/io.R — Raw-data I/O helpers (read Stata .dta, cache .rds)
#
# Single place for Stata-label handling, ID standardization, and cache
# semantics. Downstream stages use `read_raw()` not `haven::read_dta()`.
# =====================================================================

#' Read an IFLS Stata file from the raw root.
#'
#' Reads, zaps value-labels (keeps the underlying integer/numeric), and
#' casts the common identifier columns to character so cross-wave joins
#' don't silently fail on type drift.
#'
#' @param relpath Path relative to `ifls_root()`.
#' @return A tibble.
read_ifls_dta <- function(relpath) {
  path <- ifls_path(relpath)
  if (!file.exists(path)) {
    stop(sprintf("Raw file missing: %s", path))
  }
  df <- haven::read_dta(path)
  df <- haven::zap_labels(df)
  id_cols <- intersect(
    c("pidlink", "hhid93", "hhid97", "hhid00", "hhid07", "hhid14",
      "hhid93_9", "hhid97_9", "hhid00_9", "hhid07_9", "hhid14_9",
      "commid", "commid93", "commid97", "commid00", "commid07", "commid14",
      "pid", "pid93", "pid97", "pid00", "pid07", "pid14",
      "ea", "ea93", "ea97", "ea00", "ea07", "ea14"),
    names(df)
  )
  for (col in id_cols) {
    df[[col]] <- as.character(df[[col]])
  }
  tibble::as_tibble(df)
}

#' Write a tibble to `data/intermediate/<subdir>/<name>.rds`.
#'
#' Always `compress = "xz"` — IFLS intermediates are small enough that the
#' CPU cost is trivial and the on-disk footprint matters for Time Machine
#' and backups.
#'
#' @param obj Object to save.
#' @param subdir Subdirectory under `data/intermediate/`.
#' @param name File basename (no extension).
write_intermediate <- function(obj, subdir, name) {
  d <- int_subdir(subdir)
  saveRDS(obj, file.path(d, paste0(name, ".rds")), compress = "xz")
  invisible(NULL)
}

#' Read an intermediate `.rds` back.
read_intermediate <- function(subdir, name) {
  readRDS(file.path(int_subdir(subdir), paste0(name, ".rds")))
}

#' Append one timestamped line to a stage log at
#' `data/intermediate/_logs/<stage>.log`. Logs are gitignored along with
#' the rest of `data/intermediate/`.
log_stage <- function(stage, msg) {
  d <- here::here("data", "intermediate", "_logs")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  line <- sprintf("[%s] %s: %s\n", format(Sys.time(), "%F %T"), stage, msg)
  cat(line, file = file.path(d, paste0(stage, ".log")), append = TRUE)
  cat(line)
}
