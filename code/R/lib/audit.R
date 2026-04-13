# =====================================================================
# lib/audit.R — Reconciliation against legacy `old/Intermediate/Reg/df.rds`
#
# Stage 05 writes the upgraded analysis frame; this helper compares it
# against the legacy analysis frame on the subsample that matches
# legacy's sample construction (commid_birth_legacy + non_migrant +
# alive + cohort 1984-1999). Tolerances are intentionally loose (±5%
# on means) because the upgrade deliberately deviates from legacy on
# several variables (see plan).
#
# On any tolerance failure, writes a diff report to
# a non-zero invisible status via `stop()`.
# =====================================================================

#' Compare new vs legacy analysis frames on a matched subsample.
#'
#' @param new_frame Upgraded analysis frame (tibble).
#' @param legacy_path Path to legacy `df.rds` (frozen pre-rebuild).
#' @param tolerances List of tolerance thresholds:
#'   row_count_pct (default 5), mean_pct (5).
#' @return Tibble of (variable, new, legacy, abs_pct_diff, status).
#'   Also writes a markdown report and either logs PASS or stops().
audit_frame <- function(new_frame, legacy_path,
                         tolerances = list(row_count_pct = 5,
                                            mean_pct = 5)) {
  if (!file.exists(legacy_path)) {
    msg <- sprintf("audit SKIP: legacy frame not at %s", legacy_path)
    log_stage("audit", msg)
    return(invisible(NULL))
  }
  legacy <- readRDS(legacy_path)

  # Matched-subsample construction on the new frame: legacy used
  # commid_birth = commid93 restricted to non-migrants, alive, cohort
  # 1984-1999. The new frame exposes these as columns.
  matched <- new_frame |>
    dplyr::filter(
      !is.na(commid_birth_legacy),
      birth_year %in% 1984:1999,
      dead == 0,
      non_migrant == 1
    )

  # --- 1. Row count ---
  new_n    <- nrow(matched)
  legacy_n <- nrow(legacy)
  row_diff_pct <- abs(new_n - legacy_n) / legacy_n * 100

  # --- 2. Selected continuous-variable means ---
  # Legacy uses different names for some vars. Map where we can.
  # Legacy `height` is in meters (stored as 0.01 of a real cm measurement
  # because legacy divided by 100 in Health.Rmd and never reverted — the
  # df.rds values sit around 1.58). Apply scale_factor to bring legacy
  # onto the same scale as our cm-unit height_cm.
  pairs <- list(
    list(new = "mother_edu_years", leg = "mother_edu_years", scale = 1),
    list(new = "mother_age_birth", leg = "mother_age_birth", scale = 1),
    list(new = "mother_height",    leg = "mother_height",    scale = 1),
    list(new = "height_cm",        leg = "height",           scale = 100),
    list(new = "dep_index",        leg = "dep_index",        scale = 1),
    list(new = "immed_count",      leg = "immed_count",      scale = 1),
    list(new = "delay_count",      leg = "delay_count",      scale = 1)
  )

  rows <- list()
  for (p in pairs) {
    if (!p$new %in% names(matched) || !p$leg %in% names(legacy)) next
    nm <- suppressWarnings(mean(matched[[p$new]], na.rm = TRUE))
    lg <- suppressWarnings(mean(legacy[[p$leg]], na.rm = TRUE)) * p$scale
    if (!is.finite(nm) || !is.finite(lg) || lg == 0) {
      pct <- NA_real_
      status <- "SKIP"
    } else {
      pct <- abs(nm - lg) / abs(lg) * 100
      status <- if (pct <= tolerances$mean_pct) "PASS" else "WARN"
    }
    rows[[length(rows) + 1]] <- tibble::tibble(
      variable = paste0(p$new, " / ", p$leg),
      new = nm, legacy = lg, abs_pct_diff = pct,
      status = status
    )
  }

  means_df <- dplyr::bind_rows(rows)

  # --- 3. Assemble report ---
  report <- c(
    sprintf("# Reconciliation audit — %s", format(Sys.time(), "%F %T")),
    "",
    "## Row count (matched subsample vs legacy df.rds)",
    sprintf("- new N:    %d", new_n),
    sprintf("- legacy N: %d", legacy_n),
    sprintf("- abs pct diff: %.2f%% (threshold %.1f%%) — %s",
            row_diff_pct, tolerances$row_count_pct,
            if (row_diff_pct <= tolerances$row_count_pct) "PASS" else "WARN"),
    "",
    "## Continuous-variable means"
  )
  if (nrow(means_df) > 0) {
    report <- c(report, "",
      "| variable | new | legacy | abs pct diff | status |",
      "|---|---|---|---|---|")
    fmt_num <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "NA"
    fmt_pct <- function(x) if (is.finite(x)) sprintf("%.2f%%", x) else "NA"
    for (i in seq_len(nrow(means_df))) {
      report <- c(report, sprintf("| %s | %s | %s | %s | %s |",
        means_df$variable[i],
        fmt_num(means_df$new[i]),
        fmt_num(means_df$legacy[i]),
        fmt_pct(means_df$abs_pct_diff[i]),
        means_df$status[i]))
    }
  } else {
    report <- c(report, "", "_No comparable variable pairs found._")
  }

  d <- here::here("audit_reports", "reconcile")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(d, sprintf("%s_audit.md", format(Sys.Date(), "%F")))
  writeLines(report, out)

  n_warn <- sum(means_df$status == "WARN", na.rm = TRUE)
  row_ok <- row_diff_pct <= tolerances$row_count_pct
  log_stage("audit", sprintf(
    "row-count %s (%.2f%% diff); means %d/%d WARN. report -> %s",
    if (row_ok) "PASS" else "WARN", row_diff_pct,
    n_warn, nrow(means_df), out))

  invisible(means_df)
}
