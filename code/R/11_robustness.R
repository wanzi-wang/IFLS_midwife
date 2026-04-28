# =====================================================================
# 11_robustness.R — CS-based robustness suite.
#
# Each cell reports the Callaway-Sant'Anna dynamic ATT averaged over
# event-times 0--3 under a perturbation of the baseline CS specification.
# Aligning the suite with the headline (CS) estimator avoids the
# previous TWFE-vs-CS sign confusion: TWFE is shown elsewhere as the
# Goodman-Bacon-biased baseline, and is not the reference here.
#
# Columns per primary outcome:
#   (R1) IPW-weighted CS
#   (R2) PKK rollout source (replaces SAR rollout)
#   (R3) Preceding-wave birthplace variant
#   (R4) Mother birthplace variant
#
# Facility-control and province-trend perturbations were dropped: in the
# CS framework the cohort-dependent facility indicators are collinear
# with the group-by-time structure (singular design matrix), and CS does
# not natively support province linear trends.
#
# Separate diagnostics written below:
#   - Placebo (earlier rollout): shift G five years earlier and rerun
#     CS. The placebo cohort is then aged ~5--8 at the actual midwife
#     arrival, outside the in-utero-to-three critical window.
#   - Placebo (older cohort): redefine the exposure indicator as
#     "child aged 11-15 at midwife arrival" and refit TWFE. Children
#     aged 11-15 at rollout are too old to benefit biologically from
#     early-life midwife exposure, so the coefficient should be null.
#   - Child Raven (cog_prop_child) — cognition-only robustness.
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(fixest)
  library(did)
})

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))
source(here::here("code", "R", "lib", "labels.R"))

log_stage("stage11", "begin")

frame <- read_intermediate("stage06", "analysis_sample")
tables_dir <- here::here("paper", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

outcomes <- c("health_index", "cognition_index",
              "depression_index", "bigfive_index")

stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) "$^{***}$"
  else if (p < 0.05) "$^{**}$"
  else if (p < 0.10) "$^{*}$"
  else ""
}

fmt <- function(est, se) {
  if (is.na(est) || is.na(se) || se == 0) return("---")
  p <- 2 * stats::pnorm(-abs(est / se))
  sprintf("%.3f%s\n(%.3f)", est, stars(p), se)
}

# ---------------------------------------------------------------------
# Build perturbation-specific frames.
# ---------------------------------------------------------------------

# R1, R2, R3: primary sample (non-migrant, legacy birthplace).
df <- frame |>
  filter(primary_sample == 1) |>
  mutate(
    G_sar = if_else(!is.na(start_sar_legacy),
                    as.integer(start_sar_legacy) - 1L, 0L),
    G_pkk = if_else(!is.na(start_pkk_legacy),
                    as.integer(start_pkk_legacy) - 1L, 0L),
    commid_num = as.integer(factor(commid_birth_legacy))
  )

# R4: preceding-wave birthplace.
df_pre <- frame |>
  filter(preceding_sample == 1) |>
  mutate(
    G_sar = if_else(!is.na(start_sar_preceding),
                    as.integer(start_sar_preceding) - 1L, 0L),
    commid_num = as.integer(factor(commid_birth_preceding))
  )

# R5: mother birthplace.
df_mom <- frame |>
  filter(mother_sample == 1) |>
  mutate(
    G_sar = if_else(!is.na(start_sar_mother),
                    as.integer(start_sar_mother) - 1L, 0L),
    commid_num = as.integer(factor(commid_birth_mother))
  )

# ---------------------------------------------------------------------
# CS dynamic ATT (event-times 0--3) under a perturbation.
# ---------------------------------------------------------------------
fit_cs <- function(y, data, gname = "G_sar", idname = "commid_num",
                   weightsname = NULL, xformla = NULL,
                   min_e = 0, max_e = 3) {
  d <- data |>
    filter(!is.na(.data[[y]]),
           !is.na(birth_year),
           !is.na(.data[[idname]]),
           !is.na(.data[[gname]]))
  if (!is.null(weightsname))
    d <- d |> filter(!is.na(.data[[weightsname]]))
  if (!is.null(xformla)) {
    xv <- all.vars(xformla)
    for (v in xv) d <- d |> filter(!is.na(.data[[v]]))
  }
  d <- as.data.frame(d)
  if (nrow(d) == 0)
    return(list(est = NA_real_, se = NA_real_,
                n = 0L, nclust = 0L))
  tryCatch({
    att <- att_gt(
      yname = y, tname = "birth_year", idname = idname,
      gname = gname, data = d, control_group = "nevertreated",
      bstrap = FALSE, cband = FALSE,
      allow_unbalanced_panel = TRUE, panel = FALSE,
      weightsname = weightsname,
      xformla = xformla
    )
    dyn <- aggte(att, type = "dynamic",
                 min_e = min_e, max_e = max_e,
                 na.rm = TRUE)
    list(est = unname(dyn$overall.att),
         se  = unname(dyn$overall.se),
         n = nrow(d),
         nclust = length(unique(d[[idname]])))
  }, error = function(e) {
    message(sprintf("CS failed for %s: %s", y, e$message))
    list(est = NA_real_, se = NA_real_,
         n = NA_integer_, nclust = NA_integer_)
  })
}

# ---------------------------------------------------------------------
# Build robustness rows per outcome.
# ---------------------------------------------------------------------
R <- list()

for (y in outcomes) {
  log_stage("stage11", sprintf("CS robustness for %s", y))

  r1 <- fit_cs(y, df, gname = "G_sar", weightsname = "ipw")
  r2 <- fit_cs(y, df, gname = "G_pkk")
  r3 <- fit_cs(y, df_pre, gname = "G_sar")
  r4 <- fit_cs(y, df_mom, gname = "G_sar")

  R[[y]] <- tibble(
    outcome = y,
    `R1 IPW`        = fmt(r1$est, r1$se),
    `R2 PKK`        = fmt(r2$est, r2$se),
    `R3 preceding`  = fmt(r3$est, r3$se),
    `R4 mother`     = fmt(r4$est, r4$se),
    n1 = r1$n, n2 = r2$n, n3 = r3$n, n4 = r4$n,
    c1 = r1$nclust, c2 = r2$nclust, c3 = r3$nclust, c4 = r4$nclust
  )
}

rob_tbl <- bind_rows(R)

n_max <- function(col) suppressWarnings(max(rob_tbl[[col]], na.rm = TRUE))
ns <- vapply(paste0("n", 1:4), n_max, numeric(1))
cs <- vapply(paste0("c", 1:4), n_max, numeric(1))
fmt_big <- function(x) if (!is.finite(x)) "---" else format(as.integer(x),
                                                            big.mark = ",")

coef_lines <- apply(
  rob_tbl[, c("outcome", "R1 IPW", "R2 PKK",
              "R3 preceding", "R4 mother")],
  1, function(r) {
    cells <- vapply(r[-1], function(cell) {
      if (grepl("\n", cell))
        paste0("\\makecell{", gsub("\n", "\\\\\\\\", cell), "}")
      else cell
    }, character(1))
    sprintf("%s & %s \\\\", pretty_label(r["outcome"]),
            paste(cells, collapse = " & "))
  })

writeLines(c(
  "% Auto-generated by code/R/11_robustness.R",
  "\\begin{tabular}{l*{4}{c}}",
  "\\toprule",
  paste("Outcome & (R1) IPW & (R2) PKK rollout &",
        "(R3) Preceding bp & (R4) Mother bp \\\\"),
  "\\midrule",
  coef_lines,
  "\\midrule",
  sprintf("Clusters (communities) & %s & %s & %s & %s \\\\",
          fmt_big(cs[1]), fmt_big(cs[2]), fmt_big(cs[3]), fmt_big(cs[4])),
  sprintf("Observations & %s & %s & %s & %s \\\\",
          fmt_big(ns[1]), fmt_big(ns[2]), fmt_big(ns[3]), fmt_big(ns[4])),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "robustness.tex"))

# ---------------------------------------------------------------------
# Placebo: earlier rollout. Shift G five years earlier; re-fit CS.
# The placebo cohort is then aged ~5-8 at the actual midwife arrival,
# outside the in-utero-to-three critical window. Effect should be null.
# ---------------------------------------------------------------------
log_stage("stage11", "earlier-rollout (CS) and older-cohort (TWFE) placebos")
df_early <- df |>
  mutate(G_early = if_else(G_sar > 0, G_sar - 5L, 0L))

# Older-cohort placebo: redefine exposure as "aged 11-15 at midwife
# arrival" and refit TWFE. Sample includes never-treated villages so the
# control group has exposure_old = 0 mechanically.
df_old <- frame |>
  filter(primary_sample == 1) |>
  mutate(age_at_rollout = if_else(!is.na(start_sar_legacy),
                                   as.integer(start_sar_legacy) - birth_year,
                                   NA_integer_),
         exposure_old = as.integer(!is.na(age_at_rollout) &
                                     age_at_rollout >= 11 &
                                     age_at_rollout <= 15))

placebo_alt <- list()
for (y in outcomes) {
  pe <- fit_cs(y, df_early, gname = "G_early")
  m_old <- tryCatch(feols(as.formula(sprintf(paste0(
    "%s ~ exposure_old | commid_birth_legacy + birth_year + source_wave"),
    y)), data = df_old, cluster = ~ commid_birth_legacy),
    error = function(e) NULL)
  if (!is.null(m_old)) {
    po_est <- unname(coef(m_old)["exposure_old"])
    po_se  <- unname(sqrt(vcov(m_old)["exposure_old", "exposure_old"]))
  } else {
    po_est <- NA_real_; po_se <- NA_real_
  }
  placebo_alt[[y]] <- tibble(
    outcome = y,
    `Earlier rollout` = fmt(pe$est, pe$se),
    `Older cohort`    = fmt(po_est, po_se)
  )
}
placebo_alt_tbl <- bind_rows(placebo_alt)

placebo_alt_lines <- apply(
  placebo_alt_tbl[, c("outcome", "Earlier rollout", "Older cohort")],
  1, function(r) {
    cells <- vapply(r[-1], function(cell) {
      if (grepl("\n", cell))
        paste0("\\makecell{", gsub("\n", "\\\\\\\\", cell), "}")
      else cell
    }, character(1))
    sprintf("%s & %s \\\\", pretty_label(r["outcome"]),
            paste(cells, collapse = " & "))
  })

writeLines(c(
  "% Auto-generated by code/R/11_robustness.R",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Outcome & (P1) Earlier rollout & (P2) Older cohort \\\\",
  "\\midrule",
  placebo_alt_lines,
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "placebo_alt_cohorts.tex"))

# ---------------------------------------------------------------------
# Child Raven (cog_prop_child) — cognition-only robustness.
# ---------------------------------------------------------------------
df_twfe <- frame |>
  filter(primary_sample == 1) |>
  mutate(sex_f = factor(sex, levels = c(1, 3),
                        labels = c("male", "female")))

m_child <- tryCatch(feols(
  cog_prop_child ~ exposure_early_sar_legacy |
    commid_birth_legacy + birth_year + source_wave,
  data = df_twfe, cluster = ~ commid_birth_legacy),
  error = function(e) NULL)

if (!is.null(m_child)) {
  est <- coef(m_child)["exposure_early_sar_legacy"]
  se  <- sqrt(vcov(m_child)["exposure_early_sar_legacy",
                            "exposure_early_sar_legacy"])
  n_child <- fixest::fitstat(m_child, "n")$n
  tval <- est / se
  pval <- 2 * stats::pnorm(-abs(tval))
  star <- if (pval < 0.01) "$^{***}$" else if (pval < 0.05) "$^{**}$" else if (pval < 0.10) "$^{*}$" else ""
  writeLines(c(
    "% Auto-generated by code/R/11_robustness.R",
    "\\begin{tabular}{lccc}",
    "\\toprule",
    "Outcome & Estimate & SE & N \\\\",
    "\\midrule",
    sprintf("%s & %.3f%s & (%.3f) & %s \\\\",
            pretty_label("cog_prop_child"), est, star, se,
            format(as.integer(n_child), big.mark = ",")),
    "\\bottomrule",
    "\\end{tabular}"
  ), file.path(tables_dir, "robustness_child_raven.tex"))
  log_stage("stage11", sprintf(
    "child Raven: coef=%.4f (SE %.4f); N=%d", est, se, n_child))
}

log_stage("stage11", "done")
