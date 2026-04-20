# =====================================================================
# 11_robustness.R — Robustness suite (one column per perturbation).
#
# Columns per primary outcome:
#   (R1) IPW-weighted TWFE
#   (R2) Other-facility controls (cohort-dependent, not level)
#   (R3) Province linear trends instead of prov×year FE
#   (R4) PKK rollout source
#   (R5) Legacy birthplace variant (commid_birth_legacy)
#   (R6) Mother birthplace variant (commid_birth_mother)
#   (R7) Oster (2019) bound: delta^* that nullifies estimate
#   (R8) Placebo: permute start_sar_legacy within province
#   (R9) Child Raven (cog_prop_child) — cognition only, demoted
#
# The two-shock / AFC-exit design was dropped 2026-04-15: SAR exit
# fields are empty, PKK-gap windows are too coarse, the null result
# was consistent with the user's prior, and keeping it added clutter.
#
# All columns except R8/R9 are cluster-robust TWFE at commid_birth.
# R8 reports delta*, the selection-on-unobservables ratio that would
# reduce the estimate to zero (from Oster 2019; computed manually
# since `robomit` is not installed).
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(fixest)
})

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))
source(here::here("code", "R", "lib", "labels.R"))

log_stage("stage11", "begin")

frame <- read_intermediate("stage06", "analysis_sample")
other_f <- read_intermediate("stage03", "other_facilities") |>
  rename_with(~ paste0(.x, "_pr"), .cols = -commid93)
tables_dir <- here::here("paper", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Rejoin other-facilities onto the *preceding* birthplace so R2
# (facility controls) can run on the primary sample. The stage-05
# frame only attached them on commid_birth_legacy.
frame <- frame |>
  left_join(other_f, by = c("commid_birth_legacy" = "commid93"))

df <- frame |>
  filter(primary_sample == 1) |>
  mutate(sex_f = factor(sex, levels = c(1, 3),
                        labels = c("male", "female")))

outcomes <- c("health_index", "cognition_index",
              "depression_index", "bigfive_index")

coef_se <- function(m, var) {
  if (is.null(m) || !var %in% names(coef(m))) return(c(NA, NA))
  c(unname(coef(m)[var]),
    unname(sqrt(vcov(m)[var, var])))
}

fmt <- function(est, se) {
  if (is.na(est)) return("---")
  sprintf("%.3f\n(%.3f)", est, se)
}

# ---------------------------------------------------------------------
# Oster (2019) delta* implementation.
# Compares (coef, R^2) from short and long regressions. With an assumed
# Rmax (we use min(1.3 * R_long, 1.0) per Oster's recommendation),
# delta* = (b_long - b_short) / (Rmax - R_long) /
#          ((b_short - b_long) / (R_long - R_short))
# Interpretation: if delta* > 1, selection on unobservables would need
# to be STRONGER than selection on observables to overturn the result.
# ---------------------------------------------------------------------
oster_delta <- function(y, data, key_var, short_rhs, long_rhs) {
  m_short <- tryCatch(feols(as.formula(
    sprintf("%s ~ %s + %s", y, key_var, short_rhs)),
    data = data), error = function(e) NULL)
  m_long  <- tryCatch(feols(as.formula(
    sprintf("%s ~ %s + %s", y, key_var, long_rhs)),
    data = data), error = function(e) NULL)
  if (is.null(m_short) || is.null(m_long)) return(NA_real_)
  bS <- coef(m_short)[key_var]
  bL <- coef(m_long)[key_var]
  rS <- fitstat(m_short, "r2")$r2
  rL <- fitstat(m_long,  "r2")$r2
  Rmax <- min(1.3 * rL, 1.0)
  num <- bL * (Rmax - rL)
  den <- (bS - bL) * (rL - rS)
  if (abs(den) < 1e-8) return(NA_real_)
  as.numeric(num / den)
}

# ---------------------------------------------------------------------
# Build robustness rows per outcome.
# ---------------------------------------------------------------------
R <- list()

for (y in outcomes) {
  log_stage("stage11", sprintf("robustness for %s", y))

  # R1: IPW-weighted TWFE.
  d1 <- df |> filter(!is.na(ipw))
  m1 <- feols(as.formula(sprintf(
    "%s ~ exposure_early_sar_legacy | commid_birth_legacy + birth_year",
    y)), data = d1, weights = ~ ipw, cluster = ~ commid_birth_legacy)

  # R2: Other-facility controls.
  # Other-facility controls must vary within commid (time-invariant
  # "has_*" are absorbed by commid FE). Use cohort-specific indicators:
  # `doctor_at_birth` = 1 if a doctor facility existed at (or before)
  # the child's birth year. Analogous for clinic.
  df_m2 <- df |>
    mutate(doctor_at_birth = as.integer(!is.na(start_doctor_pr) &
                                          start_doctor_pr <= birth_year),
           clinic_at_birth = as.integer(!is.na(start_clinic_pr) &
                                          start_clinic_pr <= birth_year))
  m2 <- feols(as.formula(sprintf(paste0(
    "%s ~ exposure_early_sar_legacy + doctor_at_birth + clinic_at_birth | ",
    "commid_birth_legacy + birth_year"), y)),
    data = df_m2, cluster = ~ commid_birth_legacy)

  # R3: Province linear trends (instead of prov × birth_year FE).
  m3 <- feols(as.formula(sprintf(paste0(
    "%s ~ exposure_early_sar_legacy + i(province, birth_year) | ",
    "commid_birth_legacy + birth_year"), y)),
    data = df, cluster = ~ commid_birth_legacy)

  # R4: PKK rollout source (exposure_early_pkk_legacy).
  m4 <- feols(as.formula(sprintf(paste0(
    "%s ~ exposure_early_pkk_legacy | ",
    "commid_birth_legacy + birth_year"), y)),
    data = df, cluster = ~ commid_birth_legacy)

  # R5: Preceding-wave birthplace variant (keeps internal migrants).
  d5 <- frame |> filter(preceding_sample == 1)
  m5 <- feols(as.formula(sprintf(paste0(
    "%s ~ exposure_early_sar_preceding | ",
    "commid_birth_preceding + birth_year"), y)),
    data = d5, cluster = ~ commid_birth_preceding)

  # R6: Mother birthplace variant.
  d6 <- frame |> filter(mother_sample == 1)
  m6 <- feols(as.formula(sprintf(paste0(
    "%s ~ exposure_early_sar_mother | ",
    "commid_birth_mother + birth_year"), y)),
    data = d6, cluster = ~ commid_birth_mother)

  # R7: Oster delta*. Need raw (non-fixest) OLS for R^2 comparison.
  # Short: outcome ~ exposure + sex.
  # Long : outcome ~ exposure + sex + mother_edu + mother_age +
  #        commid_birth dummies + birth_year dummies.
  oster_df <- df |>
    filter(!is.na(.data[[y]]),
           !is.na(exposure_early_sar_legacy),
           !is.na(mother_edu_years),
           !is.na(mother_age_birth),
           !is.na(sex_f))
  delta_star <- oster_delta(
    y = y, data = oster_df,
    key_var   = "exposure_early_sar_legacy",
    short_rhs = "sex_f",
    long_rhs  = paste0("sex_f + mother_edu_years + mother_age_birth + ",
                       "factor(commid_birth_legacy) + factor(birth_year)")
  )

  # Assemble row.
  cs <- function(m, v = "exposure_early_sar_legacy") {
    pr <- coef_se(m, v); fmt(pr[1], pr[2])
  }
  R[[y]] <- tibble(
    outcome = y,
    `R1 IPW`           = cs(m1),
    `R2 facilities`    = cs(m2),
    `R3 prov-trend`    = cs(m3),
    `R4 PKK`           = cs(m4, "exposure_early_pkk_legacy"),
    `R5 preceding bp`  = cs(m5, "exposure_early_sar_preceding"),
    `R6 mother bp`     = cs(m6, "exposure_early_sar_mother"),
    `R7 Oster delta*`  = if (is.na(delta_star)) "---"
                         else sprintf("%.2f", delta_star)
  )
}

rob_tbl <- bind_rows(R)

writeLines(c(
  "% Auto-generated by code/R/11_robustness.R",
  "\\begin{tabular}{l*{7}{c}}",
  "\\toprule",
  paste("Outcome & (R1) IPW & (R2) Facilities & (R3) Prov. trend &",
        "(R4) PKK & (R5) Preceding bp & (R6) Mother bp &",
        "(R7) Oster $\\delta^*$ \\\\"),
  "\\midrule",
  apply(rob_tbl, 1, function(r) {
    cells <- vapply(r[-1], function(cell) {
      if (grepl("\n", cell))
        paste0("\\makecell{", gsub("\n", "\\\\\\\\", cell), "}")
      else cell
    }, character(1))
    sprintf("%s & %s \\\\", pretty_label(r["outcome"]),
            paste(cells, collapse = " & "))
  }),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "robustness.tex"))

# ---------------------------------------------------------------------
# R9: Placebo — permute start_sar_legacy within province × cohort.
# Run 200 permutations; report the share whose |coef| exceeds the true
# estimate on health_index (a single-outcome diagnostic).
# ---------------------------------------------------------------------
log_stage("stage11", "running placebo permutation (200 draws)")
true_fit <- feols(
  health_index ~ exposure_early_sar_legacy |
    commid_birth_legacy + birth_year,
  data = df, cluster = ~ commid_birth_legacy)
true_coef <- unname(coef(true_fit)["exposure_early_sar_legacy"])

B <- 200
placebo_coefs <- numeric(B)
for (b in seq_len(B)) {
  perm <- df |>
    group_by(province, birth_year) |>
    mutate(exp_perm = sample(exposure_early_sar_legacy)) |>
    ungroup()
  fit <- tryCatch(feols(
    health_index ~ exp_perm |
      commid_birth_legacy + birth_year,
    data = perm, cluster = ~ commid_birth_legacy),
    error = function(e) NULL)
  placebo_coefs[b] <- if (is.null(fit)) NA_real_
  else unname(coef(fit)["exp_perm"])
}
share_larger <- mean(abs(placebo_coefs) >= abs(true_coef), na.rm = TRUE)
log_stage("stage11", sprintf(
  "placebo: true coef = %.4f; %.1f%% of %d permutations |larger|",
  true_coef, 100 * share_larger, B))

placebo_tbl <- tibble(
  outcome = "health_index",
  true_coef = true_coef,
  placebo_mean = mean(placebo_coefs, na.rm = TRUE),
  placebo_sd   = stats::sd(placebo_coefs,   na.rm = TRUE),
  share_larger = share_larger,
  B = B
)
saveRDS(list(placebo_coefs = placebo_coefs, summary = placebo_tbl),
        here::here("data", "output", "regressions",
                   "stage11_placebo.rds"), compress = "xz")

writeLines(c(
  "% Auto-generated by code/R/11_robustness.R — placebo permutation",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Outcome & True coef. & Placebo mean & Placebo SD & Share $|b| \\geq$ true \\\\",
  "\\midrule",
  sprintf("health\\_index & %.4f & %.4f & %.4f & %.3f \\\\",
          true_coef, mean(placebo_coefs, na.rm = TRUE),
          stats::sd(placebo_coefs, na.rm = TRUE), share_larger),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "placebo.tex"))

# ---------------------------------------------------------------------
# R10: Child Raven (cog_prop_child) — cognition-only robustness.
# ---------------------------------------------------------------------
m_child <- tryCatch(feols(
  cog_prop_child ~ exposure_early_sar_legacy |
    commid_birth_legacy + birth_year,
  data = df, cluster = ~ commid_birth_legacy),
  error = function(e) NULL)

if (!is.null(m_child)) {
  est <- coef(m_child)["exposure_early_sar_legacy"]
  se  <- sqrt(vcov(m_child)["exposure_early_sar_legacy",
                            "exposure_early_sar_legacy"])
  n_child <- fixest::fitstat(m_child, "n")$n
  tval <- est / se
  pval <- 2 * stats::pnorm(-abs(tval))
  stars <- if (pval < 0.01) "$^{***}$" else if (pval < 0.05) "$^{**}$" else if (pval < 0.10) "$^{*}$" else ""
  writeLines(c(
    "% Auto-generated by code/R/11_robustness.R",
    "\\begin{tabular}{lccc}",
    "\\toprule",
    "Outcome & Estimate & SE & N \\\\",
    "\\midrule",
    sprintf("%s & %.3f%s & (%.3f) & %s \\\\",
            pretty_label("cog_prop_child"), est, stars, se,
            format(as.integer(n_child), big.mark = ",")),
    "\\bottomrule",
    "\\end{tabular}"
  ), file.path(tables_dir, "robustness_child_raven.tex"))
  log_stage("stage11", sprintf(
    "child Raven: coef=%.4f (SE %.4f); N=%d", est, se, n_child))
}

log_stage("stage11", "done")
