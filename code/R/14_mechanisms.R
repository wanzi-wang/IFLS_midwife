# =====================================================================
# 14_mechanisms.R — Mechanism / channel analysis.
#
# Question: if village-midwife exposure moves child health (as the
# headline suggests for girls and the prenatal subsample), through
# which CHANNELS does the effect travel? We test three mother-reported
# mediators extracted in stage 04b from the IFLS book-4 pregnancy
# history:
#
#   (1) pnc_any            — any prenatal check-up during pregnancy
#   (2) skilled_attendant  — birth attended by doctor/midwife/nurse
#                            (the canonical VM channel)
#   (3) ever_breastfed     — breastfeeding initiation
#
# Design:
#
#   Step 1 (first stage).    Regress each mediator M_ic on early VM
#                            exposure E_ic with the same TWFE + controls
#                            as the main specification. A non-zero
#                            first-stage estimate is necessary for the
#                            mediator to be a candidate channel.
#
#   Step 2 (reduced form).   Already reported in stage 08:
#                            Y_ic on E_ic.
#
#   Step 3 (conditional).    Regress Y_ic on E_ic AND M_ic. The
#                            attenuation of the E coefficient relative
#                            to Step 2, normalized by the first-stage,
#                            gives a suggestive decomposition: how
#                            much of the reduced-form effect is "net"
#                            of variation in the mediator.
#
# This is Baron–Kenny decomposition; it is suggestive, not causal.
# Standard errors on the indirect component are not reported because
# product-of-coefficients inference is ill-conditioned when either
# leg is imprecise; readers interested in formal mediation inference
# can consult \\textcite{imai2010}.
#
# Inputs :  data/intermediate/stage06/analysis_sample.rds
# Outputs:  paper/tables/mechanism_first_stage.tex
#           paper/tables/mechanism_decomposition.tex
#           data/output/regressions/stage14_mechanisms.rds
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

log_stage("stage14", "begin")

frame <- read_intermediate("stage06", "analysis_sample")

tables_dir <- here::here("paper", "tables")
reg_dir    <- here::here("data", "output", "regressions")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reg_dir,    recursive = TRUE, showWarnings = FALSE)

# NOTE: the mechanism analysis uses the wider preceding-wave sample
# rather than the non-migrant primary sample. Pregnancy-history
# coverage is thin enough on the non-migrant sample that the
# exposure variable becomes collinear with the community fixed
# effects after NA-dropping, so first-stage estimation is
# infeasible. The preceding-wave sample roughly doubles the
# population with mother-reported mediator data and restores
# within-community variation in exposure. This is explicitly
# flagged in paper/sections/mechanisms.tex as a separate sample
# from the reduced-form analysis.
df <- frame |>
  filter(preceding_sample == 1) |>
  mutate(
    sex_f = factor(sex, levels = c(1, 3), labels = c("male", "female"))
  )

mediators <- c("pnc_any", "skilled_attendant", "ever_breastfed")
outcomes  <- c("health_index", "cognition_index",
               "depression_index", "bigfive_index")

# ---------------------------------------------------------------------
# 1. First stage: VM exposure → each mediator.
# ---------------------------------------------------------------------
first_stage <- lapply(mediators, function(m) {
  fit <- tryCatch(
    feols(as.formula(sprintf(paste0(
      "%s ~ exposure_early_sar_preceding + sex_f + mother_edu_years + ",
      "mother_age_birth | commid_birth_preceding + birth_year"), m)),
      data = df, cluster = ~ commid_birth_preceding,
      notes = FALSE, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) ||
      !"exposure_early_sar_preceding" %in% names(coef(fit))) {
    return(tibble(mediator = m, n = NA_integer_,
                  beta = NA_real_, se = NA_real_, p = NA_real_,
                  mean = mean(df[[m]], na.rm = TRUE)))
  }
  co <- coef(fit)["exposure_early_sar_preceding"]
  se <- sqrt(vcov(fit)["exposure_early_sar_preceding",
                       "exposure_early_sar_preceding"])
  tibble(
    mediator = m,
    n     = fit$nobs,
    beta  = unname(co),
    se    = unname(se),
    p     = 2 * pnorm(-abs(unname(co / se))),
    mean  = mean(df[[m]], na.rm = TRUE)
  )
}) |> bind_rows()

first_stage$stars <- ifelse(first_stage$p < 0.01, "$^{***}$",
                     ifelse(first_stage$p < 0.05,  "$^{**}$",
                     ifelse(first_stage$p < 0.10,   "$^{*}$", "")))

med_labels <- c(
  pnc_any           = "Any prenatal check-up",
  skilled_attendant = "Skilled birth attendant",
  ever_breastfed    = "Ever breastfed"
)
first_stage$mediator_label <- med_labels[first_stage$mediator]

log_stage("stage14", sprintf(
  "first-stage estimates: %s",
  paste(sprintf("%s beta=%.3f se=%.3f p=%.3f N=%d",
                first_stage$mediator,
                first_stage$beta, first_stage$se,
                first_stage$p, first_stage$n),
        collapse = "; ")))

writeLines(c(
  "% Auto-generated by code/R/14_mechanisms.R",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Mediator & Sample mean & $\\hat{\\beta}$ & SE & $N$ \\\\",
  "\\midrule",
  apply(first_stage, 1, function(r) sprintf(
    "%s & %.3f & %.3f%s & (%.3f) & %s \\\\",
    as.character(r["mediator_label"]),
    as.numeric(r["mean"]),
    as.numeric(r["beta"]),
    as.character(r["stars"]),
    as.numeric(r["se"]),
    format(as.integer(r["n"]), big.mark = ",")
  )),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "mechanism_first_stage.tex"))

log_stage("stage14", "first-stage table written.")

# ---------------------------------------------------------------------
# 2. Baron–Kenny decomposition on the canonical channel:
#    skilled_attendant → health_index.
#    We also run the full 4-outcome × 3-mediator grid for the
#    appendix.
# ---------------------------------------------------------------------
run_cond <- function(y, m, data) {
  # Total (reduced form).
  f_total <- as.formula(sprintf(paste0(
    "%s ~ exposure_early_sar_preceding + sex_f + mother_edu_years + ",
    "mother_age_birth | commid_birth_preceding + birth_year"), y))
  # Conditional (adds mediator as regressor).
  f_cond <- as.formula(sprintf(paste0(
    "%s ~ exposure_early_sar_preceding + %s + sex_f + mother_edu_years + ",
    "mother_age_birth | commid_birth_preceding + birth_year"), y, m))
  fit_t <- feols(f_total, data = data,
                 cluster = ~ commid_birth_preceding,
                 notes = FALSE, warn = FALSE)
  fit_c <- tryCatch(feols(f_cond, data = data,
                          cluster = ~ commid_birth_preceding,
                          notes = FALSE, warn = FALSE),
                    error = function(e) NULL)
  if (is.null(fit_c)) return(NULL)
  list(total = fit_t, cond = fit_c)
}

decomp_rows <- list()
for (y in outcomes) {
  for (m in mediators) {
    # Restrict to rows with non-missing mediator so total and
    # conditional are estimated on the same sample.
    sub <- df |> filter(!is.na(.data[[m]]))
    res <- run_cond(y, m, sub)
    if (is.null(res)) next
    if (!"exposure_early_sar_preceding" %in% names(coef(res$total)) ||
        !"exposure_early_sar_preceding" %in% names(coef(res$cond))) next
    co_t <- coef(res$total)["exposure_early_sar_preceding"]
    se_t <- sqrt(vcov(res$total)["exposure_early_sar_preceding",
                                 "exposure_early_sar_preceding"])
    co_c <- coef(res$cond)["exposure_early_sar_preceding"]
    se_c <- sqrt(vcov(res$cond)["exposure_early_sar_preceding",
                                "exposure_early_sar_preceding"])
    # Share mediated = (total - direct) / total.
    share <- (co_t - co_c) / co_t
    decomp_rows[[paste0(y, "_", m)]] <- tibble(
      outcome  = y,
      mediator = m,
      beta_total = unname(co_t), se_total = unname(se_t),
      beta_cond  = unname(co_c), se_cond  = unname(se_c),
      share_mediated = unname(share),
      n = res$cond$nobs
    )
  }
}
decomp <- bind_rows(decomp_rows)

# Suppress noisy % shares when the total effect is very close to zero
# (the ratio is unstable). We report "n/a" when |beta_total| < 1e-3
# and flag it in the note.
decomp <- decomp |>
  mutate(
    share_mediated_disp = ifelse(
      abs(beta_total) < 1e-3, NA_real_, share_mediated * 100)
  )

decomp$mediator_label <- med_labels[decomp$mediator]
decomp$outcome_label  <- pretty_label(decomp$outcome)

fmt_pct <- function(x) ifelse(is.na(x), "---",
                              sprintf("%+6.1f\\%%", x))

writeLines(c(
  "% Auto-generated by code/R/14_mechanisms.R",
  "\\begin{tabular}{llcccc}",
  "\\toprule",
  paste("Outcome & Mediator & $\\hat{\\beta}_{\\text{total}}$ &",
        "$\\hat{\\beta}_{\\text{conditional}}$ &",
        "Share mediated & $N$ \\\\"),
  "\\midrule",
  unlist(lapply(unique(decomp$outcome), function(y) {
    sub <- decomp[decomp$outcome == y, , drop = FALSE]
    lines <- apply(sub, 1, function(r) sprintf(
      "%s & %s & \\makecell{%.3f\\\\(%.3f)} & \\makecell{%.3f\\\\(%.3f)} & %s & %s \\\\",
      ifelse(r["mediator"] == sub$mediator[1], as.character(r["outcome_label"]), ""),
      as.character(r["mediator_label"]),
      as.numeric(r["beta_total"]), as.numeric(r["se_total"]),
      as.numeric(r["beta_cond"]),  as.numeric(r["se_cond"]),
      fmt_pct(as.numeric(r["share_mediated_disp"])),
      format(as.integer(r["n"]), big.mark = ",")
    ))
    c(lines, "\\midrule")
  })) |> head(-1),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "mechanism_decomposition.tex"))

log_stage("stage14", "decomposition table written.")

# ---------------------------------------------------------------------
# 3. Persist.
# ---------------------------------------------------------------------
saveRDS(list(
  first_stage   = first_stage,
  decomposition = decomp
), file.path(reg_dir, "stage14_mechanisms.rds"))

log_stage("stage14", "done")
