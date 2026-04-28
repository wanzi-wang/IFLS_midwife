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
#   First-stage.   Regress each mediator M_ic on early VM exposure
#                  E_ic with the same TWFE + controls as the main
#                  specification. A non-zero first-stage estimate is
#                  necessary for the mediator to be a candidate channel.
#
# A Baron-Kenny channel-share decomposition was dropped 2026-04-28:
# on the pregnancy-history sub-sample the reduced-form effect is small
# and imprecise, so any share-mediated quantity inherits the imprecision
# of the numerator and does not support a credible channel-share
# inference.
#
# Inputs :  data/intermediate/stage06/analysis_sample.rds
# Outputs:  paper/tables/mechanism_first_stage.tex
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
  filter(birth_year >= 1989) |>   # match Ahsan: post-program cohorts only
  mutate(
    sex_f = factor(sex, levels = c(1, 3), labels = c("male", "female"))
  )

# Mediators measured at or shortly after birth.  The cleanest first-
# stage exposure for these is "in-utero exposure": the program was
# present in the village at the time the pregnancy was occurring.
# Children whose program arrived after birth (within the broader
# "exposure_early" window) cannot have their birth-attendant or birth-
# weight outcomes moved by the program, so using the wider exposure
# variable would attenuate the first-stage mechanically.
exposure_var <- "inutero_sar_preceding"
mediators    <- c("pnc_any", "skilled_attendant",
                  "attendant_village_midwife",
                  "exclusive_bf_6m", "birth_weight")
outcomes     <- c("health_index", "cognition_index",
                  "depression_index", "bigfive_index")

# ---------------------------------------------------------------------
# 1. First stage: VM in-utero exposure -> each mediator.
# ---------------------------------------------------------------------
first_stage <- lapply(mediators, function(m) {
  fit <- tryCatch(
    feols(as.formula(sprintf(paste0(
      "%s ~ %s + sex_f + mother_edu_years + ",
      "mother_age_birth | commid_birth_preceding + birth_year + source_wave"),
      m, exposure_var)),
      data = df, cluster = ~ commid_birth_preceding,
      notes = FALSE, warn = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || !exposure_var %in% names(coef(fit))) {
    return(tibble(mediator = m, n = NA_integer_,
                  beta = NA_real_, se = NA_real_, p = NA_real_,
                  mean = mean(df[[m]], na.rm = TRUE)))
  }
  co <- coef(fit)[exposure_var]
  se <- sqrt(vcov(fit)[exposure_var, exposure_var])
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
  pnc_any                    = "Any prenatal check-up",
  skilled_attendant          = "Any skilled attendant at birth",
  attendant_village_midwife  = "Village midwife at birth",
  exclusive_bf_6m            = "Exclusive breastfeeding $\\geq 6$ mo.",
  birth_weight               = "Birth weight (kg)"
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
# 2. Persist first-stage results only.  The Baron-Kenny decomposition
# was dropped 2026-04-28: on the pregnancy-history sub-sample the
# reduced-form effect is small and imprecise, so any share-mediated
# quantity inherits the imprecision of the numerator and does not
# support a credible channel-share inference.  See discussion in
# paper/sections/mechanisms.tex.
# ---------------------------------------------------------------------
saveRDS(list(
  first_stage = first_stage
), file.path(reg_dir, "stage14_mechanisms.rds"))

log_stage("stage14", "done")
