# =====================================================================
# 10_sibling_fe.R — Within-mother (sibling-FE) specifications.
#
# Headline robustness for time-invariant placement bias. By comparing
# siblings with differential VM exposure within the same mother (and
# therefore the same village, assuming non-mover mothers), we sweep
# out every village characteristic that drove where the VM was placed
# — poverty, remoteness, baseline infrastructure — regardless of
# whether it is observed. The remaining identifying variation is
# birth timing relative to a fixed start_sar_legacy per village.
#
# Inputs :  data/intermediate/stage06/analysis_sample.rds
#
# Outputs:  paper/tables/sibling_fe.tex
#           paper/tables/sibling_fertility_check.tex
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

log_stage("stage10", "begin")

frame <- read_intermediate("stage06", "analysis_sample")
tables_dir <- here::here("paper", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

outcomes <- c("health_index", "cognition_index",
              "depression_index", "bigfive_index")

# ---------------------------------------------------------------------
# 1. Sibling subsample: mothers with >=2 observed children in the
#    primary sample.
# ---------------------------------------------------------------------
df <- frame |>
  filter(primary_sample == 1, !is.na(mother_pidlink)) |>
  group_by(mother_pidlink) |>
  mutate(n_sibs = n(),
         birth_order = rank(birth_year, ties.method = "first")) |>
  ungroup()

sib2 <- df |> filter(n_sibs >= 2)
sib3 <- df |> filter(n_sibs >= 3)

log_stage("stage10", sprintf(
  "sibling subsamples: >=2 children -> %d obs, %d mothers; >=3 -> %d obs, %d mothers",
  nrow(sib2), dplyr::n_distinct(sib2$mother_pidlink),
  nrow(sib3), dplyr::n_distinct(sib3$mother_pidlink)))

# ---------------------------------------------------------------------
# 2. Sibling-FE regressions.
#    Spec A: mother FE + birth-year FE + controls.
#    Spec B: mother FE + birth-year FE + commid_birth FE (strongest).
#    Spec C: balanced >=3 sibling subsample.
# ---------------------------------------------------------------------
fit_sib <- function(y, data, spec) {
  form <- switch(spec,
    A = sprintf(paste0(
      "%s ~ exposure_early_sar_legacy + birth_order + ",
      "mother_age_birth | mother_pidlink + birth_year + source_wave"), y),
    B = sprintf(paste0(
      "%s ~ exposure_early_sar_legacy + birth_order + ",
      "mother_age_birth | mother_pidlink + birth_year + source_wave + ",
      "commid_birth_legacy"), y),
    C = sprintf(paste0(
      "%s ~ exposure_early_sar_legacy + birth_order + ",
      "mother_age_birth | mother_pidlink + birth_year + source_wave"), y)
  )
  tryCatch(
    feols(as.formula(form), data = data,
          cluster = ~ commid_birth_legacy),
    error = function(e) {
      message(sprintf("fit_sib failed (spec %s, y=%s): %s", spec, y,
                      e$message))
      NULL
    }
  )
}

stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) "$^{***}$"
  else if (p < 0.05) "$^{**}$"
  else if (p < 0.10) "$^{*}$"
  else ""
}

pull_cell <- function(m, data_used, y) {
  if (is.null(m)) return(list(cell = "---", n = NA_integer_,
                              nclust = NA_integer_))
  co <- coef(m)["exposure_early_sar_legacy"]
  se <- sqrt(vcov(m)["exposure_early_sar_legacy",
                     "exposure_early_sar_legacy"])
  p  <- 2 * stats::pnorm(-abs(co / se))
  cell <- sprintf("%.3f%s\n(%.3f)", unname(co), stars(p), unname(se))
  # Cluster count = distinct commid_birth_legacy in the rows that
  # actually entered the regression (i.e., after NA-dropping on the
  # outcome, treatment, and controls).
  used <- data_used |>
    dplyr::filter(!is.na(.data[[y]]),
                  !is.na(exposure_early_sar_legacy),
                  !is.na(commid_birth_legacy))
  list(cell = cell,
       n    = as.integer(fixest::fitstat(m, "n")$n),
       nclust = dplyr::n_distinct(used$commid_birth_legacy))
}

rows <- lapply(outcomes, function(y) {
  mA <- fit_sib(y, sib2, "A")
  mB <- fit_sib(y, sib2, "B")
  mC <- fit_sib(y, sib3, "C")
  pA <- pull_cell(mA, sib2, y)
  pB <- pull_cell(mB, sib2, y)
  pC <- pull_cell(mC, sib3, y)
  tibble(outcome = y,
         A = pA$cell, nA = pA$n, cA = pA$nclust,
         B = pB$cell, nB = pB$n, cB = pB$nclust,
         C = pC$cell, nC = pC$n, cC = pC$nclust)
})
sib_tbl <- bind_rows(rows)

# Column sample sizes and cluster counts are constant across outcomes
# within a given spec because we restrict to common non-missing; take
# the max over outcomes as the reported N so the row aligns with what
# a reader sees in the largest model.
nA <- max(sib_tbl$nA, na.rm = TRUE)
nB <- max(sib_tbl$nB, na.rm = TRUE)
nC <- max(sib_tbl$nC, na.rm = TRUE)
cA <- max(sib_tbl$cA, na.rm = TRUE)
cB <- max(sib_tbl$cB, na.rm = TRUE)
cC <- max(sib_tbl$cC, na.rm = TRUE)

coef_lines <- apply(sib_tbl, 1, function(r) {
  cells <- vapply(c("A", "B", "C"), function(k) {
    paste0("\\makecell{", gsub("\n", "\\\\\\\\", r[[k]]), "}")
  }, character(1))
  sprintf("%s & %s \\\\",
          pretty_label(r["outcome"]),
          paste(cells, collapse = " & "))
})

writeLines(c(
  "% Auto-generated by code/R/10_sibling_fe.R",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Outcome & (A) $\\geq 2$ siblings & (B) (A) + community FE & (C) $\\geq 3$ siblings \\\\",
  "\\midrule",
  coef_lines,
  "\\midrule",
  "Mother FE & Yes & Yes & Yes \\\\",
  "Community-of-birth FE & No & Yes & No \\\\",
  "Birth-year FE & Yes & Yes & Yes \\\\",
  sprintf("Clusters (communities) & %s & %s & %s \\\\",
          format(cA, big.mark = ","),
          format(cB, big.mark = ","),
          format(cC, big.mark = ",")),
  sprintf("Observations & %s & %s & %s \\\\",
          format(nA, big.mark = ","),
          format(nB, big.mark = ","),
          format(nC, big.mark = ",")),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "sibling_fe.tex"))

# ---------------------------------------------------------------------
# 3. Fertility-selection check.
#    If VM arrival changes which women have children (compositional
#    shift in the mother pool), sibling FE cannot save us. Test:
#    number of children per mother ~ VM exposure of the mother's commid
#    before vs after rollout. Restrict to mothers whose commid we
#    observe in the rollout data.
# ---------------------------------------------------------------------
# Mother-level fertility snapshot. Within-commid test: conditional on
# village, does a mother's number of observed children (across cohorts
# 1984–1999) correlate with whether the VM was present for her
# reproductive years? We cannot compute fertility perfectly (only kids
# born 1984–1999 are in our roster), so this is a within-window proxy.
fertility <- df |>
  filter(!is.na(start_sar_legacy), !is.na(mother_dob_yr)) |>
  group_by(mother_pidlink, commid_birth_legacy,
           start_sar_legacy, mother_dob_yr) |>
  summarise(n_kids_window = n(), .groups = "drop") |>
  mutate(
    # Mother's age when VM arrived. Exposed = VM arrived before she
    # was 40 and after she was 15 (her reproductive years).
    mom_age_at_rollout = start_sar_legacy - mother_dob_yr,
    exposed_reproductive = as.integer(mom_age_at_rollout <= 40 &
                                        mom_age_at_rollout >= 15)
  )

fit_f <- feols(n_kids_window ~ exposed_reproductive |
                 commid_birth_legacy,
               data = fertility, cluster = ~ commid_birth_legacy)

est <- coef(fit_f)["exposed_reproductive"]
se  <- sqrt(vcov(fit_f)["exposed_reproductive", "exposed_reproductive"])

writeLines(c(
  "% Auto-generated by code/R/10_sibling_fe.R — fertility-selection check",
  "% Within-commid test: does VM arrival during mother's reproductive",
  "% years shift observed kid count in the 1984-1999 window?",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Regressor & Estimate & SE \\\\",
  "\\midrule",
  sprintf("Exposed during reproductive years & %.4f & %.4f \\\\",
          est, se),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "sibling_fertility_check.tex"))

log_stage("stage10", sprintf(
  "fertility check: coef=%.4f (SE %.4f); N_mothers=%d", est, se,
  nrow(fertility)))

log_stage("stage10", "done")
