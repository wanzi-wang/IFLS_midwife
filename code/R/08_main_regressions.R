# =====================================================================
# 08_main_regressions.R — Main specifications on the primary sample.
#
# Runs five specifications per primary outcome:
#   (1) TWFE baseline
#   (2) TWFE + mother covariates + sex
#   (3) TWFE + province × birth-year FE
#   (4) Callaway-Sant'Anna (did::att_gt, event-time 0-3 average)
#   (5) de Chaisemartin-D'Haultfoeuille (DIDmultiplegtDYN)
#
# Primary outcomes (Anderson indices built in stage 06):
#   health_index, cognition_index, depression_index, bigfive_index.
#
# Treatment:    exposure_early_sar_preceding
# Fixed effects: commid_birth_preceding + birth_year
# Clustering:   commid_birth_preceding
#
# For staggered-DiD-robust methods (CS, dCDH), treatment at the
# COMMUNITY level. Group G = first-treated cohort = start_sar_preceding
# - 1 (the in-utero cohort when the VM arrived). Never-treated villages
# act as the control group. The ATT is averaged over event-times 0 to 3
# to match the "age 0-3 or in utero" exposure definition.
#
# Outputs:  paper/tables/main_results.tex
#           paper/tables/main_results_by_sex.tex
#           data/output/regressions/stage08_models.rds
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(fixest)
  library(modelsummary)
  library(did)
  library(DIDmultiplegtDYN)
  if (requireNamespace("polars", quietly = TRUE)) library(polars)
})

options(modelsummary_factory_latex = "kableExtra",
        modelsummary_format_numeric_latex = "plain")

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))

log_stage("stage08", "begin")

frame <- read_intermediate("stage06", "analysis_sample")

tables_dir <- here::here("paper", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
reg_dir    <- here::here("data", "output", "regressions")
dir.create(reg_dir,    recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 1. Primary-sample frame.
#    G (first-treated cohort) = start_sar_preceding - 1. Never-treated
#    villages have G = 0 (sentinel used by did::att_gt).
# ---------------------------------------------------------------------
df <- frame |>
  filter(primary_sample == 1) |>
  mutate(
    G = if_else(!is.na(start_sar_preceding),
                as.integer(start_sar_preceding) - 1L, 0L),
    # Numeric village IDs for did::att_gt (requires integer `idname`).
    commid_num = as.integer(factor(commid_birth_preceding)),
    sex_f = factor(sex, levels = c(1, 3), labels = c("male", "female"))
  )

log_stage("stage08", sprintf(
  "primary-sample N=%d; %d villages; %d never-treated (G=0); birth_year %d-%d.",
  nrow(df), dplyr::n_distinct(df$commid_birth_preceding),
  sum(df$G == 0), min(df$birth_year), max(df$birth_year)))

outcomes <- c("health_index", "cognition_index",
              "depression_index", "bigfive_index")

# ---------------------------------------------------------------------
# 2. TWFE specifications (models 1-3 per outcome).
# ---------------------------------------------------------------------
fit_twfe_set <- function(y, data) {
  m1 <- feols(
    as.formula(sprintf(
      "%s ~ exposure_early_sar_preceding | commid_birth_preceding + birth_year",
      y)),
    data = data, cluster = ~ commid_birth_preceding
  )
  m2 <- feols(
    as.formula(sprintf(
      paste0("%s ~ exposure_early_sar_preceding + sex_f + mother_edu_years + ",
             "mother_age_birth | commid_birth_preceding + birth_year"),
      y)),
    data = data, cluster = ~ commid_birth_preceding
  )
  m3 <- feols(
    as.formula(sprintf(
      paste0("%s ~ exposure_early_sar_preceding + sex_f + mother_edu_years + ",
             "mother_age_birth | commid_birth_preceding + birth_year + ",
             "province^birth_year"),
      y)),
    data = data, cluster = ~ commid_birth_preceding
  )
  list(TWFE = m1, TWFE_X = m2, `TWFE_X+prov-yr` = m3)
}

twfe_models <- setNames(vector("list", length(outcomes)), outcomes)
for (y in outcomes) {
  log_stage("stage08", sprintf("fitting TWFE for %s", y))
  twfe_models[[y]] <- fit_twfe_set(y, df)
}

# ---------------------------------------------------------------------
# 3. Callaway-Sant'Anna: ATT averaged over event-time 0:3.
# ---------------------------------------------------------------------
fit_cs <- function(y, data) {
  d <- data |>
    filter(!is.na(.data[[y]]), !is.na(birth_year), !is.na(commid_num)) |>
    as.data.frame()
  tryCatch({
    att <- att_gt(
      yname = y, tname = "birth_year", idname = "commid_num",
      gname = "G", data = d, control_group = "nevertreated",
      bstrap = FALSE, cband = FALSE,
      allow_unbalanced_panel = TRUE, panel = FALSE
    )
    dyn <- aggte(att, type = "dynamic", min_e = 0, max_e = 3,
                 na.rm = TRUE)
    list(att = att, dyn = dyn,
         coef = dyn$overall.att, se = dyn$overall.se)
  }, error = function(e) {
    message(sprintf("CS failed for %s: %s", y, e$message))
    list(att = NULL, dyn = NULL, coef = NA_real_, se = NA_real_)
  })
}

cs_models <- setNames(vector("list", length(outcomes)), outcomes)
for (y in outcomes) {
  log_stage("stage08", sprintf("fitting CS for %s", y))
  cs_models[[y]] <- fit_cs(y, df)
}

# ---------------------------------------------------------------------
# 4. de Chaisemartin-D'Haultfoeuille (dCDH).
#    `did_multiplegt_dyn` needs a treatment indicator that is 1 for
#    cohorts in the "post" period per village. We use a "year-of-birth
#    in-or-after G" variable as the instantaneous treatment.
# ---------------------------------------------------------------------
fit_dcdh <- function(y, data) {
  d <- data |>
    filter(!is.na(.data[[y]])) |>
    mutate(treat_indicator = as.integer(!is.na(G) & G > 0 &
                                          birth_year >= G)) |>
    select(commid_num, birth_year, all_of(y), treat_indicator) |>
    as.data.frame()
  tryCatch({
    fit <- did_multiplegt_dyn(
      df = d, outcome = y, group = "commid_num",
      time = "birth_year", treatment = "treat_indicator",
      effects = 4, placebo = 0, graph_off = TRUE,
      cluster = "commid_num"
    )
    # Overall dynamic effect (avg of event-times 0 to 3).
    eff <- fit$results$Effects
    list(fit = fit,
         coef = mean(eff[, "Estimate"], na.rm = TRUE),
         se   = sqrt(mean(eff[, "SE"]^2, na.rm = TRUE)))
  }, error = function(e) {
    message(sprintf("dCDH failed for %s: %s", y, e$message))
    list(fit = NULL, coef = NA_real_, se = NA_real_)
  })
}

dcdh_models <- setNames(vector("list", length(outcomes)), outcomes)
for (y in outcomes) {
  log_stage("stage08", sprintf("fitting dCDH for %s", y))
  dcdh_models[[y]] <- fit_dcdh(y, df)
}

# ---------------------------------------------------------------------
# 5. Assemble main-results table.
#    Rows = outcomes; columns = TWFE (1-3), CS, dCDH. Report exposure-
#    early coefficient + cluster SE.
# ---------------------------------------------------------------------
# Pull coef+SE from TWFE models via fixest.
twfe_summary <- function(mods) {
  lapply(mods, function(m) {
    co <- coef(m)["exposure_early_sar_preceding"]
    se <- sqrt(vcov(m)["exposure_early_sar_preceding",
                       "exposure_early_sar_preceding"])
    c(est = unname(co), se = unname(se))
  })
}

main_rows <- lapply(outcomes, function(y) {
  ts <- twfe_summary(twfe_models[[y]])
  tibble(
    outcome = y,
    TWFE    = sprintf("%.3f\n(%.3f)", ts$TWFE["est"], ts$TWFE["se"]),
    TWFE_X  = sprintf("%.3f\n(%.3f)", ts$TWFE_X["est"], ts$TWFE_X["se"]),
    `TWFE_X+prov-yr` = sprintf("%.3f\n(%.3f)",
                               ts$`TWFE_X+prov-yr`["est"],
                               ts$`TWFE_X+prov-yr`["se"]),
    CS      = if (!is.na(cs_models[[y]]$coef))
      sprintf("%.3f\n(%.3f)",
              cs_models[[y]]$coef, cs_models[[y]]$se) else "---",
    dCDH    = if (!is.na(dcdh_models[[y]]$coef))
      sprintf("%.3f\n(%.3f)",
              dcdh_models[[y]]$coef, dcdh_models[[y]]$se) else "---"
  )
})
main_tbl <- bind_rows(main_rows)

# Write as simple LaTeX (newline inside cell rendered via \makecell).
writeLines(c(
  "% Auto-generated by code/R/08_main_regressions.R",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Outcome & (1) TWFE & (2) +covariates & (3) +prov$\\times$yr & (4) CS & (5) dCDH \\\\",
  "\\midrule",
  apply(main_tbl, 1, function(r) {
    cells <- vapply(r[-1], function(cell) {
      paste0("\\makecell{", gsub("\n", "\\\\\\\\", cell), "}")
    }, character(1))
    sprintf("%s & %s \\\\",
            gsub("_", "\\\\_", r["outcome"]),
            paste(cells, collapse = " & "))
  }),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "main_results.tex"))

# ---------------------------------------------------------------------
# 6. Sex heterogeneity — run TWFE_X by sex, coefficient-only table.
# ---------------------------------------------------------------------
sex_panel <- lapply(outcomes, function(y) {
  mm <- list(
    male   = feols(
      as.formula(sprintf(paste0(
        "%s ~ exposure_early_sar_preceding + mother_edu_years + ",
        "mother_age_birth | commid_birth_preceding + birth_year"), y)),
      data = df |> filter(sex == 1),
      cluster = ~ commid_birth_preceding
    ),
    female = feols(
      as.formula(sprintf(paste0(
        "%s ~ exposure_early_sar_preceding + mother_edu_years + ",
        "mother_age_birth | commid_birth_preceding + birth_year"), y)),
      data = df |> filter(sex == 3),
      cluster = ~ commid_birth_preceding
    )
  )
  ts <- lapply(mm, function(m) {
    co <- coef(m)["exposure_early_sar_preceding"]
    se <- sqrt(vcov(m)["exposure_early_sar_preceding",
                       "exposure_early_sar_preceding"])
    c(est = unname(co), se = unname(se))
  })
  tibble(
    outcome = y,
    Male    = sprintf("%.3f\n(%.3f)", ts$male["est"],   ts$male["se"]),
    Female  = sprintf("%.3f\n(%.3f)", ts$female["est"], ts$female["se"])
  )
})
sex_tbl <- bind_rows(sex_panel)

writeLines(c(
  "% Auto-generated by code/R/08_main_regressions.R",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Outcome & Male & Female \\\\",
  "\\midrule",
  apply(sex_tbl, 1, function(r) {
    cells <- vapply(r[-1], function(cell) {
      paste0("\\makecell{", gsub("\n", "\\\\\\\\", cell), "}")
    }, character(1))
    sprintf("%s & %s \\\\",
            gsub("_", "\\\\_", r["outcome"]),
            paste(cells, collapse = " & "))
  }),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "main_results_by_sex.tex"))

# ---------------------------------------------------------------------
# 7. Persist fitted objects for downstream stages (09 event study).
# ---------------------------------------------------------------------
saveRDS(list(twfe = twfe_models, cs = cs_models, dcdh = dcdh_models),
        file.path(reg_dir, "stage08_models.rds"), compress = "xz")

log_stage("stage08", "main tables written.")
log_stage("stage08", "done")
