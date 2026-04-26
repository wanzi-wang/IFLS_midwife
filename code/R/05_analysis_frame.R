# =====================================================================
# 05_analysis_frame.R — Final merge + exposure timing + diagnostics
#
# Produces `data/intermediate/stage05/analysis_frame.rds`: one row per
# pidlink in the unrestricted birth roster (cohorts 1984-1999), with:
#   - identifiers + demographics + mother covariates (stage 02)
#   - start_year_sar, start_year_pkk (stage 03) joined on three
#     commid_birth variants
#   - exposure timing vars computed for each (rollout source) ×
#     (commid_birth variant) combination
#   - outcomes from stage 04
#   - other-facility controls from stage 03
#
# Exposure timing variables (per SAR/PKK × legacy/mother/preceding):
#   - age_at_rollout_{src}_{var} = start_year - birth_year
#   - inutero_{src}_{var}        = start_year <= birth_year - 1
#   - within_3_{src}_{var}       = !inutero & age ∈ {0..3}
#   - within_46_{src}_{var}      = !inutero & !within_3 & age ∈ {4..6}
#   - exposure_early_{src}_{var} = (age ≤ 3) | inutero
#
# Diagnostics logged (feed regression plan):
#   - siblings-per-mother distribution (prof comment #9: power check for
#     sibling FE)
#   - mortality × exposure cross-tab (selection check)
#   - SAR-vs-PKK start_year disagreement on attached commids
#
# Ends with reconciliation gate (lib/audit.R) vs legacy df.rds.
# =====================================================================

set.seed(20260412)

library(here)
library(dplyr)
library(tidyr)
library(tibble)

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))
source(here::here("code", "R", "lib", "audit.R"))

log_stage("stage05", "begin")

# ---------------------------------------------------------------------
# 1. Load upstream intermediates.
# ---------------------------------------------------------------------
br       <- read_intermediate("stage02", "birth_roster")
rollout  <- read_intermediate("stage03", "midwife_rollout") |>
  select(commid93, start_year_sar, start_year_pkk,
         present_any_sar, present_any_pkk)
other_fac <- read_intermediate("stage03", "other_facilities")
outcomes  <- read_intermediate("stage04", "outcomes")
mediators <- read_intermediate("stage04b", "mediators") |>
  select(-wave_source)

# ---------------------------------------------------------------------
# 2. Attach rollout (SAR + PKK) for each birthplace variant.
#    Each variant gets its own pair of start-year columns.
# ---------------------------------------------------------------------
attach_rollout <- function(df, birth_col, suffix) {
  sar_col <- paste0("start_sar_", suffix)
  pkk_col <- paste0("start_pkk_", suffix)
  df |>
    left_join(
      rollout |>
        rename(!!sar_col := start_year_sar,
               !!pkk_col := start_year_pkk) |>
        select(commid93, !!sar_col, !!pkk_col),
      by = c(stats::setNames("commid93", birth_col))
    ) |>
    # Coalesce SAR -> PKK so the main exposure variable uses the best
    # available rollout year per community. SAR is preferred (it tags
    # facility-level rollout dates); PKK fills in communities where the
    # SAR year is missing but the PKK year is available. The PKK-only
    # variant `start_pkk_<suffix>` is preserved for the R4 robustness
    # check in stage 11. Without this coalesce, ~740 individuals in
    # 62 communities are dropped from the main analysis because their
    # SAR rollout year is missing even though PKK records it.
    mutate(!!sar_col := dplyr::coalesce(.data[[sar_col]], .data[[pkk_col]]))
}

frame <- br |>
  attach_rollout("commid_birth_legacy",    "legacy") |>
  attach_rollout("commid_birth_mother",    "mother") |>
  attach_rollout("commid_birth_preceding", "preceding")

# ---------------------------------------------------------------------
# 3. Exposure-timing variables per (src, variant).
# ---------------------------------------------------------------------
build_exposure <- function(df, src, variant) {
  sy <- paste0("start_", src, "_", variant)
  df |>
    mutate(
      !!paste0("age_at_rollout_", src, "_", variant) :=
        .data[[sy]] - birth_year,
      !!paste0("inutero_", src, "_", variant) :=
        as.integer(.data[[sy]] <= birth_year - 1),
      !!paste0("exposure_early_", src, "_", variant) :=
        as.integer(
          (.data[[sy]] - birth_year) <= 3 |
          (.data[[sy]] <= birth_year - 1)
        ),
      !!paste0("within_3_", src, "_", variant) :=
        as.integer(
          !(.data[[sy]] <= birth_year - 1) &
          (.data[[sy]] - birth_year) %in% 0:3
        ),
      !!paste0("within_46_", src, "_", variant) :=
        as.integer(
          !(.data[[sy]] <= birth_year - 1) &
          !((.data[[sy]] - birth_year) %in% 0:3) &
          (.data[[sy]] - birth_year) %in% 4:6
        )
    )
}

for (src in c("sar", "pkk")) {
  for (variant in c("legacy", "mother", "preceding")) {
    frame <- build_exposure(frame, src, variant)
  }
}

# ---------------------------------------------------------------------
# 4. Attach outcomes and other-facility controls.
# ---------------------------------------------------------------------
frame <- frame |>
  left_join(outcomes, by = "pidlink") |>
  left_join(other_fac, by = c("commid_birth_legacy" = "commid93")) |>
  left_join(mediators, by = c("mother_pidlink", "birth_year"))

# ---------------------------------------------------------------------
# 5. Diagnostics.
# ---------------------------------------------------------------------
# Siblings-per-mother: informs prof comment #9 (sibling-FE feasibility).
sib_dist <- frame |>
  filter(!is.na(mother_pidlink)) |>
  count(mother_pidlink, name = "n_sibs") |>
  count(n_sibs, name = "n_mothers") |>
  arrange(n_sibs)

log_stage("stage05", sprintf(
  "siblings-per-mother: %d mothers; %d ≥2 sibs; %d ≥3 sibs",
  sum(sib_dist$n_mothers),
  sum(sib_dist$n_mothers[sib_dist$n_sibs >= 2]),
  sum(sib_dist$n_mothers[sib_dist$n_sibs >= 3])
))

# Mortality × exposure_early (SAR×legacy variant, the thesis's baseline).
mort_x_exp <- frame |>
  filter(birth_year %in% 1984:1999, non_migrant == 1,
         !is.na(exposure_early_sar_legacy)) |>
  count(dead, exposure_early_sar_legacy) |>
  arrange(dead, exposure_early_sar_legacy)

log_stage("stage05", paste0(
  "mortality × exposure_early_sar_legacy: ",
  paste(sprintf("dead=%d,expo=%d:n=%d", mort_x_exp$dead,
                mort_x_exp$exposure_early_sar_legacy, mort_x_exp$n),
        collapse = "; ")
))

# Analytic-cohort counts per variant for the plan's sample-tracking log.
cohort_counts <- c(
  legacy    = sum(!is.na(frame$commid_birth_legacy) &
                    frame$birth_year %in% 1984:1999 & frame$dead == 0),
  mother    = sum(!is.na(frame$commid_birth_mother) &
                    frame$birth_year %in% 1984:1999 & frame$dead == 0),
  preceding = sum(!is.na(frame$commid_birth_preceding) &
                    frame$birth_year %in% 1984:1999 & frame$dead == 0)
)
log_stage("stage05", sprintf(
  "analytic-cohort n per birthplace variant: legacy=%d mother=%d preceding=%d",
  cohort_counts["legacy"], cohort_counts["mother"], cohort_counts["preceding"]
))

# Integrity: no NA in commid_birth_legacy among the thesis-matched sample.
stopifnot(sum(is.na(frame$commid_birth_legacy) &
              frame$non_migrant == 1) == 0)

# ---------------------------------------------------------------------
# 6. Write analysis frame + diagnostics.
# ---------------------------------------------------------------------
write_intermediate(frame,     "stage05", "analysis_frame")
write_intermediate(sib_dist,  "stage05", "siblings_per_mother")
write_intermediate(mort_x_exp,"stage05", "mortality_by_exposure")

# Manifest of columns for quick reference.
manifest <- tibble::tibble(
  column = names(frame),
  class  = vapply(frame, function(x) class(x)[1], character(1)),
  n_non_na = vapply(frame, function(x) sum(!is.na(x)), integer(1))
)
write_intermediate(manifest, "stage05", "analysis_frame_manifest")

# ---------------------------------------------------------------------
# 7. Reconciliation gate vs legacy df.rds.
# ---------------------------------------------------------------------
legacy_path <- "/Users/wwz/Academic/Thesis_Data_IFLS/old/Intermediate/Reg/df.rds"
audit_frame(frame, legacy_path)

log_stage("stage05", "done")
