# =====================================================================
# 06_sample_weights_indices.R — Regression sample + Anderson indices +
#   IPW weights + VM exit-year (second-shock construction).
#
# Inputs :  data/intermediate/stage05/analysis_frame.rds
#           data/intermediate/stage03/midwife_rollout.rds
#           data/intermediate/stage03/other_facilities.rds
#
# Outputs:  data/intermediate/stage06/analysis_sample.rds
#             — all 11,117 rows of analysis_frame, augmented with:
#                 primary_sample (0/1), Anderson indices
#                 (health_index, cognition_index, depression_index,
#                  bigfive_index), ipw (trimmed).
#           data/intermediate/stage06/ipw_diagnostics.rds
#             — propensity-model coefficients + weight summary.
#
# Identification notes:
#   * Primary sample: birth_year ∈ 1984–1999, non-missing
#     commid_birth_preceding, dead == 0. Legacy and mother variants are
#     kept for robustness specs (stage 11) — all three flags computed.
#   * IPW is DEMOTED from identification strategy to a robustness column
#     `early_treated` (start_sar_preceding ≤ median among treated)
#     regressed on province dummies + post-period other-facility
#     presence (imperfect proxy for pre-program village traits —
#     IFLS W1 has no usable pre-1989 community module).
#   * The two-shock / AFC-exit design was explored here originally
#     (see git history). Dropped 2026-04-15: SAR exit-year fields are
#     empty across all waves, leaving only 3-year-resolution PKK
#     presence-gap windows; identification on ~36 coarse events was
#     too weak to support a headline. Code removed to simplify the
#     pipeline; vm_exits.rds is no longer emitted.
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))
source(here::here("code", "R", "lib", "indices.R"))

log_stage("stage06", "begin")

# ---------------------------------------------------------------------
# 1. Load upstream.
# ---------------------------------------------------------------------
frame   <- read_intermediate("stage05", "analysis_frame")
rollout <- read_intermediate("stage03", "midwife_rollout") |>
  select(commid93, start_year_sar, start_year_pkk,
         present_any_sar, present_any_pkk)
other_f <- read_intermediate("stage03", "other_facilities")

# ---------------------------------------------------------------------
# 2. Primary-sample flag (kept as a column; no rows dropped).
# ---------------------------------------------------------------------
frame <- frame |>
  mutate(
    # Headline sample = non-migrant children (legacy birthplace variant).
    # Restricts to children whose own community identifier is unchanged
    # across all four IFLS waves, so birthplace is cleanly identified.
    primary_sample = as.integer(
      birth_year >= 1984 & birth_year <= 1999 &
      !is.na(commid_birth_legacy) &
      dead == 0
    ),
    # Companion flags for robustness in stage 11.
    preceding_sample = as.integer(
      birth_year >= 1984 & birth_year <= 1999 &
      !is.na(commid_birth_preceding) & dead == 0
    ),
    mother_sample = as.integer(
      birth_year >= 1984 & birth_year <= 1999 &
      !is.na(commid_birth_mother) & dead == 0
    )
  )

log_stage("stage06", sprintf(
  "sample flags: primary=%d legacy=%d mother=%d (of %d rows)",
  sum(frame$primary_sample), sum(frame$legacy_sample),
  sum(frame$mother_sample), nrow(frame)))

# ---------------------------------------------------------------------
# 3. Anderson composite indices.
#    Sign convention: higher = better. Flip where needed.
# ---------------------------------------------------------------------
# Health — general adult health + early-life health retrospection.
#   height_cm, no_underweight (binary, higher = better), lung_log,
#   no_hypertension (binary, higher = better), good_early_health
#   (binary), ghs (general health self-rating, 1–4 where 1=very healthy;
#   flip so higher = better).
health_inputs <- orient_components(frame, list(
  list(var = "height_cm",        flip = FALSE),
  list(var = "no_underweight",   flip = FALSE),
  list(var = "lung_log",         flip = FALSE),
  list(var = "no_hypertension",  flip = FALSE),
  list(var = "good_early_health",flip = FALSE),
  list(var = "ghs",              flip = TRUE)
))
frame$health_index <- anderson_index(
  as.data.frame(health_inputs), names(health_inputs)
)

# Cognition — adult IRT primary + supporting recall/serial-7 counts.
#   w_abil (continuous, higher = better), serial7_correct (0–5),
#   immed_count (0–10 memory recall), delay_count (0–10 delayed recall).
cog_inputs <- orient_components(frame, list(
  list(var = "w_abil",          flip = FALSE),
  list(var = "serial7_correct", flip = FALSE),
  list(var = "immed_count",     flip = FALSE),
  list(var = "delay_count",     flip = FALSE)
))
frame$cognition_index <- anderson_index(
  as.data.frame(cog_inputs), names(cog_inputs)
)

# Depression — use CESD-10 dep_index directly (already a sum of items);
# z-score AND flip so higher = less depressed (better).
frame$depression_index <- zscore(-frame$dep_index)

# Big Five — composite over the 5 trait z-scores already in the frame.
# These are separate constructs so the Anderson weighting will mostly
# equal-weight them (covariance is modest); index is reported alongside
# trait-by-trait results in stage 08.
bf_inputs <- orient_components(frame, list(
  list(var = "agreeableness_z",      flip = FALSE),
  list(var = "conscientiousness_z",  flip = FALSE),
  list(var = "extraversion_z",       flip = FALSE),
  list(var = "emotional_stability_z",flip = FALSE),
  list(var = "openness_z",           flip = FALSE)
))
frame$bigfive_index <- anderson_index(
  as.data.frame(bf_inputs), names(bf_inputs)
)

log_stage("stage06", sprintf(
  "Anderson indices — health:%d nonNA cognition:%d depression:%d bigfive:%d",
  sum(!is.na(frame$health_index)),
  sum(!is.na(frame$cognition_index)),
  sum(!is.na(frame$depression_index)),
  sum(!is.na(frame$bigfive_index))))

# ---------------------------------------------------------------------
# 4. IPW weights (demoted, one robustness column).
#    Propensity model at the individual level (exposure assignment is
#    effectively at the commid × cohort level but weights attach to
#    individual rows). Model: binary "early-treated" ~ province dummies
#    + post-period other-facility presence (imperfect proxy for pre-
#    program village traits — IFLS W1 has no pre-1989 community module).
#    Weights trimmed at 1st/99th pct and centered so the weighted sample
#    size equals the unweighted N.
# ---------------------------------------------------------------------
# Define early-treated = exposure_early_sar_preceding == 1 (exposed
# within age 3 or in utero). The propensity is P(early | X).
ipw_df <- frame |>
  filter(primary_sample == 1,
         !is.na(exposure_early_sar_preceding),
         !is.na(province)) |>
  mutate(
    early = as.integer(exposure_early_sar_preceding),
    prov  = factor(province),
    # Imperfect pre-program proxies (post-period facility presence):
    doc   = as.integer(coalesce(has_doctor, 0)),
    clin  = as.integer(coalesce(has_clinic, 0)),
    pm    = as.integer(coalesce(has_private_midwife, 0)),
    tm    = as.integer(coalesce(has_traditional_midwife, 0))
  )

pscore_fit <- stats::glm(
  early ~ prov + doc + clin + pm + tm,
  data = ipw_df, family = stats::binomial(link = "logit")
)
ipw_df$pscore <- stats::predict(pscore_fit, type = "response")
ipw_df$ipw_raw <- with(ipw_df, ifelse(early == 1, 1 / pscore,
                                      1 / (1 - pscore)))
# Trim tails at 1st and 99th percentile.
q <- stats::quantile(ipw_df$ipw_raw, probs = c(0.01, 0.99), na.rm = TRUE)
ipw_df$ipw <- pmin(pmax(ipw_df$ipw_raw, q[1]), q[2])
# Center so mean(ipw) = 1 (stabilises weighted N).
ipw_df$ipw <- ipw_df$ipw / mean(ipw_df$ipw, na.rm = TRUE)

frame <- frame |>
  left_join(ipw_df |> select(pidlink, ipw), by = "pidlink")

log_stage("stage06", sprintf(
  "IPW: n=%d mean=%.3f sd=%.3f trimmed at [%.2f, %.2f]",
  sum(!is.na(frame$ipw)), mean(frame$ipw, na.rm = TRUE),
  stats::sd(frame$ipw, na.rm = TRUE), q[1], q[2]))

ipw_diag <- list(
  coefs = broom::tidy(pscore_fit),
  trim  = c(p01 = unname(q[1]), p99 = unname(q[2])),
  n     = nrow(ipw_df)
)

# ---------------------------------------------------------------------
# 5. Persist.
# ---------------------------------------------------------------------
write_intermediate(frame,    "stage06", "analysis_sample")
write_intermediate(ipw_diag, "stage06", "ipw_diagnostics")

log_stage("stage06", "done")
