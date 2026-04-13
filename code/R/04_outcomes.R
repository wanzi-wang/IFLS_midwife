# =====================================================================
# 04_outcomes.R — W5 outcome variables per pidlink (raw only, no indices)
#
# Pulls health, cognition, depression, and Big-5 personality outcomes
# from W5 modules, joins all on `pidlink`, writes a single tibble.
#
# Explicit non-goals (deferred to regression plan):
#   * Anderson composite indices (depend on control-group covariance,
#     which depends on sample restrictions).
#   * Within-subset z-scoring (same).
#   * Sample restrictions (cohort, non-migrant, alive). This file emits
#     every pidlink we can build outcomes for; stage 05 & regression
#     stage subset to the analysis cohort.
#
# Audit notes (upgrade over legacy):
#   * Legacy used Book EK (child module, ages 7-14) for `cog_prop` and
#     `math_prop`, which truncates our 1984-1999 cohort (mostly aged
#     15-30 in 2014) to only the 1999-born at the boundary — plausibly
#     the source of the prof's "weak power?" comment on cognition
#     results. We compute BOTH:
#       - `cog_prop_child` / `math_prop_child`: from EK (1999-born)
#       - `w_abil`: IFLS-precomputed IRT ability score from Book 3B
#         (adaptive number-series test), available for all adults
#     `w_abil` is the appropriate primary cognition measure for the
#     15-30 adult cohort.
#   * Depression: legacy labelled this "Kessler-10" but the CESD-10
#     items A..J (reverse-coded E,H) are actually the CESD-10 short
#     form per IFLS documentation. We keep the legacy reverse-code
#     mapping but rename the metric accordingly.
# =====================================================================

set.seed(20260412)

library(here)
library(dplyr)
library(tidyr)
library(tibble)

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))

log_stage("stage04", "begin")

# ---------------------------------------------------------------------
# 1. Health (bus_us): height, weight, BMI, BP, lung, GHS.
# ---------------------------------------------------------------------
health <- read_intermediate("raw", "W5__health") |>
  transmute(
    pidlink = as.character(pidlink),
    height_cm   = as.numeric(us04),
    weight_kg   = as.numeric(us06),
    bp_syst_1   = as.numeric(us07a1),
    bp_syst_2   = as.numeric(us07b1),
    bp_syst_3   = as.numeric(us07c1),
    bp_diast_1  = as.numeric(us07a2),
    bp_diast_2  = as.numeric(us07b2),
    bp_diast_3  = as.numeric(us07c2),
    lung_1      = as.numeric(us09a),
    lung_2      = as.numeric(us09b),
    ghs         = as.numeric(us14)
  ) |>
  mutate(
    # Sanity bounds — implausible measurements drop to NA.
    height_cm  = if_else(height_cm >= 100 & height_cm <= 210,
                         height_cm, NA_real_),
    weight_kg  = if_else(weight_kg >= 15 & weight_kg <= 200,
                         weight_kg, NA_real_),
    bmi            = weight_kg / (height_cm / 100)^2,
    no_underweight = as.integer(bmi >= 18.5),
    # Average of 2nd and 3rd BP readings (first is discarded in
    # clinical practice to avoid white-coat bias; matches legacy).
    bp_syst_mean   = rowMeans(cbind(bp_syst_2, bp_syst_3), na.rm = TRUE),
    bp_diast_mean  = rowMeans(cbind(bp_diast_2, bp_diast_3), na.rm = TRUE),
    no_hypertension = as.integer(bp_syst_mean <= 140 & bp_diast_mean <= 90),
    # Lung function: log mean of the two peak-expiratory-flow readings.
    lung_mean = rowMeans(cbind(lung_1, lung_2), na.rm = TRUE),
    lung_log  = if_else(lung_mean > 0, log(lung_mean), NA_real_),
    # GHS is 1-9 self-report; sanity bounds.
    ghs = if_else(ghs %in% 1:9, ghs, NA_real_)
  ) |>
  select(pidlink, height_cm, weight_kg, bmi, no_underweight,
         bp_syst_mean, bp_diast_mean, no_hypertension, lung_log, ghs) |>
  distinct(pidlink, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 2. Early-life health (b3b_eh): eh01 self-rated; eh02-04 missed school/
#    confined/hospitalized for ≥1 month.
# ---------------------------------------------------------------------
early_health <- read_intermediate("raw", "W5__early_health") |>
  transmute(
    pidlink = as.character(pidlink),
    eh01    = as.numeric(eh01),
    eh02    = as.numeric(eh02),
    eh03    = as.numeric(eh03),
    eh04    = as.numeric(eh04)
  ) |>
  mutate(
    good_early_health = as.integer(eh01 %in% c(1, 2, 3)),
    # eh02/03/04: 3 = no, 1 = yes per IFLS codebook.
    notsick_a_month = as.integer(eh02 == 3 & eh03 == 3 & eh04 == 3)
  ) |>
  select(pidlink, eh01, eh02, eh03, eh04,
         good_early_health, notsick_a_month) |>
  distinct(pidlink, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 3. Adult cognition (b3b_co1): date-awareness, serial-7, word recall.
# ---------------------------------------------------------------------
cog_adult1 <- read_intermediate("raw", "W5__cog_adult1") |>
  transmute(
    pidlink     = as.character(pidlink),
    date_all    = as.numeric(co02),      # 1=all components correct
    day_correct = as.numeric(co04),      # 1=correct
    self_memory = as.numeric(co04aa),
    s7_1        = as.numeric(co04a),
    s7_2        = as.numeric(co04b),
    s7_3        = as.numeric(co04c),
    s7_4        = as.numeric(co04d),
    s7_5        = as.numeric(co04e),
    immed_count = as.numeric(co07count),
    delay_count = as.numeric(co10count)
  ) |>
  mutate(
    memory_days = as.integer(date_all == 1 & day_correct == 1),
    # Serial-7 perfect: 100,93,86,79,72,65 — all five answers correct.
    serial7_correct = as.integer(
      s7_1 == 93 & s7_2 == 86 & s7_3 == 79 & s7_4 == 72 & s7_5 == 65
    ),
    self_memory = if_else(self_memory %in% 1:5, self_memory, NA_real_)
  ) |>
  select(pidlink, memory_days, serial7_correct, self_memory,
         immed_count, delay_count) |>
  distinct(pidlink, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 4. Adult adaptive number series (b3b_cob): IRT ability `w_abil`.
#    This is the IFLS-precomputed fluid-intelligence score (~500 mean,
#    SD ~50) — cognition measure appropriate for the 15-30 adult cohort.
# ---------------------------------------------------------------------
cog_adultb <- read_intermediate("raw", "W5__cog_adultb") |>
  transmute(
    pidlink = as.character(pidlink),
    w_abil  = as.numeric(w_abil)
  ) |>
  filter(!is.na(w_abil)) |>
  distinct(pidlink, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 5. Child cognition (ek_ek2): Raven-style items for ages 7-14.
#    Only covers the 1999-born at the analytic cohort boundary. We also
#    pull the EK math items for comparison with legacy.
# ---------------------------------------------------------------------
cog_items_child <- c("ek1_ans", "ek2_ans", "ek3_ans", "ek4_ans",
                     "ek5_ans", "ek6_ans", "ek11_ans", "ek12_ans")
math_items_child <- c("ek18_ans", "ek19_ans", "ek20_ans",
                      "ek21_ans", "ek22_ans")

ek2 <- read_intermediate("raw", "W5__cog_child2") |>
  mutate(across(all_of(c(cog_items_child, math_items_child)),
                as.numeric))

cog_child <- ek2 |>
  transmute(
    pidlink = as.character(pidlink),
    cog_prop_child = rowSums(across(all_of(cog_items_child)),
                              na.rm = TRUE) / length(cog_items_child),
    math_prop_child = rowSums(across(all_of(math_items_child)),
                               na.rm = TRUE) / length(math_items_child)
  ) |>
  distinct(pidlink, .keep_all = TRUE)

# ---------------------------------------------------------------------
# 6. Depression (b3b_kp, CESD-10): long → wide → reverse-code → sum.
#    kptype letters A..J are the 10 items. E & H are "positive-mood"
#    items reverse-scored on a 1-4 Likert scale (rarely..most-days).
# ---------------------------------------------------------------------
depression <- read_intermediate("raw", "W5__depression") |>
  transmute(
    pidlink = as.character(pidlink),
    kptype  = as.character(kptype),
    kp02    = as.numeric(kp02)
  ) |>
  filter(!is.na(kp02), kp02 %in% 1:4) |>
  mutate(
    kp_score = case_when(
      kptype %in% c("A", "B", "C", "D", "F", "G", "I", "J") ~ kp02 - 1,
      kptype %in% c("E", "H")                               ~ 4 - kp02,
      TRUE ~ NA_real_
    )
  ) |>
  group_by(pidlink) |>
  summarise(
    dep_index = sum(kp_score, na.rm = TRUE),
    n_dep_items = sum(!is.na(kp_score)),
    .groups = "drop"
  ) |>
  # Only trust the score if all 10 items are present.
  mutate(
    dep_index = if_else(n_dep_items == 10, dep_index, NA_real_),
    no_depression = as.integer(dep_index < 10)
  )

# Also pivot items wide so the regression plan can compute alternative
# indices (Rasch, 2-factor) without re-reading raw.
depression_items <- read_intermediate("raw", "W5__depression") |>
  transmute(
    pidlink = as.character(pidlink),
    kptype  = paste0("kp02_", as.character(kptype)),
    kp02    = as.numeric(kp02)
  ) |>
  filter(!is.na(kp02), kp02 %in% 1:4) |>
  distinct(pidlink, kptype, .keep_all = TRUE) |>
  pivot_wider(names_from = kptype, values_from = kp02)

depression <- depression |> left_join(depression_items, by = "pidlink")

# ---------------------------------------------------------------------
# 7. Big-5 personality (b3b_psn, 15 items → 5 traits, 3 items each).
#    Legacy reverse-codes psntype 4, 7, 9, 14, 15. Trait → items:
#       agreeableness       = 6, 11, 14
#       conscientiousness   = 2, 9, 12
#       extraversion        = 1, 4, 13
#       emotional_stability = 5, 7, 15 (reverse of neuroticism)
#       openness            = 3, 8, 10
# ---------------------------------------------------------------------
psn_raw <- read_intermediate("raw", "W5__personality") |>
  transmute(
    pidlink = as.character(pidlink),
    psntype = as.integer(psntype),
    psn01   = as.numeric(psn01)
  ) |>
  filter(!is.na(psn01), psn01 %in% 1:5) |>
  mutate(
    score = if_else(psntype %in% c(4, 7, 9, 14, 15),
                    6 - psn01, psn01),
    trait = case_when(
      psntype %in% c(6, 11, 14) ~ "agreeableness",
      psntype %in% c(2, 9, 12)  ~ "conscientiousness",
      psntype %in% c(1, 4, 13)  ~ "extraversion",
      psntype %in% c(5, 7, 15)  ~ "emotional_stability",
      psntype %in% c(3, 8, 10)  ~ "openness",
      TRUE ~ NA_character_
    )
  )

personality <- psn_raw |>
  filter(!is.na(trait)) |>
  group_by(pidlink, trait) |>
  summarise(trait_sum = mean(score, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = trait, values_from = trait_sum)

# Full-sample z-score per trait (regression plan may re-z within subset).
traits <- intersect(
  c("agreeableness", "conscientiousness", "extraversion",
    "emotional_stability", "openness"),
  names(personality)
)
for (tr in traits) {
  z_col <- paste0(tr, "_z")
  personality[[z_col]] <- as.numeric(scale(personality[[tr]]))
}

# Preserve 15 raw items wide as well.
psn_items_wide <- read_intermediate("raw", "W5__personality") |>
  transmute(
    pidlink = as.character(pidlink),
    psntype = paste0("psn_", sprintf("%02d", as.integer(psntype))),
    psn01   = as.numeric(psn01)
  ) |>
  filter(!is.na(psn01), psn01 %in% 1:5) |>
  distinct(pidlink, psntype, .keep_all = TRUE) |>
  pivot_wider(names_from = psntype, values_from = psn01)

personality <- personality |> left_join(psn_items_wide, by = "pidlink")

# ---------------------------------------------------------------------
# 8. Merge everything on pidlink.
# ---------------------------------------------------------------------
outcomes <- health |>
  full_join(early_health, by = "pidlink") |>
  full_join(cog_adult1,   by = "pidlink") |>
  full_join(cog_adultb,   by = "pidlink") |>
  full_join(cog_child,    by = "pidlink") |>
  full_join(depression,   by = "pidlink") |>
  full_join(personality,  by = "pidlink") |>
  distinct(pidlink, .keep_all = TRUE)

log_stage("stage04", sprintf(
  "outcomes rows=%d cols=%d; coverage — height=%d dep_index=%d w_abil=%d",
  nrow(outcomes), ncol(outcomes),
  sum(!is.na(outcomes$height_cm)),
  sum(!is.na(outcomes$dep_index)),
  sum(!is.na(outcomes$w_abil))))

write_intermediate(outcomes, "stage04", "outcomes")

log_stage("stage04", "done")
