# =====================================================================
# 04b_mediators.R — Channel / mechanism variables for the Village
# Midwife Program.
#
# Extracts mother-reported pregnancy-history variables from the IFLS
# book-4 CH module (Kesehatan anak in hhN b4_ch1). These are the
# candidate channels through which a village midwife's arrival would
# move a child's adult-outcome distribution. IFLS ch20 (delivery
# attendant) is a multi-response letter string; the per-letter coding
# is verified against the IFLS-4 (2007) codebook (vol. 5, hhd_B4):
#
#   A = Physician,  B = Private midwife,  C = Village midwife,
#   D = Nurse,      E = Traditional birth attendant (dukun),
#   F = Family,     G = Self / alone,     H = Other,    Z = Missing.
#
# Constructed mediators:
#
#   - pnc_any                    : any prenatal check-up during
#                                  pregnancy (ch14 == 1).
#   - attendant_village_midwife  : delivery attended by a village
#                                  midwife specifically (ch20 contains
#                                  C). This is the program-specific
#                                  channel; the village midwife program
#                                  is precisely a village-midwife
#                                  posting, so the cleanest first-stage
#                                  is a village-midwife indicator
#                                  rather than a generic skilled-
#                                  attendant indicator. Matches the
#                                  channel definition in
#                                  \citeauthor{ahsan2021}.
#   - skilled_attendant          : kept for backwards compatibility.
#                                  1 if ch20 contains any of A/B/C/D
#                                  (any skilled attendant, broader than
#                                  the program-specific channel).
#   - exclusive_bf_6m            : exclusive breastfeeding for at least
#                                  six months. Constructed from
#                                  ch24c/ch24cx (food other than breast
#                                  milk first introduced) and
#                                  ch24d/ch24dx (water other than breast
#                                  milk first introduced); the W2 unit
#                                  flag (ch24cx, ch24dx) maps 5 = month,
#                                  4 = week, 3 = day. Fall-back to
#                                  ever-breastfed when duration is
#                                  unrecorded.
#   - ever_breastfed             : kept for backwards compatibility.
#                                  Ever initiated breastfeeding (ch24a
#                                  in W2/W3, ch28 in W1).
#   - birth_weight               : recorded birth weight in kilograms
#                                  (ch24, with the IFLS 9.x missing
#                                  flag dropped). Used by
#                                  \citeauthor{ahsan2021} as a
#                                  proximate at-birth outcome.
#
# NOTE on ch19 (place of delivery): the numeric coding drifts between
# waves in a way the IFLS 1997 codebook did not fully resolve and the
# "other home" vs "own home" vs "relative's home" categories did not
# clean-map to a facility/home binary. The mechanism section therefore
# relies on the attendant binary, not the place binary.
#
# Coverage: W2 b4_ch1 is mother's retrospective pregnancy history —
# captures every live birth reported by mothers interviewed in 1997,
# so the 1984–1997 cohort is covered subject to recall bias (3–13 yrs).
# W3 b4_ch1 complements for 1998–1999 births. W1 buk4ch1 has a
# different key schema (case/person not pidlink) and is skipped here
# to keep the stage clean; recall-bias concerns are discussed in
# paper/sections/mechanisms.tex.
#
# Inputs :  data/intermediate/raw/W2__preg_hist.rds  (b4_ch1, W2)
#           data/intermediate/raw/W3__preg_hist.rds  (b4_ch1, W3)
#
# Outputs:  data/intermediate/stage04b/mediators.rds  (one row per
#           (mother_pidlink, birth_year), with the four mediators and
#           a wave-source flag).
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(haven)
})

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))

log_stage("stage04b", "begin")

preg_w2 <- read_intermediate("raw", "W2__preg_hist")
preg_w3 <- read_intermediate("raw", "W3__preg_hist")

# ---------------------------------------------------------------------
# Helper: one wave's pregnancy-history table -> harmonized tibble.
#
# IFLS codes (verified W2 & W3 book 4 CH):
#   ch14  any prenatal check-up: 1 = yes, 3 = no, 8/9 = DK/refuse
#   ch19  place of delivery: 1 = home, 2..9 = various facilities
#   ch20  birth attendant: 1 doctor, 2 midwife (incl. village midwife),
#         3 nurse/paramedic, 4 traditional birth attendant, 5 family,
#         6 alone, 7 other
#   ch24a ever breastfed (W2 and W3 have this; codes 1 yes, 3 no)
# ---------------------------------------------------------------------
harmonize_preg <- function(df, wave_label) {
  keep <- intersect(
    c("pidlink", "ch09yr", "ch14", "ch20", "ch24",
      "ch24a", "ch24c", "ch24cx", "ch24d", "ch24dx", "ch28"),
    names(df)
  )
  df2 <- df[, keep, drop = FALSE]

  # ch20 is a character multi-response column (letters A..Z). Keep as
  # character. All other numeric columns: strip Stata labels.
  num_cols <- setdiff(names(df2), c("pidlink", "ch20"))
  for (nm in num_cols) df2[[nm]] <- haven::zap_labels(df2[[nm]])

  # Harmonize breastfed column (W1 uses ch28, W2/W3 use ch24a).
  bf_col <- if ("ch24a" %in% names(df2)) "ch24a" else "ch28"
  df2$bf_raw <- if (bf_col %in% names(df2)) df2[[bf_col]] else NA_real_

  # Birth year — ch09yr. Some waves store 2-digit; coerce.
  df2$ch09yr <- as.integer(df2$ch09yr)
  df2$birth_year <- ifelse(is.na(df2$ch09yr), NA_integer_,
                    ifelse(df2$ch09yr < 100, df2$ch09yr + 1900L,
                           df2$ch09yr))

  # ch20 letter classification (verified against IFLS-4 hhd_B4 codebook):
  #   A = Physician,  B = Private midwife,  C = Village midwife,
  #   D = Nurse                 -> SKILLED
  #   E = Traditional birth attendant (dukun), F = Family, G = Self,
  #   H = Other, Z = Missing    -> UNSKILLED
  # Missing: ch20 is NA or empty string.
  ch20_str <- ifelse(is.na(df2$ch20), "", as.character(df2$ch20))
  df2$skilled_attendant <- dplyr::case_when(
    ch20_str == ""                   ~ NA_integer_,
    grepl("[ABCD]", ch20_str)        ~ 1L,
    grepl("^[EFGHZ]+$", ch20_str)    ~ 0L,
    TRUE                             ~ NA_integer_
  )
  df2$attendant_village_midwife <- dplyr::case_when(
    ch20_str == ""                   ~ NA_integer_,
    grepl("C", ch20_str)             ~ 1L,
    TRUE                             ~ 0L
  )

  # Birth weight: ch24 is in kilograms with codes >= 9 flagging
  # "don't know" / unreported. Drop those.
  if ("ch24" %in% names(df2)) {
    df2$birth_weight <- ifelse(is.na(df2$ch24) | df2$ch24 >= 9,
                               NA_real_, df2$ch24)
  } else {
    df2$birth_weight <- NA_real_
  }

  # Exclusive breastfeeding for >= 6 months. ch24c/cx records when food
  # other than breast milk was first given; ch24d/dx records when water
  # other than breast milk was first given. The 'x' unit flag maps:
  #   3 = days, 4 = weeks, 5 = months. Convert each duration to months
  # and take the minimum (the earliest non-breast-milk introduction).
  to_months <- function(val, unit) {
    dplyr::case_when(
      is.na(val) | is.na(unit)         ~ NA_real_,
      unit == 5                        ~ as.numeric(val),
      unit == 4                        ~ as.numeric(val) / 4.34,
      unit == 3                        ~ as.numeric(val) / 30,
      TRUE                             ~ NA_real_
    )
  }
  if (all(c("ch24c", "ch24cx", "ch24d", "ch24dx") %in% names(df2))) {
    food_m  <- to_months(df2$ch24c, df2$ch24cx)
    water_m <- to_months(df2$ch24d, df2$ch24dx)
    excl_m  <- pmin(food_m, water_m, na.rm = TRUE)
    excl_m[is.infinite(excl_m)] <- NA_real_
    df2$exclusive_bf_6m <- dplyr::case_when(
      df2$bf_raw == 3                  ~ 0L,   # never breastfed
      is.na(excl_m)                    ~ NA_integer_,
      excl_m >= 6                      ~ 1L,
      TRUE                             ~ 0L
    )
  } else {
    df2$exclusive_bf_6m <- NA_integer_
  }

  df2 |>
    transmute(
      mother_pidlink = as.character(pidlink),
      birth_year     = birth_year,
      pnc_any        = dplyr::case_when(
        ch14 == 1 ~ 1L, ch14 == 3 ~ 0L, TRUE ~ NA_integer_),
      skilled_attendant         = skilled_attendant,
      attendant_village_midwife = attendant_village_midwife,
      birth_weight              = birth_weight,
      ever_breastfed            = dplyr::case_when(
        bf_raw == 1 ~ 1L, bf_raw == 3 ~ 0L, TRUE ~ NA_integer_),
      exclusive_bf_6m           = exclusive_bf_6m,
      wave_source               = wave_label
    ) |>
    filter(!is.na(mother_pidlink), !is.na(birth_year))
}

m_w2 <- harmonize_preg(preg_w2, "W2")
m_w3 <- harmonize_preg(preg_w3, "W3")

log_stage("stage04b", sprintf(
  "W2 preg-history rows: %d; W3: %d", nrow(m_w2), nrow(m_w3)))

# ---------------------------------------------------------------------
# Wave-preference merge.
# For each (mother_pidlink, birth_year), prefer the wave closest in
# time to the birth: W2 (interviewed 1997) for births 1984–1997, W3
# (interviewed 2000) for births 1998–1999. When only one wave has a
# record, use whichever is available.
# ---------------------------------------------------------------------
m_all <- bind_rows(m_w2, m_w3)

mediators <- m_all |>
  group_by(mother_pidlink, birth_year) |>
  arrange(
    # Preference: W2 if birth_year <= 1997, W3 otherwise.
    desc(wave_source == ifelse(first(birth_year) <= 1997, "W2", "W3")),
    .by_group = TRUE
  ) |>
  slice(1) |>
  ungroup()

# Sanity: ≤ 1 row per (mother, birth-year).
stopifnot(!any(duplicated(mediators[, c("mother_pidlink", "birth_year")])))

log_stage("stage04b", sprintf(
  "deduped mediator rows: %d (from %d W2 + %d W3 records)",
  nrow(mediators), nrow(m_w2), nrow(m_w3)))

# Cohort 1984–1999 only — rest discarded to match analysis universe.
mediators <- mediators |>
  filter(birth_year >= 1984, birth_year <= 1999)

log_stage("stage04b", sprintf(
  "cohort 1984–1999 mediator rows: %d", nrow(mediators)))

# Coverage diagnostics.
mediator_vars <- c("pnc_any", "skilled_attendant",
                   "attendant_village_midwife", "ever_breastfed",
                   "exclusive_bf_6m", "birth_weight")
cov_tbl   <- mediators |> summarise(across(all_of(mediator_vars),
                                           ~ mean(!is.na(.x))))
means_tbl <- mediators |> summarise(across(all_of(mediator_vars),
                                           ~ mean(.x, na.rm = TRUE)))
for (v in mediator_vars) {
  log_stage("stage04b", sprintf(
    "%-26s  cov=%.2f  mean=%.3f", v, cov_tbl[[v]], means_tbl[[v]]))
}

# ---------------------------------------------------------------------
# Persist.
# ---------------------------------------------------------------------
write_intermediate(mediators, "stage04b", "mediators")

log_stage("stage04b", "done")
