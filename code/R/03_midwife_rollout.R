# =====================================================================
# 03_midwife_rollout.R — Village-midwife rollout year per community
#
# Produces THREE intermediates:
#   * midwife_rollout.rds     — commid93-level start_year (SAR + PKK)
#   * rollout_comparison.rds  — SAR vs PKK diagnostic per commid93
#   * other_facilities.rds    — commid93-level presence + start years of
#                               doctor / clinic / private-midwife /
#                               traditional-midwife facilities, for use
#                               as regression-stage confounder controls.
#
# Two data-source channels:
#   (1) SAR (Service Availability Roster) — facility-level enumeration.
#       Modern schema (W4, W5): x14b1 == 8 for village midwife, x12yr
#       for open year, x12b for years-in-operation, x07 for still-open,
#       x08yr / x08b for closure.  Old schema (W2): factpint == 8,
#       calyropn for open year, calyrend for close, operates for status.
#       W3 SAR has the old schema but legacy code did not use it (W3's
#       final `calyropn` coverage was insufficient). We follow legacy:
#       SAR channels come from W2, W4, W5.
#   (2) PKK (Village Administration) — community self-report.
#       Columns identical across W2-W5: j10a (midwife present), j10c
#       (year first arrived), j10d (years in practice). All four waves
#       used.
#
# Both channels are mapped to commid93 anchors via htrack (W5 htrack
# carries commid93/97/00/07/14 jointly). Per commid93 the final
# start_year is min() across waves, clipped to [1989, survey_year].
#
# Legacy note (`midwife.Rmd`): the author observed the SAR/PKK numbers
# "cannot match" and deferred the decision. We surface the mismatch as
# a first-class diagnostic (rollout_comparison) so the regression plan
# can pick one with evidence.
# =====================================================================

set.seed(20260412)

library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))

log_stage("stage03", "begin")

# ---------------------------------------------------------------------
# 1. Load cached SAR/PKK + htrack for commid mapping.
# ---------------------------------------------------------------------
sar_97   <- read_intermediate("raw", "W2__sar")
sar_07   <- read_intermediate("raw", "W4__sar")
sar_14   <- read_intermediate("raw", "W5__sar")

pkk_97   <- read_intermediate("raw", "W2__pkk")
pkk_00   <- read_intermediate("raw", "W3__pkk")
pkk_07   <- read_intermediate("raw", "W4__pkk")
pkk_14   <- read_intermediate("raw", "W5__pkk")

htrack   <- read_intermediate("raw", "W5__htrack")

# ---------------------------------------------------------------------
# 2. commidNN -> commid93 mappings (many-to-one by construction).
# ---------------------------------------------------------------------
make_map <- function(col_wave) {
  htrack |>
    select(all_of(c(col_wave, "commid93"))) |>
    filter(!is.na(.data[[col_wave]]), !is.na(commid93)) |>
    distinct() |>
    # Collapse ambiguous cases (rare): pick the modal commid93 per wave commid.
    group_by(.data[[col_wave]]) |>
    slice(1) |>
    ungroup()
}

map_97_to_93 <- make_map("commid97")
map_00_to_93 <- make_map("commid00")
map_07_to_93 <- make_map("commid07")
map_14_to_93 <- make_map("commid14")

log_stage("stage03", sprintf(
  "commid mappings — 97:%d 00:%d 07:%d 14:%d",
  nrow(map_97_to_93), nrow(map_00_to_93),
  nrow(map_07_to_93), nrow(map_14_to_93)))

# ---------------------------------------------------------------------
# 3. SAR: per-wave village-midwife rollout years.
# ---------------------------------------------------------------------

#' Modern-schema SAR (W4, W5). Expects columns x14b1, x12yr, x12b, x07,
#' x08yr, x08b plus a `commid_col` with the wave-specific commid.
sar_modern <- function(sar, commid_col, survey_year) {
  sar |>
    filter(x14b1 == 8, !is.na(.data[[commid_col]])) |>
    mutate(
      open_year_fac = coalesce(
        as.integer(x12yr),
        if_else(!is.na(x12b), survey_year - as.integer(x12b), NA_integer_)
      ),
      close_year_fac = case_when(
        x07 == 1      ~ NA_integer_,
        !is.na(x08yr) ~ as.integer(x08yr),
        !is.na(x08b)  ~ survey_year - as.integer(x08b),
        TRUE          ~ NA_integer_
      ),
      open_year_fac = if_else(
        open_year_fac >= 1989 & open_year_fac <= survey_year,
        open_year_fac, NA_integer_
      )
    ) |>
    group_by(.data[[commid_col]]) |>
    summarise(
      start_sar = suppressWarnings(min(open_year_fac, na.rm = TRUE)),
      present_sar = as.integer(any(
        x07 == 1 |
        (!is.na(open_year_fac) & open_year_fac <= survey_year &
           (is.na(close_year_fac) | close_year_fac >= survey_year))
      )),
      n_fac = n(),
      .groups = "drop"
    ) |>
    mutate(start_sar = if_else(is.infinite(start_sar), NA_real_, start_sar))
}

#' Old-schema SAR (W2). Uses factpint, calyropn, calyrend, operates.
sar_old_97 <- function(sar, survey_year = 1997L) {
  sar |>
    filter(factpint == 8, !is.na(commid97)) |>
    mutate(
      open_year_fac  = as.integer(calyropn),
      close_year_fac = as.integer(calyrend),
      open_year_fac  = if_else(
        open_year_fac >= 1989 & open_year_fac <= survey_year,
        open_year_fac, NA_integer_
      )
    ) |>
    group_by(commid97) |>
    summarise(
      start_sar = suppressWarnings(min(open_year_fac, na.rm = TRUE)),
      present_sar = as.integer(any(
        operates == 1 |
        (!is.na(open_year_fac) & open_year_fac <= survey_year &
           (is.na(close_year_fac) | close_year_fac >= survey_year))
      )),
      n_fac = n(),
      .groups = "drop"
    ) |>
    mutate(start_sar = if_else(is.infinite(start_sar), NA_real_, start_sar))
}

sar_97_commid <- sar_old_97(sar_97)                           # keyed on commid97
sar_07_commid <- sar_modern(sar_07, "commid07", 2007L)        # keyed on commid07
sar_14_commid <- sar_modern(sar_14, "commid14", 2015L)        # keyed on commid14

# Anchor to commid93.
sar_97_c93 <- sar_97_commid |>
  left_join(map_97_to_93, by = "commid97") |>
  filter(!is.na(commid93)) |>
  group_by(commid93) |>
  summarise(start_sar_97 = suppressWarnings(min(start_sar, na.rm = TRUE)),
            present_sar_97 = suppressWarnings(max(present_sar, na.rm = TRUE)),
            .groups = "drop") |>
  mutate(start_sar_97 = if_else(is.infinite(start_sar_97),
                                NA_real_, start_sar_97))

sar_07_c93 <- sar_07_commid |>
  left_join(map_07_to_93, by = "commid07") |>
  filter(!is.na(commid93)) |>
  group_by(commid93) |>
  summarise(start_sar_07 = suppressWarnings(min(start_sar, na.rm = TRUE)),
            present_sar_07 = suppressWarnings(max(present_sar, na.rm = TRUE)),
            .groups = "drop") |>
  mutate(start_sar_07 = if_else(is.infinite(start_sar_07),
                                NA_real_, start_sar_07))

sar_14_c93 <- sar_14_commid |>
  left_join(map_14_to_93, by = "commid14") |>
  filter(!is.na(commid93)) |>
  group_by(commid93) |>
  summarise(start_sar_14 = suppressWarnings(min(start_sar, na.rm = TRUE)),
            present_sar_14 = suppressWarnings(max(present_sar, na.rm = TRUE)),
            .groups = "drop") |>
  mutate(start_sar_14 = if_else(is.infinite(start_sar_14),
                                NA_real_, start_sar_14))

# ---------------------------------------------------------------------
# 4. PKK: per-wave village-midwife starting year (j10c when j10a == 1).
# ---------------------------------------------------------------------
pkk_wave <- function(pkk, commid_col, survey_year) {
  pkk |>
    filter(!is.na(.data[[commid_col]])) |>
    mutate(
      start_pkk = if_else(!is.na(j10a) & j10a == 1 &
                            as.integer(j10c) >= 1989 &
                            as.integer(j10c) <= survey_year,
                          as.numeric(j10c), NA_real_),
      midwife_present = as.integer(!is.na(j10a) & j10a == 1)
    ) |>
    select(all_of(commid_col), start_pkk, midwife_present)
}

pkk_97_c <- pkk_wave(pkk_97, "commid97", 1997L)
pkk_00_c <- pkk_wave(pkk_00, "commid00", 2000L)
pkk_07_c <- pkk_wave(pkk_07, "commid07", 2007L)
pkk_14_c <- pkk_wave(pkk_14, "commid14", 2015L)

anchor_pkk <- function(df, commid_col, map, new_start_col, new_pres_col) {
  df |>
    left_join(map, by = commid_col) |>
    filter(!is.na(commid93)) |>
    group_by(commid93) |>
    summarise(
      !!new_start_col := suppressWarnings(min(start_pkk, na.rm = TRUE)),
      !!new_pres_col  := max(midwife_present, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(!!new_start_col := if_else(is.infinite(.data[[new_start_col]]),
                                      NA_real_, .data[[new_start_col]]))
}

pkk_97_c93 <- anchor_pkk(pkk_97_c, "commid97", map_97_to_93,
                         "start_pkk_97", "present_pkk_97")
pkk_00_c93 <- anchor_pkk(pkk_00_c, "commid00", map_00_to_93,
                         "start_pkk_00", "present_pkk_00")
pkk_07_c93 <- anchor_pkk(pkk_07_c, "commid07", map_07_to_93,
                         "start_pkk_07", "present_pkk_07")
pkk_14_c93 <- anchor_pkk(pkk_14_c, "commid14", map_14_to_93,
                         "start_pkk_14", "present_pkk_14")

# ---------------------------------------------------------------------
# 5. Union to commid93 level: start_year_sar = min across SAR waves;
#    start_year_pkk = min across PKK waves. present_* = any-wave presence.
# ---------------------------------------------------------------------
rollout <- list(sar_97_c93, sar_07_c93, sar_14_c93,
                pkk_97_c93, pkk_00_c93, pkk_07_c93, pkk_14_c93) |>
  purrr::reduce(full_join, by = "commid93")

rollout <- rollout |>
  rowwise() |>
  mutate(
    start_year_sar = suppressWarnings(min(
      c(start_sar_97, start_sar_07, start_sar_14), na.rm = TRUE
    )),
    start_year_pkk = suppressWarnings(min(
      c(start_pkk_97, start_pkk_00, start_pkk_07, start_pkk_14), na.rm = TRUE
    ))
  ) |>
  ungroup() |>
  mutate(
    start_year_sar = if_else(is.infinite(start_year_sar),
                             NA_real_, start_year_sar),
    start_year_pkk = if_else(is.infinite(start_year_pkk),
                             NA_real_, start_year_pkk),
    # present_* columns may be -Inf if every wave was NA for that commid;
    # clamp -Inf (unobserved) to NA and treat a real 0 as "observed, no VM".
    across(starts_with("present_"),
           ~ if_else(is.infinite(.x), NA_real_, .x)),
    present_any_sar = suppressWarnings(pmax(present_sar_97, present_sar_07,
                                            present_sar_14, na.rm = TRUE)),
    present_any_pkk = suppressWarnings(pmax(present_pkk_97, present_pkk_00,
                                            present_pkk_07, present_pkk_14,
                                            na.rm = TRUE)),
    present_any_sar = if_else(is.infinite(present_any_sar),
                              NA_real_, present_any_sar),
    present_any_pkk = if_else(is.infinite(present_any_pkk),
                              NA_real_, present_any_pkk)
  )

# ---------------------------------------------------------------------
# 6. SAR vs PKK comparison diagnostic.
# ---------------------------------------------------------------------
comparison <- rollout |>
  transmute(
    commid93, start_year_sar, start_year_pkk,
    delta = start_year_pkk - start_year_sar,
    has_sar = as.integer(!is.na(start_year_sar)),
    has_pkk = as.integer(!is.na(start_year_pkk)),
    has_both = as.integer(!is.na(start_year_sar) & !is.na(start_year_pkk))
  )

both <- comparison |> filter(has_both == 1)
n_commids <- nrow(rollout)
n_sar     <- sum(comparison$has_sar)
n_pkk     <- sum(comparison$has_pkk)
n_both    <- sum(comparison$has_both)
cor_sp    <- if (n_both > 3) cor(both$start_year_sar,
                                 both$start_year_pkk,
                                 use = "complete.obs") else NA_real_
share_disagree <- if (n_both > 0)
  mean(both$start_year_sar != both$start_year_pkk, na.rm = TRUE) else NA_real_

log_stage("stage03", sprintf(
  "communities: total=%d with SAR=%d with PKK=%d with both=%d; rho(SAR,PKK)=%.3f share-disagree=%.3f",
  n_commids, n_sar, n_pkk, n_both, cor_sp, share_disagree))

# ---------------------------------------------------------------------
# 7. Other facilities (modern W4, W5 SAR only). x14b1 codes:
#       1 = puskesmas, 2 = puskesmas pembantu, 3 = polindes,
#       4 = posyandu, 5 = other govt, 6 = private hospital,
#       7 = private midwife, 8 = village midwife, 9 = private doctor,
#       10 = nurse paramedic, 11 = traditional,
#       12 = traditional birth attendant
#    We capture presence of four confounder groups that can spuriously
#    correlate with VM rollout. Anchored to commid93.
# ---------------------------------------------------------------------
other_facility_wave <- function(sar, commid_col, survey_year) {
  sar |>
    filter(!is.na(.data[[commid_col]])) |>
    mutate(
      open_year = coalesce(
        as.integer(x12yr),
        if_else(!is.na(x12b), survey_year - as.integer(x12b), NA_integer_)
      ),
      open_year = if_else(open_year >= 1970 & open_year <= survey_year,
                          open_year, NA_integer_),
      type_doctor              = as.integer(x14b1 %in% c(1, 6, 9)),
      type_clinic              = as.integer(x14b1 %in% c(2, 3, 4, 5)),
      type_private_midwife     = as.integer(x14b1 == 7),
      type_traditional_midwife = as.integer(x14b1 %in% c(11, 12))
    ) |>
    group_by(.data[[commid_col]]) |>
    summarise(
      has_doctor              = max(type_doctor, 0L),
      has_clinic              = max(type_clinic, 0L),
      has_private_midwife     = max(type_private_midwife, 0L),
      has_traditional_midwife = max(type_traditional_midwife, 0L),
      start_doctor = suppressWarnings(min(
        open_year[type_doctor == 1], na.rm = TRUE)),
      start_clinic = suppressWarnings(min(
        open_year[type_clinic == 1], na.rm = TRUE)),
      start_private_midwife = suppressWarnings(min(
        open_year[type_private_midwife == 1], na.rm = TRUE)),
      start_traditional_midwife = suppressWarnings(min(
        open_year[type_traditional_midwife == 1], na.rm = TRUE)),
      .groups = "drop"
    ) |>
    mutate(across(starts_with("start_"),
                  ~ if_else(is.infinite(.x), NA_real_, .x)))
}

other_07 <- other_facility_wave(sar_07, "commid07", 2007L) |>
  left_join(map_07_to_93, by = "commid07") |>
  filter(!is.na(commid93)) |>
  select(-commid07)

other_14 <- other_facility_wave(sar_14, "commid14", 2015L) |>
  left_join(map_14_to_93, by = "commid14") |>
  filter(!is.na(commid93)) |>
  select(-commid14)

other_facilities <- bind_rows(other_07, other_14) |>
  filter(!is.na(commid93), commid93 != "") |>
  group_by(commid93) |>
  summarise(
    has_doctor = suppressWarnings(max(has_doctor, na.rm = TRUE)),
    has_clinic = suppressWarnings(max(has_clinic, na.rm = TRUE)),
    has_private_midwife = suppressWarnings(max(has_private_midwife,
                                               na.rm = TRUE)),
    has_traditional_midwife = suppressWarnings(max(has_traditional_midwife,
                                                   na.rm = TRUE)),
    start_doctor = suppressWarnings(min(start_doctor, na.rm = TRUE)),
    start_clinic = suppressWarnings(min(start_clinic, na.rm = TRUE)),
    start_private_midwife = suppressWarnings(min(
      start_private_midwife, na.rm = TRUE)),
    start_traditional_midwife = suppressWarnings(min(
      start_traditional_midwife, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(
    across(starts_with("start_"),
           ~ if_else(is.infinite(.x), NA_real_, .x)),
    across(starts_with("has_"),
           ~ if_else(is.infinite(.x), NA_real_, .x))
  )

log_stage("stage03", sprintf(
  "other-facility coverage: doctor=%d clinic=%d private-midwife=%d traditional=%d",
  sum(other_facilities$has_doctor == 1, na.rm = TRUE),
  sum(other_facilities$has_clinic == 1, na.rm = TRUE),
  sum(other_facilities$has_private_midwife == 1, na.rm = TRUE),
  sum(other_facilities$has_traditional_midwife == 1, na.rm = TRUE)))

# ---------------------------------------------------------------------
# 8. Persist outputs.
# ---------------------------------------------------------------------
write_intermediate(rollout,          "stage03", "midwife_rollout")
write_intermediate(comparison,       "stage03", "rollout_comparison")
write_intermediate(other_facilities, "stage03", "other_facilities")

log_stage("stage03", "done")
