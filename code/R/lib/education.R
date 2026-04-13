# =====================================================================
# lib/education.R — IFLS education code → years-of-schooling decoder
#
# IFLS `ar16` (highest level attained) and `ar17` (highest grade completed)
# follow a multi-wave rubric. Primary codes (SD), junior/senior high
# (SMP/SMA/SMK), Islamic Madrasah schools (Ibtidaiyah/Tsanawiyah/Aliyah),
# diplomas (D1/D2/D3), and university (S1+) all appear. This helper
# harmonizes them into a single `years` integer.
#
# Audit note (for the upgrade): The legacy project's Maternal_edu.Rmd
# used a slightly different rubric in its `Birth_info copy 3.Rmd` draft
# vs. the later Maternal_edu.Rmd. The version below follows the more
# complete `Birth_info copy 3.Rmd` rubric (it distinguishes D1 vs D2 by
# grade, D3 tiers, and university progress). Islamic-school codes
# (11, 12, 13, 14, 15) are mapped to the same years as their secular
# equivalents (SD, SMP, SMA, etc.) per IFLS documentation.
# =====================================================================

#' Convert IFLS (level, grade) education codes to years of schooling.
#'
#' @param lvl Integer vector, IFLS `ar16` codes. 90-99 → missing.
#' @param grd Integer vector, IFLS `ar17` codes. 90-99 → missing.
#' @return Integer vector of years (0..17), NA where not decodable.
edu_years <- function(lvl, grd) {
  stopifnot(length(lvl) == length(grd))
  lvl <- ifelse(lvl %in% 90:99, NA_real_, as.numeric(lvl))
  grd <- ifelse(grd %in% 90:99, NA_real_, as.numeric(grd))

  dplyr::case_when(
    # No school / kindergarten (secular code 1; W3+ Islamic pre-primary 11)
    lvl %in% c(1, 11) ~ 0,

    # Primary (SD, Madrasah Ibtidaiyah): levels 2, 72 (W3+ extension code)
    # (Islamic primary code 12 handled as MI → SD equivalent.)
    lvl %in% c(2, 72) & !is.na(grd) & grd <= 6 ~ grd,
    lvl %in% c(2, 72) & grd == 7              ~ 6,
    lvl %in% c(2, 72)                         ~ 6,
    lvl == 12 & !is.na(grd) & grd <= 6        ~ grd,
    lvl == 12                                 ~ 6,

    # Junior High (SMP, MTs): levels 3, 4, 73
    lvl %in% c(3, 4, 73) & !is.na(grd) & grd <= 3 ~ 6 + grd,
    lvl %in% c(3, 4, 73) & grd == 7               ~ 9,
    lvl %in% c(3, 4, 73)                          ~ 9,

    # Senior High (SMA/SMK, MA): levels 5, 6, 15, 74
    lvl %in% c(5, 6, 15, 74) & !is.na(grd) & grd <= 3 ~ 9 + grd,
    lvl %in% c(5, 6, 15, 74) & grd == 7               ~ 12,
    lvl %in% c(5, 6, 15, 74)                          ~ 12,

    # Diploma D1/D2 (level 7) — grade distinguishes 13 vs 14
    lvl == 7 & grd %in% 1:2 ~ 12 + grd,
    lvl == 7 & grd == 7     ~ 14,
    lvl == 7                ~ 13,

    # Diploma D3 (level 8)
    lvl == 8 & grd %in% 1:3 ~ 12 + grd,
    lvl == 8 & grd == 7     ~ 15,
    lvl == 8                ~ 15,

    # University (S1, level 9) — progress via grade, graduation at 16
    lvl == 9 & grd %in% 1:4 ~ 12 + grd,
    lvl == 9 & grd == 7     ~ 16,
    lvl == 9                ~ 16,

    # Post-graduate / later-wave tertiary levels (13, 60-63)
    lvl %in% c(13, 60, 61) ~ 16,  # M-equivalent / S2 entry
    lvl %in% c(62, 63)     ~ 18,  # S2 and beyond

    # "Other" (10) → unknown
    lvl == 10 ~ NA_real_,

    # Level missing but grade looks like primary
    is.na(lvl) & !is.na(grd) & grd <= 6 ~ grd,
    is.na(lvl) & grd == 7               ~ 6,

    TRUE ~ NA_real_
  )
}
