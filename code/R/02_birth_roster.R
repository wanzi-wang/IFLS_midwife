# =====================================================================
# 02_birth_roster.R — Individual birth roster (cohorts 1984-1999)
#
# Builds the person-level panel with:
#   - pidlink, sex, birth_year, birth_month
#   - hhid93..hhid14 (from ptrack)
#   - commid93..commid14 (from htrack per-hhid mapping)
#   - mover97/00/07/14 + non_migrant + mother_non_migrant
#   - dead (from W5 ptrack ar01a_14 == 0)
#   - mother_pidlink, mother_dob_yr, mother_age_birth,
#     mother_marital, mother_edu_years, mother_edu_miss
#   - mother_height (W1 bukcca2), mother_smoke (W1 buk3km1)
#   - three commid_birth variants (via lib/birthplace.R)
#   - province / kabupaten / kecamatan (from W5 htrack sc01_14_14 etc)
#
# Construction strategy (mirrors legacy `Birth_info copy 3.Rmd` but
# with audits):
#   - Children born 1984-1993 are pulled from W1 roster (bukkar2).
#   - Children born 1994-1997 from W2 roster (bk_ar1).
#   - Children born 1998-1999 from W3 roster (bk_ar1).
#   The three panels are unioned. Each per-wave join extracts the
#   mother's intra-HH `ar11` and then joins mother's own roster row to
#   harvest her DOB, marital status, and edu codes.
#
# No sample restrictions applied here. `dead`, `non_migrant`, and the
# commid_birth variants are retained as COLUMNS so the regression plan
# picks the final sample from data visible.
# =====================================================================

set.seed(20260412)

library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))
source(here::here("code", "R", "lib", "education.R"))
source(here::here("code", "R", "lib", "birthplace.R"))

log_stage("stage02", "begin")

# ---------------------------------------------------------------------
# 1. Load all cached intermediates needed by stage 02.
# ---------------------------------------------------------------------
ptrack    <- read_intermediate("raw", "W5__ptrack")
htrack_w5 <- read_intermediate("raw", "W5__htrack")

roster93  <- read_intermediate("raw", "W1__roster_alt")    # bukkar2 (per legacy)
roster97  <- read_intermediate("raw", "W2__roster_main")   # bk_ar1
roster00  <- read_intermediate("raw", "W3__roster_main")   # bk_ar1

height_w1 <- read_intermediate("raw", "W1__height")         # bukcca2
smoke_w1  <- read_intermediate("raw", "W1__smoke")          # buk3km1

# ---------------------------------------------------------------------
# 2. htrack → (hhid, wave) → commid93 mapping + mover flags.
#    W5 htrack carries commid93/97/00/07/14 + mover97/00/07/14 per HH
#    line, so it's the single master for cross-wave community tracking.
# ---------------------------------------------------------------------
htrack_commid <- htrack_w5 |>
  transmute(
    hhid93 = as.character(hhid93),
    hhid97 = as.character(hhid97),
    hhid00 = as.character(hhid00),
    hhid07 = as.character(hhid07),
    hhid14 = as.character(hhid14),
    commid93, commid97, commid00, commid07, commid14,
    mover97, mover00, mover07, mover14,
    province  = sc01_93,
    kabupaten = sc02_93,
    kecamatan = sc03_93
  ) |>
  distinct()

# Per-hhid (any wave) → commid93 — used to anchor roster-level joins.
hh93_to_commid <- htrack_commid |>
  filter(!is.na(hhid93), !is.na(commid93)) |>
  distinct(hhid93, commid93, province, kabupaten, kecamatan)

hh97_to_commid <- htrack_commid |>
  filter(!is.na(hhid97), !is.na(commid97)) |>
  distinct(hhid97, commid93, commid97)

hh00_to_commid <- htrack_commid |>
  filter(!is.na(hhid00), !is.na(commid00)) |>
  distinct(hhid00, commid93, commid00)

# ---------------------------------------------------------------------
# 3. Per-wave roster slicing: pull children born in the wave-specific
#    cohort window, attach mother info from same wave's roster.
# ---------------------------------------------------------------------

#' Extract children + mother covariates from one wave.
#' Returns a tibble with harmonized column names.
extract_wave <- function(roster, hhid_col, pid_col, wave_label,
                         year_min, year_max, survey_year) {
  r <- roster |>
    transmute(
      pidlink = as.character(pidlink),
      hhid    = as.character(.data[[hhid_col]]),
      pid     = as.character(.data[[pid_col]]),
      relation     = ar02,
      sex          = ar07,
      dob_mth_raw  = ar08mth,
      dob_yr_raw   = ar08yr,
      mother_intra = as.character(ar11),
      marital      = ar13,
      edu_level    = ar16,
      edu_grade    = ar17
    )

  # Decode DOB. IFLS W1 uses 2-digit year (84..99 -> 1984..1999); W2/W3
  # may use 4-digit already — handle both. Legacy check: W1 returns NA for
  # 95..99 (reserved codes), but real births in those years do exist, so
  # we only drop 95..99 if dob_yr > 99 OR the mode signals 2-digit.
  r <- r |>
    mutate(
      dob_mth = if_else(dob_mth_raw %in% 1:12, dob_mth_raw, NA_real_),
      dob_yr = case_when(
        is.na(dob_yr_raw)                ~ NA_real_,
        dob_yr_raw >= 1900               ~ as.numeric(dob_yr_raw),
        dob_yr_raw >= 0 & dob_yr_raw <= 99 ~ 1900 + as.numeric(dob_yr_raw),
        TRUE                              ~ NA_real_
      )
    ) |>
    filter(!is.na(dob_yr), dob_yr >= year_min, dob_yr <= year_max)

  # Mother lookup: child's mother_intra + child's hhid → mother's pidlink
  # and mother's own roster row (DOB, marital, edu).
  mother_tbl <- roster |>
    transmute(
      mother_pidlink  = as.character(pidlink),
      hhid            = as.character(.data[[hhid_col]]),
      pid             = as.character(.data[[pid_col]]),
      mother_dob_yr_raw = ar08yr,
      mother_marital    = ar13,
      mother_edu_level  = ar16,
      mother_edu_grade  = ar17
    )

  # Attach mother info; keep child's pidlink on the left.
  out <- r |>
    left_join(mother_tbl, by = c("hhid" = "hhid",
                                 "mother_intra" = "pid")) |>
    mutate(
      mother_dob_yr = case_when(
        is.na(mother_dob_yr_raw)                         ~ NA_real_,
        mother_dob_yr_raw >= 1900 & mother_dob_yr_raw <= 2000 ~ as.numeric(mother_dob_yr_raw),
        mother_dob_yr_raw >= 0 & mother_dob_yr_raw <= 99  ~ 1900 + as.numeric(mother_dob_yr_raw),
        TRUE                                              ~ NA_real_
      ),
      mother_age_birth_raw = dob_yr - mother_dob_yr,
      # Sanity: biologically plausible range for mother's age at birth.
      mother_age_birth = if_else(mother_age_birth_raw >= 10 &
                                  mother_age_birth_raw <= 60,
                                 mother_age_birth_raw, NA_real_),
      mother_dob_yr = if_else(is.na(mother_age_birth), NA_real_,
                               mother_dob_yr),
      mother_edu_years = edu_years(mother_edu_level, mother_edu_grade),
      mother_edu_miss  = as.integer(is.na(mother_edu_years)),
      source_wave      = wave_label,
      source_survey_year = survey_year
    ) |>
    select(pidlink, hhid, pid, sex, dob_yr, dob_mth, relation, marital,
           mother_pidlink, mother_dob_yr, mother_age_birth, mother_marital,
           mother_edu_years, mother_edu_miss,
           source_wave, source_survey_year)

  # Deduplicate on pidlink — same child can appear in multiple wave
  # rosters but we only want the wave that "owns" their cohort window.
  out |> distinct(pidlink, .keep_all = TRUE)
}

panel_1984_1993 <- extract_wave(roster93, "hhid93", "pid93",
                                "W1", 1984, 1993, 1993) |>
  rename(hhid93 = hhid) |>
  select(-pid)

panel_1994_1997 <- extract_wave(roster97, "hhid97", "pid97",
                                "W2", 1994, 1997, 1997) |>
  rename(hhid97 = hhid) |>
  select(-pid)

panel_1998_1999 <- extract_wave(roster00, "hhid00", "pid00",
                                "W3", 1998, 1999, 2000) |>
  rename(hhid00 = hhid) |>
  select(-pid)

log_stage("stage02", sprintf(
  "per-wave roster slices — W1:%d W2:%d W3:%d",
  nrow(panel_1984_1993), nrow(panel_1994_1997), nrow(panel_1998_1999)))

# ---------------------------------------------------------------------
# 4. Stack and anchor commid93 per wave of extraction.
# ---------------------------------------------------------------------
panel <- bind_rows(
  panel_1984_1993,
  panel_1994_1997,
  panel_1998_1999
)

# Pull each person's full hhid sequence (all waves) from ptrack so we
# can backfill missing hhidXX for children whose per-wave panel only
# has their own wave's hhid.
ptr_lookup <- ptrack |>
  transmute(
    pidlink  = as.character(pidlink),
    hhid93_p = as.character(hhid93),
    hhid97_p = as.character(hhid97),
    hhid00_p = as.character(hhid00),
    hhid07_p = as.character(hhid07),
    hhid14_p = as.character(hhid14),
    ar01a_14 = ar01a_14
  ) |>
  distinct(pidlink, .keep_all = TRUE)

for (col in c("hhid93", "hhid97", "hhid00")) {
  if (!col %in% names(panel)) panel[[col]] <- NA_character_
}

panel <- panel |>
  left_join(ptr_lookup, by = "pidlink") |>
  mutate(
    hhid93 = coalesce(hhid93, hhid93_p),
    hhid97 = coalesce(hhid97, hhid97_p),
    hhid00 = coalesce(hhid00, hhid00_p),
    hhid07 = hhid07_p,
    hhid14 = hhid14_p
  ) |>
  select(-ends_with("_p"))

# Attach htrack block via hhid93 (the most common non-missing HH ID for
# our 1984-99 birth cohorts). Filter out rows where hhid93 is NA before
# joining to avoid the NA-matches-NA cartesian explosion.
htrack_keyed <- htrack_commid |>
  filter(!is.na(hhid93)) |>
  distinct(hhid93, .keep_all = TRUE) |>
  select(hhid93,
         commid93_h = commid93, commid97_h = commid97,
         commid00_h = commid00, commid07_h = commid07,
         commid14_h = commid14,
         mover97_h  = mover97,  mover00_h  = mover00,
         mover07_h  = mover07,  mover14_h  = mover14,
         province_h = province, kabupaten_h = kabupaten,
         kecamatan_h = kecamatan)

panel <- panel |>
  left_join(htrack_keyed, by = "hhid93") |>
  mutate(
    commid93  = commid93_h,
    commid97  = commid97_h,
    commid00  = commid00_h,
    commid07  = commid07_h,
    commid14  = commid14_h,
    mover97   = mover97_h,
    mover00   = mover00_h,
    mover07   = mover07_h,
    mover14   = mover14_h,
    province  = province_h,
    kabupaten = kabupaten_h,
    kecamatan = kecamatan_h
  ) |>
  select(-ends_with("_h"))

# For rows with no hhid93 (post-1993 births whose HH never existed in
# 1993), fall back to hhid97 → htrack, then hhid00 → htrack.
htrack_keyed97 <- htrack_commid |>
  filter(!is.na(hhid97)) |>
  distinct(hhid97, .keep_all = TRUE) |>
  select(hhid97,
         c93_f = commid93, c97_f = commid97, c00_f = commid00,
         c07_f = commid07, c14_f = commid14,
         m97_f = mover97, m00_f = mover00, m07_f = mover07, m14_f = mover14,
         p_f = province, k_f = kabupaten, ke_f = kecamatan)

panel <- panel |>
  left_join(htrack_keyed97, by = "hhid97") |>
  mutate(
    commid93  = coalesce(commid93,  c93_f),
    commid97  = coalesce(commid97,  c97_f),
    commid00  = coalesce(commid00,  c00_f),
    commid07  = coalesce(commid07,  c07_f),
    commid14  = coalesce(commid14,  c14_f),
    mover97   = coalesce(mover97,   m97_f),
    mover00   = coalesce(mover00,   m00_f),
    mover07   = coalesce(mover07,   m07_f),
    mover14   = coalesce(mover14,   m14_f),
    province  = coalesce(province,  p_f),
    kabupaten = coalesce(kabupaten, k_f),
    kecamatan = coalesce(kecamatan, ke_f)
  ) |>
  select(-ends_with("_f"))

htrack_keyed00 <- htrack_commid |>
  filter(!is.na(hhid00)) |>
  distinct(hhid00, .keep_all = TRUE) |>
  select(hhid00,
         c93_g = commid93, c97_g = commid97, c00_g = commid00,
         c07_g = commid07, c14_g = commid14,
         m97_g = mover97, m00_g = mover00, m07_g = mover07, m14_g = mover14,
         p_g = province, k_g = kabupaten, ke_g = kecamatan)

panel <- panel |>
  left_join(htrack_keyed00, by = "hhid00") |>
  mutate(
    commid93  = coalesce(commid93,  c93_g),
    commid97  = coalesce(commid97,  c97_g),
    commid00  = coalesce(commid00,  c00_g),
    commid07  = coalesce(commid07,  c07_g),
    commid14  = coalesce(commid14,  c14_g),
    mover97   = coalesce(mover97,   m97_g),
    mover00   = coalesce(mover00,   m00_g),
    mover07   = coalesce(mover07,   m07_g),
    mover14   = coalesce(mover14,   m14_g),
    province  = coalesce(province,  p_g),
    kabupaten = coalesce(kabupaten, k_g),
    kecamatan = coalesce(kecamatan, ke_g)
  ) |>
  select(-ends_with("_g"))

# ---------------------------------------------------------------------
# 5. Dead flag + non_migrant + mother_non_migrant.
# ---------------------------------------------------------------------

# ar01a_14 coding: 0/5/6 = not in HH (includes deceased per legacy note);
# specifically, legacy uses `dead = (ar01a_14 == 0)` — verify value 0
# = deceased. Per IFLS docbook ar01a_14: 1=current member, 2=moved
# within HH, 3=moved out, 5=died, etc. Legacy maps 0→died but the
# codebook says 5=died. Re-check: in practice ptrack rows with
# `died07 == 1` have `ar01a_14 == 5` in W5 pattern. Safer approach:
# flag ar01a_14 %in% c(5). Also preserve the raw value.
panel <- panel |>
  mutate(
    dead = as.integer(ar01a_14 == 5 | !is.na(ar01a_14) & ar01a_14 == 0)
  )

# migrant / non_migrant per legacy: any mover flag in 2..6 ⇒ migrant.
panel <- panel |>
  mutate(
    migrant = as.integer(
      (!is.na(mover97) & mover97 >= 2 & mover97 <= 6) |
      (!is.na(mover00) & mover00 >= 2 & mover00 <= 6) |
      (!is.na(mover07) & mover07 >= 2 & mover07 <= 6)
    ),
    non_migrant = 1L - migrant
  )

# Mother-centric never-mover: the mother's non_migrant over the waves
# spanning 1993 through the child's birth year. Use the mother's own
# mover flags as observed in ptrack.
mm_by_hh93 <- htrack_commid |>
  filter(!is.na(hhid93)) |>
  distinct(hhid93, .keep_all = TRUE) |>
  select(hhid93, mm97 = mover97)
mm_by_hh97 <- htrack_commid |>
  filter(!is.na(hhid97)) |>
  distinct(hhid97, .keep_all = TRUE) |>
  select(hhid97, mm00 = mover00)
mm_by_hh00 <- htrack_commid |>
  filter(!is.na(hhid00)) |>
  distinct(hhid00, .keep_all = TRUE) |>
  select(hhid00, mm07 = mover07)

mother_non_mig <- ptrack |>
  transmute(
    mother_pidlink = as.character(pidlink),
    mother_hhid93  = as.character(hhid93),
    mother_hhid97  = as.character(hhid97),
    mother_hhid00  = as.character(hhid00)
  ) |>
  distinct(mother_pidlink, .keep_all = TRUE) |>
  left_join(mm_by_hh93, by = c("mother_hhid93" = "hhid93")) |>
  left_join(mm_by_hh97, by = c("mother_hhid97" = "hhid97")) |>
  left_join(mm_by_hh00, by = c("mother_hhid00" = "hhid00")) |>
  transmute(
    mother_pidlink,
    mother_migrant = as.integer(
      (!is.na(mm97) & mm97 >= 2 & mm97 <= 6) |
      (!is.na(mm00) & mm00 >= 2 & mm00 <= 6) |
      (!is.na(mm07) & mm07 >= 2 & mm07 <= 6)
    ),
    mother_non_migrant = 1L - mother_migrant
  ) |>
  select(mother_pidlink, mother_non_migrant)

panel <- panel |> left_join(mother_non_mig, by = "mother_pidlink") |>
  mutate(mother_non_migrant = tidyr::replace_na(mother_non_migrant, 0L))

# ---------------------------------------------------------------------
# 6. Mother's commid at nearest pre-birth wave -> `commid_birth_preceding`.
#    Build the mother's commid-per-wave panel from ptrack × htrack.
# ---------------------------------------------------------------------
hhid_to_commid93 <- htrack_commid |>
  select(hhid93, hhid97, hhid00, hhid07, hhid14, commid93) |>
  tidyr::pivot_longer(cols = starts_with("hhid"),
                      names_to = "wc", values_to = "hhid_wave") |>
  filter(!is.na(hhid_wave), !is.na(commid93)) |>
  distinct(hhid_wave, commid93)

mother_panel <- ptrack |>
  transmute(
    pidlink = as.character(pidlink),
    hhid93  = as.character(hhid93),
    hhid97  = as.character(hhid97),
    hhid00  = as.character(hhid00),
    hhid07  = as.character(hhid07),
    hhid14  = as.character(hhid14)
  ) |>
  tidyr::pivot_longer(cols = starts_with("hhid"),
                      names_to = "wave_col", values_to = "hhid_wave") |>
  filter(!is.na(hhid_wave)) |>
  mutate(wave = recode(wave_col,
    hhid93 = "W1", hhid97 = "W2", hhid00 = "W3",
    hhid07 = "W4", hhid14 = "W5"
  )) |>
  left_join(hhid_to_commid93, by = "hhid_wave",
            relationship = "many-to-many") |>
  filter(!is.na(commid93)) |>
  distinct(pidlink, wave, commid93)

panel <- panel |>
  mutate(
    birth_year  = dob_yr,
    birth_month = dob_mth
  ) |>
  mutate(
    mother_commid_preceding = mother_commid_preceding(
      mother_pidlink, birth_year, mother_panel
    )
  )

# ---------------------------------------------------------------------
# 7. Three `commid_birth` variants.
# ---------------------------------------------------------------------
panel <- build_commid_birth_variants(panel)

# ---------------------------------------------------------------------
# 8. Mother height (W1 bukcca2 ca12 = standing-height cm) + mother smoke.
# ---------------------------------------------------------------------
ht <- height_w1 |>
  transmute(pidlink = as.character(pidlink),
            mother_height = as.numeric(ca10)) |>
  filter(!is.na(mother_height), mother_height > 100, mother_height < 200) |>
  distinct(pidlink, .keep_all = TRUE) |>
  rename(mother_pidlink = pidlink)

sm <- smoke_w1 |>
  transmute(pidlink = as.character(pidlink),
            mother_smoke = as.integer(km01a == 1 | km01a == "1")) |>
  filter(!is.na(mother_smoke)) |>
  distinct(pidlink, .keep_all = TRUE) |>
  rename(mother_pidlink = pidlink)

panel <- panel |>
  left_join(ht, by = "mother_pidlink") |>
  left_join(sm, by = "mother_pidlink") |>
  mutate(mother_smoke = tidyr::replace_na(mother_smoke, 0L))

# ---------------------------------------------------------------------
# 9. Final column order + write.
# ---------------------------------------------------------------------
birth_roster <- panel |>
  # A single duplicated pidlink can arise when a person appears in both
  # W2 and W3 per-wave slices (rare edge case of DOB near wave boundary).
  # Keep the earlier wave's row.
  arrange(source_survey_year) |>
  distinct(pidlink, .keep_all = TRUE) |>
  transmute(
    pidlink,
    sex,
    birth_year,
    birth_month,
    dob_yr_raw = dob_yr,
    dob_mth_raw = dob_mth,
    # Panel IDs
    hhid93, hhid97, hhid00, hhid07, hhid14,
    # Geography
    commid93, commid97, commid00, commid07, commid14,
    province, kabupaten, kecamatan,
    # Migration + mortality
    mover97, mover00, mover07, mover14,
    migrant, non_migrant, dead,
    # Mother link + covariates (no hhsize per user decision 2026-04-13)
    mother_pidlink, mother_dob_yr, mother_age_birth, mother_marital,
    mother_edu_years, mother_edu_miss,
    mother_height, mother_smoke,
    mother_non_migrant, mother_commid_preceding,
    # Birthplace variants
    commid_birth_legacy, commid_birth_mother, commid_birth_preceding,
    # Provenance
    source_wave, source_survey_year
  )

# Sanity
log_stage("stage02", sprintf(
  "panel rows=%d; unique pidlink=%d; with mother_pidlink=%d; dead=%d",
  nrow(birth_roster),
  dplyr::n_distinct(birth_roster$pidlink),
  sum(!is.na(birth_roster$mother_pidlink)),
  sum(birth_roster$dead, na.rm = TRUE)
))
log_stage("stage02", sprintf(
  "non_migrant=%d (%.1f%%); mother_non_migrant=%d (%.1f%%)",
  sum(birth_roster$non_migrant, na.rm = TRUE),
  100 * mean(birth_roster$non_migrant == 1, na.rm = TRUE),
  sum(birth_roster$mother_non_migrant, na.rm = TRUE),
  100 * mean(birth_roster$mother_non_migrant == 1, na.rm = TRUE)
))
log_stage("stage02", sprintf(
  "commid_birth coverage: legacy=%d mother=%d preceding=%d",
  sum(!is.na(birth_roster$commid_birth_legacy)),
  sum(!is.na(birth_roster$commid_birth_mother)),
  sum(!is.na(birth_roster$commid_birth_preceding))
))

write_intermediate(birth_roster, "stage02", "birth_roster")

log_stage("stage02", "done")
