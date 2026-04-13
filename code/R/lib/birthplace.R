# =====================================================================
# lib/birthplace.R — Three `commid_birth` variants
#
# Birthplace in IFLS is not directly observed. Each variant below is one
# defensible proxy; stage 02 emits all three so the regression plan can
# pick and compare without re-running the pipeline.
#
#   commid_birth_legacy    — the child's own commid93 if the CHILD has
#                            never moved (`non_migrant == 1`). Matches
#                            Wanzi's master's thesis construction.
#
#   commid_birth_mother    — the child's commid93 if the MOTHER has never
#                            moved. Stronger birthplace claim, especially
#                            for post-1993 births where "child didn't
#                            move" is only informative conditional on the
#                            mother's stability.
#
#   commid_birth_preceding — the MOTHER's commid at the wave closest
#                            preceding the child's birth (W1 1993,
#                            W2 1997, W3 2000). For children born
#                            before 1993 or mothers not observed
#                            pre-birth, falls back to commid93.
#
# All three target `commid93` (1993 community classification) so SAR/PKK
# rollout joined on commid93 in stage 03 is directly usable.
# =====================================================================

#' Build the three commid_birth variants on a birth-roster tibble.
#'
#' Expects input columns: `pidlink`, `mother_pidlink`, `birth_year`,
#' `commid93`, `non_migrant` (0/1, child-centric), `mother_non_migrant`
#' (0/1), plus a `mother_commid_preceding` column — the mother's
#' commid93-equivalent at the wave closest preceding `birth_year`.
#'
#' @return Input tibble with three new columns:
#'   `commid_birth_legacy`, `commid_birth_mother`, `commid_birth_preceding`.
build_commid_birth_variants <- function(df) {
  stopifnot(all(
    c("pidlink", "mother_pidlink", "birth_year",
      "commid93", "non_migrant", "mother_non_migrant",
      "mother_commid_preceding") %in% names(df)
  ))

  df |>
    dplyr::mutate(
      commid_birth_legacy     = dplyr::if_else(non_migrant == 1,
                                               commid93, NA_character_),
      commid_birth_mother     = dplyr::if_else(mother_non_migrant == 1,
                                               commid93, NA_character_),
      commid_birth_preceding  = dplyr::coalesce(mother_commid_preceding,
                                                commid93)
    )
}

#' Build each person's mother's commid at the nearest wave preceding the
#' child's birth year.
#'
#' Vectorized: builds per-wave wide lookup tables and coalesces by
#' wave-anchor priority. For a child born in `birth_year`:
#'   * birth_year <= 1993 -> no pre-birth IFLS wave -> NA
#'   * 1994..1997 -> mother's W1 (1993) commid
#'   * 1998..2000 -> mother's W2 (1997), fallback W1
#'   * >= 2001    -> mother's W3 (2000), fallback W2 then W1
#'
#' @param mother_pidlink Character vector.
#' @param birth_year Integer/numeric vector.
#' @param mother_commid_panel Tibble with (pidlink, wave, commid93) where
#'        `wave` is one of "W1","W2","W3","W4","W5".
#' @return Character vector of commid93 codes, NA where unresolvable.
mother_commid_preceding <- function(mother_pidlink, birth_year,
                                     mother_commid_panel) {
  df <- tibble::tibble(
    mother_pidlink = mother_pidlink,
    birth_year     = as.numeric(birth_year)
  )

  wide <- mother_commid_panel |>
    dplyr::filter(!is.na(pidlink), !is.na(commid93),
                  wave %in% c("W1", "W2", "W3")) |>
    dplyr::distinct(pidlink, wave, .keep_all = TRUE) |>
    tidyr::pivot_wider(names_from = wave, values_from = commid93,
                       names_prefix = "c_") |>
    dplyr::rename(mother_pidlink = pidlink)

  # Ensure all three columns exist even if a wave had no rows.
  for (col in c("c_W1", "c_W2", "c_W3")) {
    if (!col %in% names(wide)) wide[[col]] <- NA_character_
  }

  joined <- dplyr::left_join(df, wide, by = "mother_pidlink")

  dplyr::case_when(
    joined$birth_year >= 2001 ~ dplyr::coalesce(joined$c_W3,
                                                joined$c_W2,
                                                joined$c_W1),
    joined$birth_year >= 1998 ~ dplyr::coalesce(joined$c_W2, joined$c_W1),
    joined$birth_year >= 1994 ~ joined$c_W1,
    TRUE                      ~ NA_character_
  )
}
