# =====================================================================
# lib/harmonize.R — Cross-wave identifier harmonization
#
# IFLS uses `commid93` as the stable community anchor across waves
# (1993 community classification). Each wave also emits its own commidXX
# which maps back via a wave-to-wave mapping table. This file provides
# a single `harmonize_commid93()` helper that resolves any wave's commid
# to the 1993 anchor.
#
# Stage 01 loads the raw `ptrack` tables which carry commid93 directly
# for individuals; stage 02 uses those to backfill commid93 onto wave-
# specific rosters. Community modules (SAR, PKK) carry their own wave
# commid and must be mapped via the community-panel anchor in stage 03.
# =====================================================================

#' Extract (commidNN -> commid93) mapping from any ptrack table.
#'
#' Ptrack in later waves (W3+) contains both the person's original 1993
#' community (`commid93`) AND the survey-wave community. For W2 ptrack,
#' `commid93` equals the wave commid. This mapping is used in stage 03
#' to anchor SAR/PKK community modules to the 1993 classification.
#'
#' @param ptrack A ptrack tibble loaded via `read_ifls_dta()`.
#' @param wave_commid Name of the wave's commid column, e.g. "commid07".
#' @return Distinct (wave_commid, commid93) tibble.
commid_mapping <- function(ptrack, wave_commid) {
  if (!all(c(wave_commid, "commid93") %in% names(ptrack))) {
    stop(sprintf(
      "ptrack missing columns: need both '%s' and 'commid93'", wave_commid))
  }
  ptrack |>
    dplyr::select(dplyr::all_of(c(wave_commid, "commid93"))) |>
    dplyr::filter(!is.na(.data[[wave_commid]]), !is.na(commid93)) |>
    dplyr::distinct()
}

#' Backfill `commid93` onto a wave-specific roster using a mapping.
#'
#' @param df Data frame with `wave_commid` but no `commid93`.
#' @param mapping Output of `commid_mapping()` for the matching wave.
#' @param wave_commid Name of the wave's commid column in `df`.
#' @return `df` with a `commid93` column appended.
attach_commid93 <- function(df, mapping, wave_commid) {
  by <- stats::setNames(wave_commid, wave_commid)
  dplyr::left_join(df, mapping, by = by)
}
