# =====================================================================
# lib/indices.R — Composite-index helpers (Anderson 2008, ICW-weighted)
#
# Used by stage 06 to collapse outcome families (health, cognition,
# depression, Big-5) into a single index per family. The index is
# inverse-covariance-weighted per Anderson (2008, JASA) so components
# with redundant information are downweighted.
#
# Convention: every input component is sign-oriented so that HIGHER
# values mean BETTER outcome. Callers are responsible for flipping
# signs before passing (e.g. negate `dep_index` so higher = less depressed).
# =====================================================================

#' z-score a vector using complete-case mean and sd.
#'
#' NA entries are preserved and do not enter the moments. Constant
#' columns return zeros (with a warning).
#' @param x numeric vector
zscore <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sd <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sd) || sd == 0) {
    warning("zscore: zero/NA sd; returning zeros.")
    return(ifelse(is.na(x), NA_real_, 0))
  }
  (x - mu) / sd
}

#' Anderson (2008) inverse-covariance-weighted index.
#'
#' Steps:
#'   1. z-score each component on rows with ≥1 non-NA (per-column moments).
#'   2. Replace remaining NAs with 0 (neutral on the standardized scale).
#'   3. Compute weights w = rowSums(inv(S)) where S is the covariance of
#'      the completed matrix. This downweights redundant components.
#'   4. Index = z %*% w, then re-standardized (mean 0, sd 1 on non-missing).
#'
#' Rows with ALL components NA return NA.
#'
#' @param df A data frame or tibble.
#' @param vars Character vector of column names in `df`. Each column must
#'   be sign-oriented so higher = better.
#' @return A numeric vector of length `nrow(df)`.
anderson_index <- function(df, vars) {
  stopifnot(all(vars %in% names(df)))
  X <- as.matrix(df[, vars, drop = FALSE])
  storage.mode(X) <- "double"

  # Per-column z-scoring (NA-aware).
  Z <- apply(X, 2, zscore)
  if (!is.matrix(Z)) Z <- matrix(Z, ncol = length(vars))
  colnames(Z) <- vars

  all_na <- rowSums(!is.na(Z)) == 0

  # Imputation for weight computation only: NA -> 0 on z-scale (group mean).
  Z_filled <- Z
  Z_filled[is.na(Z_filled)] <- 0

  # Covariance + inverse. Guard against singularity with a tiny ridge.
  S <- stats::cov(Z_filled)
  ridge <- 1e-6 * mean(diag(S))
  S_inv <- tryCatch(
    solve(S + diag(ridge, nrow(S))),
    error = function(e) {
      warning("anderson_index: cov matrix near-singular; using pseudo-inverse.")
      MASS::ginv(S + diag(ridge, nrow(S)))
    }
  )

  # Per-component weight = row-sum of inverse covariance.
  w <- rowSums(S_inv)
  names(w) <- vars

  # Weighted sum, re-scored.
  idx <- as.numeric(Z_filled %*% w)
  idx[all_na] <- NA_real_
  as.numeric(zscore(idx))
}

#' Construct a signed outcome frame from raw columns.
#'
#' For each spec element `list(var = "<name>", flip = FALSE)`, pulls the
#' column from `df` and optionally flips its sign. Returns a data frame
#' with the same row count and one column per spec element named after
#' `var`. Useful for building the input matrix to `anderson_index()`.
#'
#' @param df Source data frame.
#' @param spec Named list of `list(var, flip)` elements.
orient_components <- function(df, spec) {
  out <- lapply(spec, function(s) {
    x <- df[[s$var]]
    if (isTRUE(s$flip)) -x else x
  })
  names(out) <- vapply(spec, function(s) s$var, character(1))
  as.data.frame(out)
}
