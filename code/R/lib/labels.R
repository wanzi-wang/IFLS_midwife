# =====================================================================
# lib/labels.R — Central map from R identifiers to publication labels.
#
# Every table / figure that reaches the paper routes variable names
# through `pretty_label()` so R snake_case never appears in the PDF.
# If a new variable shows up in a published table, add it to
# `.label_map` below — the fallback is a best-effort cleanup of
# underscores, but fallback output is not publication-grade.
# =====================================================================

.label_map <- c(
  # --- Composite indices ----------------------------------------------
  health_index      = "Health index",
  cognition_index   = "Cognition index",
  depression_index  = "Depression index",
  bigfive_index     = "Socioemotional index",

  # --- Raw outcome components -----------------------------------------
  dep_index         = "CESD-10 depression score",
  w_abil            = "Word-recall score",
  height_cm         = "Height (cm)",
  raven_score       = "Raven score",
  math_score        = "Math score",
  cog_prop_child    = "Child Raven score",

  # --- Village covariates (balance table) -----------------------------
  has_doctor              = "Has doctor",
  has_clinic              = "Has clinic",
  has_private_midwife     = "Has private midwife",
  has_traditional_midwife = "Has traditional birth attendant",

  # --- Individual / household covariates ------------------------------
  mother_edu_years = "Mother's education (years)",
  mother_age_birth = "Mother's age at birth (years)",
  birth_year       = "Birth year",
  age2014          = "Age in 2014 (years)",
  Age2014          = "Age in 2014 (years)",
  Female           = "Female",
  female           = "Female",

  # --- Design / treatment ---------------------------------------------
  exposure_early_sar_preceding = "Early exposure (age \\(\\leq 3\\))",
  commid_birth_preceding       = "Community of birth",

  # --- Heterogeneity subgroups ----------------------------------------
  `Full sample` = "Full sample",
  Male          = "Male",
  Female_       = "Female",
  `Low-edu`     = "Mother: low education",
  `High-edu`    = "Mother: high education"
)

#' Map an R identifier to a publication-ready label.
#'
#' Lookup is exact in `.label_map`. If the key is missing, falls back
#' to replacing underscores with spaces (not publication-grade —
#' always add the key to `.label_map` when this appears in a final
#' table).
#'
#' @param x character vector
#' @return character vector of same length
#' @examples
#' pretty_label("health_index")           # "Health index"
#' pretty_label(c("cognition_index", "x")) # "Cognition index", "x"
pretty_label <- function(x) {
  x <- as.character(x)
  out <- unname(.label_map[x])
  miss <- is.na(out)
  if (any(miss)) out[miss] <- gsub("_", " ", x[miss])
  out
}
