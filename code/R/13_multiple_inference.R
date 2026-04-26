# =====================================================================
# 13_multiple_inference.R — Multiple-hypothesis-testing corrections.
#
# When testing k outcomes at alpha = 0.05 each, the family-wise false-
# positive rate blows up to 1 - (1-alpha)^k (≈19% for k=4). This stage
# reports two standard corrections alongside the naive per-outcome
# p-values from stage 08:
#
#   1. Romano-Wolf (2005) stepdown p-values across the four composite
#      indices (FWER control). Pair cluster bootstrap at commid_birth
#      with B draws; center t-stats under the null via the standard
#      (beta*_b - beta_obs)/se*_b transformation; apply stepdown.
#
#   2. Anderson (2008) sharpened FDR q-values within each outcome
#      domain's component-level tests (health: 6, cognition: 4,
#      bigfive: 5). Depression is a single composite (CESD-10 sum) so
#      it has no within-domain components to sharpen across; it is
#      included in the cross-domain Romano-Wolf correction only.
#
# Inputs :  data/intermediate/stage06/analysis_sample.rds
# Outputs:  paper/tables/romano_wolf.tex
#           paper/tables/anderson_qvalues.tex
#           data/output/regressions/stage13_mult_inf.rds
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(fixest)
})

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))
source(here::here("code", "R", "lib", "labels.R"))

log_stage("stage13", "begin")

# ---------------------------------------------------------------------
# 0. Configuration.
# ---------------------------------------------------------------------
B <- 1999L   # Romano-Wolf bootstrap draws. Publication-standard.
ALPHA <- 0.05

# Significance stars from a p- or q-value vector (standard econ convention).
add_stars <- function(p) {
  ifelse(is.na(p),       "",
  ifelse(p < 0.01,  "$^{***}$",
  ifelse(p < 0.05,   "$^{**}$",
  ifelse(p < 0.10,    "$^{*}$", ""))))
}

frame <- read_intermediate("stage06", "analysis_sample")

tables_dir <- here::here("paper", "tables")
reg_dir    <- here::here("data", "output", "regressions")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reg_dir,    recursive = TRUE, showWarnings = FALSE)

df <- frame |>
  filter(primary_sample == 1) |>
  mutate(
    sex_f = factor(sex, levels = c(1, 3), labels = c("male", "female"))
  )

outcomes <- c("health_index", "cognition_index",
              "depression_index", "bigfive_index")

# ---------------------------------------------------------------------
# 1. Observed coefficients & t-stats from the main TWFE + controls spec
#    (column (2) of stage 08's main_results.tex).
# ---------------------------------------------------------------------
fml <- function(y) {
  as.formula(sprintf(paste0(
    "%s ~ exposure_early_sar_legacy + sex_f + mother_edu_years + ",
    "mother_age_birth | commid_birth_legacy + birth_year + source_wave"), y))
}

obs_stats <- lapply(outcomes, function(y) {
  m <- feols(fml(y), data = df, cluster = ~ commid_birth_legacy)
  co <- coef(m)["exposure_early_sar_legacy"]
  se <- sqrt(vcov(m)["exposure_early_sar_legacy",
                     "exposure_early_sar_legacy"])
  tibble(outcome = y,
         beta   = unname(co),
         se     = unname(se),
         t      = unname(co / se),
         p_raw  = 2 * pnorm(-abs(unname(co / se))))
}) |> bind_rows()

log_stage("stage13", "observed t-stats computed; starting Romano-Wolf bootstrap")

# ---------------------------------------------------------------------
# 2. Romano-Wolf stepdown bootstrap.
#    Pair cluster bootstrap at commid_birth: for each draw, sample
#    clusters with replacement, re-estimate all 4 outcomes, compute
#    centered t*_k,b = (beta*_k,b - beta_obs_k) / se*_k,b.
# ---------------------------------------------------------------------
cluster_col <- "commid_birth_legacy"
cluster_ids <- unique(df[[cluster_col]])
n_clust     <- length(cluster_ids)
cluster_rows <- split(seq_len(nrow(df)), df[[cluster_col]])

t_boot <- matrix(NA_real_, nrow = B, ncol = length(outcomes),
                 dimnames = list(NULL, outcomes))

t0 <- Sys.time()
for (b in seq_len(B)) {
  # Resample clusters with replacement; build row index & fresh
  # cluster labels so feols's cluster-SE sees each draw as unique.
  sampled <- sample(cluster_ids, n_clust, replace = TRUE)
  row_lists <- cluster_rows[as.character(sampled)]
  row_idx   <- unlist(row_lists, use.names = FALSE)
  new_cl    <- rep(seq_along(sampled),
                   vapply(row_lists, length, integer(1)))

  boot_df <- df[row_idx, , drop = FALSE]
  boot_df$.boot_cluster <- new_cl

  for (k in seq_along(outcomes)) {
    y <- outcomes[k]
    fit <- tryCatch(
      feols(as.formula(sprintf(paste0(
        "%s ~ exposure_early_sar_legacy + sex_f + mother_edu_years + ",
        "mother_age_birth | .boot_cluster + birth_year"), y)),
        data = boot_df, cluster = ~ .boot_cluster,
        notes = FALSE, warn = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    co_b <- coef(fit)["exposure_early_sar_legacy"]
    se_b <- sqrt(vcov(fit)["exposure_early_sar_legacy",
                           "exposure_early_sar_legacy"])
    if (!is.finite(se_b) || se_b == 0) next
    # Centered t-stat (null-imposing): shift by the observed point estimate.
    t_boot[b, k] <- (co_b - obs_stats$beta[k]) / se_b
  }

  if (b %% 200 == 0) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    log_stage("stage13",
              sprintf("RW bootstrap: %d / %d draws (%.1fs)", b, B, elapsed))
  }
}

# ---------------------------------------------------------------------
# 3. Romano-Wolf stepdown p-values.
#    Order hypotheses by observed |t| from largest to smallest.
#    p_{(i)} = mean_b[max over remaining k >= i of |t*_k,b| >= |t_{(i)}|]
#    then enforce monotonicity: p_{(i)} = max(p_{(i-1)}, p_{(i)})
# ---------------------------------------------------------------------
t_abs_obs <- abs(obs_stats$t)
ord <- order(t_abs_obs, decreasing = TRUE)     # order of hypotheses

# Drop bootstrap rows with any missing t (should be rare; logged if any).
n_miss_any <- sum(!complete.cases(t_boot))
if (n_miss_any > 0) {
  log_stage("stage13", sprintf(
    "RW bootstrap: %d / %d draws had a failed re-fit; dropping those rows.",
    n_miss_any, B))
  t_boot <- t_boot[complete.cases(t_boot), , drop = FALSE]
}
B_eff <- nrow(t_boot)

p_rw_ordered <- numeric(length(outcomes))
remaining <- ord
for (i in seq_along(ord)) {
  t_max_remaining <- apply(
    abs(t_boot[, remaining, drop = FALSE]), 1, max
  )
  p_rw_ordered[i] <- mean(t_max_remaining >= t_abs_obs[ord[i]])
  remaining <- remaining[-1]
}
# Enforce monotonicity in the ordered sequence.
for (i in seq_along(p_rw_ordered)[-1]) {
  p_rw_ordered[i] <- max(p_rw_ordered[i - 1], p_rw_ordered[i])
}

# Map back to outcome order.
p_rw <- numeric(length(outcomes))
p_rw[ord] <- p_rw_ordered

obs_stats$p_rw <- p_rw
obs_stats$stars_raw <- add_stars(obs_stats$p_raw)
obs_stats$stars_rw  <- add_stars(obs_stats$p_rw)

log_stage("stage13", sprintf(
  "Romano-Wolf complete. B_eff=%d. p_RW: %s",
  B_eff,
  paste(sprintf("%s=%.3f", pretty_label(obs_stats$outcome), obs_stats$p_rw),
        collapse = "; ")))

# ---------------------------------------------------------------------
# 4. Romano-Wolf table output.
# ---------------------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("---")
  if (p < 0.001) "$<0.001$" else sprintf("%.3f", p)
}

writeLines(c(
  "% Auto-generated by code/R/13_multiple_inference.R",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  paste("Outcome & $\\hat{\\beta}$ & SE &",
        "Naive $p$ & Romano--Wolf $p$ \\\\"),
  "\\midrule",
  apply(obs_stats, 1, function(r) sprintf(
    "%s & %.3f%s & (%.3f) & %s & %s%s \\\\",
    pretty_label(r["outcome"]),
    as.numeric(r["beta"]),
    as.character(r["stars_raw"]),
    as.numeric(r["se"]),
    fmt_p(as.numeric(r["p_raw"])),
    fmt_p(as.numeric(r["p_rw"])),
    as.character(r["stars_rw"])
  )),
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "romano_wolf.tex"))

log_stage("stage13", "Romano-Wolf table written.")

# ---------------------------------------------------------------------
# 5. Anderson sharpened FDR q-values on component-level tests.
#    For each domain, run TWFE + controls on every component and
#    apply Benjamini-Hochberg (= Anderson 2008 sharpening).
#    Depression is a single composite; no within-domain components.
# ---------------------------------------------------------------------
domain_components <- list(
  "Health"         = c("height_cm", "no_underweight", "lung_log",
                       "no_hypertension", "good_early_health", "ghs"),
  "Cognition"      = c("w_abil", "serial7_correct",
                       "immed_count", "delay_count"),
  "Socioemotional" = c("agreeableness_z", "conscientiousness_z",
                       "extraversion_z", "emotional_stability_z",
                       "openness_z")
)

# The `ghs` component is sign-flipped (high = bad health) in the
# Anderson index; for component-level FDR we want each test oriented
# so that a POSITIVE beta means a program BENEFIT. Flip ghs before
# regressing.
df$ghs_flip <- -df$ghs
component_flip <- list("Health" = c(ghs = "ghs_flip"))

component_rows <- lapply(names(domain_components), function(domain) {
  comps <- domain_components[[domain]]
  flips <- component_flip[[domain]]
  if (!is.null(flips)) comps[comps %in% names(flips)] <- flips

  rows <- lapply(comps, function(v) {
    fit <- tryCatch(
      feols(as.formula(sprintf(paste0(
        "%s ~ exposure_early_sar_legacy + sex_f + mother_edu_years + ",
        "mother_age_birth | commid_birth_legacy + birth_year + source_wave"), v)),
        data = df, cluster = ~ commid_birth_legacy,
        notes = FALSE, warn = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(tibble(domain = domain, component = v,
                    beta = NA_real_, se = NA_real_,
                    t = NA_real_, p_raw = NA_real_))
    }
    co <- coef(fit)["exposure_early_sar_legacy"]
    se <- sqrt(vcov(fit)["exposure_early_sar_legacy",
                         "exposure_early_sar_legacy"])
    tt <- unname(co / se)
    tibble(domain = domain, component = v,
           beta = unname(co), se = unname(se),
           t = tt, p_raw = 2 * pnorm(-abs(tt)))
  }) |> bind_rows()

  # Anderson sharpening = Benjamini-Hochberg (monotonicity in sorted
  # p-values after scaling by K/i). p.adjust(method = "BH") implements
  # exactly this.
  rows$q_anderson <- stats::p.adjust(rows$p_raw, method = "BH")
  rows
}) |> bind_rows()

# Pretty names for the component column.
component_pretty <- c(
  height_cm             = "Height (cm)",
  no_underweight        = "Not underweight (BMI)",
  lung_log              = "Log peak flow",
  no_hypertension       = "Not hypertensive (BP)",
  good_early_health     = "Good childhood health (recall)",
  ghs_flip              = "General health (self-rated, flipped)",
  w_abil                = "Word ability (IRT)",
  serial7_correct       = "Serial-7 subtractions",
  immed_count           = "Immediate word recall (10)",
  delay_count           = "Delayed word recall (10)",
  agreeableness_z       = "Agreeableness",
  conscientiousness_z   = "Conscientiousness",
  extraversion_z        = "Extraversion",
  emotional_stability_z = "Emotional stability",
  openness_z            = "Openness"
)
component_rows$component_label <-
  unname(component_pretty[component_rows$component])
component_rows$component_label[is.na(component_rows$component_label)] <-
  gsub("_", " ", component_rows$component[is.na(component_rows$component_label)])

component_rows$stars_raw <- add_stars(component_rows$p_raw)
component_rows$stars_q   <- add_stars(component_rows$q_anderson)

writeLines(c(
  "% Auto-generated by code/R/13_multiple_inference.R",
  "\\begin{tabular}{llcccc}",
  "\\toprule",
  paste("Domain & Component & $\\hat{\\beta}$ & SE &",
        "Naive $p$ & Anderson $q$ \\\\"),
  "\\midrule",
  unlist(lapply(unique(component_rows$domain), function(d) {
    sub <- component_rows[component_rows$domain == d, , drop = FALSE]
    lines <- apply(sub, 1, function(r) sprintf(
      "%s & %s & %.3f%s & (%.3f) & %s & %s%s \\\\",
      ifelse(r["component"] == sub$component[1], d, ""),
      as.character(r["component_label"]),
      as.numeric(r["beta"]),
      as.character(r["stars_raw"]),
      as.numeric(r["se"]),
      fmt_p(as.numeric(r["p_raw"])),
      fmt_p(as.numeric(r["q_anderson"])),
      as.character(r["stars_q"])
    ))
    c(lines, "\\midrule")
  })) |> head(-1),   # drop trailing \midrule
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "anderson_qvalues.tex"))

log_stage("stage13", "Anderson q-value table written.")

# ---------------------------------------------------------------------
# 6. Persist full stats object for downstream use + reproducibility.
# ---------------------------------------------------------------------
saveRDS(list(
  romano_wolf     = obs_stats,
  anderson_fdr    = component_rows,
  bootstrap_draws = B_eff,
  alpha           = ALPHA
), file.path(reg_dir, "stage13_mult_inf.rds"))

log_stage("stage13", "done")
