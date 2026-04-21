# =====================================================================
# 15_honestdid_sensitivity.R — Rambachan-Roth PTA sensitivity for the
# health-index event study.
#
# The CS event study in stage 09 shows a significant negative pre-
# period coefficient at event-time -3 on the health index (coef -0.326,
# SE 0.120, t -2.72). This is a parallel-trends-assumption (PTA)
# violation at a specific lead, and on its own it prevents a clean
# reading of the post-treatment CS coefficients.
#
# The Rambachan-Roth (2023) "Honest Approach to Parallel Trends"
# framework provides a robust alternative: it asks how large a PTA
# violation the post-treatment estimate can tolerate while remaining
# informative. We impose the "relative magnitudes" restriction: post-
# treatment deviations from parallel trends are at most Mbar times the
# largest pre-treatment deviation. Mbar = 0 recovers exact parallel
# trends; Mbar = 1 says post-deviations can be as large as the worst
# observed pre-deviation; larger Mbar relaxes further. The "breakdown"
# Mbar is the largest value at which the robust confidence set still
# excludes zero.
#
# Inputs :  data/intermediate/stage06/analysis_sample.rds
# Outputs:  paper/figures/honestdid_health_sensitivity.pdf
#           paper/tables/honestdid_health_breakdown.tex
#           paper/tables/honestdid_trimmed_pre.tex
#
# Two specifications are reported: (i) the baseline pre-period window
# (e in {-3, -2, -1}) that matches the CS event study, and (ii) a
# trimmed pre-period window (e in {-2, -1}) that drops the extreme
# lead where cohort-village cells become thin. Dropping e = -3
# observations does not change the post-treatment CS point estimate
# because CS group-time ATTs are independent across (G, t) pairs; it
# only changes the pre-period evidence HonestDiD uses to bound the
# admissible PTA violation.
# =====================================================================

set.seed(20260412)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(fixest)
  library(ggplot2)
  library(HonestDiD)
})

source(here::here("code", "R", "lib", "paths.R"))
source(here::here("code", "R", "lib", "io.R"))
source(here::here("code", "R", "lib", "labels.R"))

log_stage("stage15", "begin")

frame <- read_intermediate("stage06", "analysis_sample")
tables_dir  <- here::here("paper", "tables")
figures_dir <- here::here("paper", "figures")

palette_oi <- c(blue = "#0072B2", vermilion = "#D55E00",
                black = "#000000")

# ---------------------------------------------------------------------
# 1. Pull the CS event-time coefficients and vcov.
#    Refit att_gt on the non-migrant sample with bootstrap SEs so the
#    influence-function-based vcov for the event-time aggregation is
#    available. The CS event study is the one that shows the t=-3
#    pre-period violation the sensitivity is designed to address.
# ---------------------------------------------------------------------
suppressPackageStartupMessages(library(did))

df <- frame |>
  filter(primary_sample == 1) |>
  mutate(
    G = if_else(!is.na(start_sar_legacy),
                as.integer(start_sar_legacy) - 1L, 0L),
    commid_num = as.integer(factor(commid_birth_legacy))
  ) |>
  filter(!is.na(health_index), !is.na(birth_year),
         !is.na(commid_num)) |>
  as.data.frame()

att <- att_gt(
  yname = "health_index", tname = "birth_year",
  idname = "commid_num", gname = "G", data = df,
  control_group = "nevertreated",
  bstrap = TRUE, biters = 1000, cband = FALSE,
  allow_unbalanced_panel = TRUE, panel = FALSE,
  clustervars = NULL
)

ET_MIN <- -3L
ET_MAX <-  4L

dyn <- aggte(att, type = "dynamic",
             min_e = ET_MIN, max_e = ET_MAX,
             na.rm = TRUE)

# Influence-function-based vcov for the event-time aggregated ATTs.
# `dyn$inf.function$dynamic.inf.func.e` is N x n_event-times.
inf_e <- dyn$inf.function$dynamic.inf.func.e
n_obs <- nrow(inf_e)
sigma  <- crossprod(inf_e) / n_obs / n_obs

et_vec  <- dyn$egt
betahat <- dyn$att.egt

# Order by event time and drop -1 (reference), since HonestDiD treats
# -1 as the omitted reference period implicitly through numPrePeriods.
keep <- et_vec != -1L
et_vec  <- et_vec[keep]
betahat <- betahat[keep]
sigma   <- sigma[keep, keep, drop = FALSE]
ord <- order(et_vec)
et_vec  <- et_vec[ord]
betahat <- betahat[ord]
sigma   <- sigma[ord, ord, drop = FALSE]

log_stage("stage15", sprintf(
  "CS event-time estimates (non-migrant sample): %s",
  paste(sprintf("t=%d: %.3f (SE %.3f)", et_vec, betahat,
                sqrt(diag(sigma))),
        collapse = "; ")))

# Pre-period indices (event_time <= -2) and post-period (>= 0).
num_pre  <- sum(et_vec <= -2L)
num_post <- sum(et_vec >=  0L)
stopifnot(num_pre > 0, num_post > 0)

# ---------------------------------------------------------------------
# 2. Sensitivity analysis on the average of the first four post-
#    treatment event times (event-times 0..3), which matches the
#    headline CS aggregation window.
# ---------------------------------------------------------------------
n_post_avg <- min(4L, num_post)
l_vec <- rep(0, num_post)
l_vec[seq_len(n_post_avg)] <- 1 / n_post_avg

orig <- constructOriginalCS(
  betahat = betahat, sigma = sigma,
  numPrePeriods  = num_pre,
  numPostPeriods = num_post,
  l_vec = l_vec, alpha = 0.05
)

Mbarvec <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.50, 1, 2)
sens <- createSensitivityResults_relativeMagnitudes(
  betahat = betahat, sigma = sigma,
  numPrePeriods  = num_pre,
  numPostPeriods = num_post,
  l_vec = l_vec, alpha = 0.05,
  Mbarvec = Mbarvec
)

# Breakdown: largest Mbar at which the robust CI still excludes zero.
sens <- sens |>
  mutate(excl_zero = !(lb <= 0 & ub >= 0))
breakdown <- if (any(sens$excl_zero)) {
  max(sens$Mbar[sens$excl_zero])
} else {
  0
}

log_stage("stage15", sprintf(
  "HonestDiD: original CS = [%.3f, %.3f]; breakdown Mbar = %.2f",
  orig$lb, orig$ub, breakdown))

# ---------------------------------------------------------------------
# 3. Plot the sensitivity curve. Style mirrors the event-study plot:
#    theme_bw, blue error bars on the robust bounds, vermilion
#    reference line for the original CS.
# ---------------------------------------------------------------------
sens_plot_df <- sens |>
  transmute(Mbar = as.numeric(Mbar),
            lb   = as.numeric(lb),
            ub   = as.numeric(ub),
            kind = "Robust CI")

orig_df <- tibble(Mbar = -0.1, lb = orig$lb, ub = orig$ub,
                  kind = "Original CS CI")

p_sens <- ggplot(mapping = aes(x = Mbar, ymin = lb, ymax = ub)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = unname(palette_oi["black"]), alpha = 0.5) +
  geom_errorbar(data = orig_df, width = 0.05, linewidth = 0.6,
                colour = unname(palette_oi["vermilion"])) +
  geom_errorbar(data = sens_plot_df, width = 0.05, linewidth = 0.6,
                colour = unname(palette_oi["blue"])) +
  scale_x_continuous(
    breaks = c(-0.1, Mbarvec),
    labels = c("Original", sprintf("%.2f", Mbarvec))
  ) +
  labs(x = expression("Original CS CI vs. robust CI at"~bar(M)),
       y = "Avg. ATT on health index over event-times 0--3 (95% CI)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9))

ggsave(file.path(figures_dir, "honestdid_health_sensitivity.pdf"),
       plot = p_sens, width = 6.5, height = 4.0, units = "in",
       device = cairo_pdf, bg = "white")

# ---------------------------------------------------------------------
# 4. Breakdown table.
# ---------------------------------------------------------------------
orig_excl <- if (orig$lb > 0 || orig$ub < 0) "Yes" else "No"

sens_rows <- vapply(seq_len(nrow(sens)), function(i) {
  mark <- if (sens$excl_zero[i]) "Yes" else "No"
  sprintf("%.2f & %.3f & %.3f & %s \\\\",
          sens$Mbar[i], sens$lb[i], sens$ub[i], mark)
}, character(1))

writeLines(c(
  "% Auto-generated by code/R/15_honestdid_sensitivity.R",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Specification & Lower bound & Upper bound & Excludes zero? \\\\",
  "\\midrule",
  sprintf("Original CS CI & %.3f & %.3f & %s \\\\",
          orig$lb, orig$ub, orig_excl),
  "\\midrule",
  "\\multicolumn{4}{l}{\\emph{Robust CI at relative-magnitudes parameter $\\bar{M}$}} \\\\",
  sens_rows,
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "honestdid_health_breakdown.tex"))

# ---------------------------------------------------------------------
# 5. Trimmed pre-period (only e = -2 retained as the pre-period
#    baseline). This addresses the concern that the t = -3 pre-period
#    spike is driven by thin cohort-village cells at the extreme lead.
#    The CS post-treatment ATTs are identical by construction; only
#    the HonestDiD bound tightens because the pre-period dispersion
#    (which the relative-magnitudes restriction caps against) is
#    smaller on this window.
# ---------------------------------------------------------------------
dyn_trim <- aggte(att, type = "dynamic",
                  min_e = -2L, max_e = ET_MAX,
                  na.rm = TRUE)
inf_e_t  <- dyn_trim$inf.function$dynamic.inf.func.e
sigma_t  <- crossprod(inf_e_t) / n_obs / n_obs
et_t     <- dyn_trim$egt
betahat_t <- dyn_trim$att.egt
keep_t   <- et_t != -1L
et_t     <- et_t[keep_t]
betahat_t <- betahat_t[keep_t]
sigma_t  <- sigma_t[keep_t, keep_t, drop = FALSE]
ord_t    <- order(et_t)
et_t     <- et_t[ord_t]
betahat_t <- betahat_t[ord_t]
sigma_t  <- sigma_t[ord_t, ord_t, drop = FALSE]

num_pre_t  <- sum(et_t <= -2L)
num_post_t <- sum(et_t >=  0L)
stopifnot(num_pre_t >= 1, num_post_t > 0)

n_post_avg_t <- min(4L, num_post_t)
l_vec_t <- rep(0, num_post_t)
l_vec_t[seq_len(n_post_avg_t)] <- 1 / n_post_avg_t

orig_t <- constructOriginalCS(
  betahat = betahat_t, sigma = sigma_t,
  numPrePeriods  = num_pre_t,
  numPostPeriods = num_post_t,
  l_vec = l_vec_t, alpha = 0.05
)
sens_t <- createSensitivityResults_relativeMagnitudes(
  betahat = betahat_t, sigma = sigma_t,
  numPrePeriods  = num_pre_t,
  numPostPeriods = num_post_t,
  l_vec = l_vec_t, alpha = 0.05,
  Mbarvec = Mbarvec
)
sens_t <- sens_t |>
  mutate(excl_zero = !(lb <= 0 & ub >= 0))
breakdown_t <- if (any(sens_t$excl_zero)) {
  max(sens_t$Mbar[sens_t$excl_zero])
} else {
  0
}

log_stage("stage15", sprintf(
  "HonestDiD (trimmed pre-period, e=-2): original CS = [%.3f, %.3f]; breakdown Mbar = %.2f",
  orig_t$lb, orig_t$ub, breakdown_t))

# Table: baseline vs. trimmed side-by-side.
orig_excl_t <- if (orig_t$lb > 0 || orig_t$ub < 0) "Yes" else "No"
rows_trim <- vapply(seq_len(nrow(sens_t)), function(i) {
  mark <- if (sens_t$excl_zero[i]) "Yes" else "No"
  sprintf("%.2f & %.3f & %.3f & %s \\\\",
          sens_t$Mbar[i], sens_t$lb[i], sens_t$ub[i], mark)
}, character(1))
writeLines(c(
  "% Auto-generated by code/R/15_honestdid_sensitivity.R",
  "\\begin{tabular}{lccc}",
  "\\toprule",
  "Specification & Lower bound & Upper bound & Excludes zero? \\\\",
  "\\midrule",
  sprintf("Original CS CI & %.3f & %.3f & %s \\\\",
          orig_t$lb, orig_t$ub, orig_excl_t),
  "\\midrule",
  "\\multicolumn{4}{l}{\\emph{Robust CI at relative-magnitudes parameter $\\bar{M}$}} \\\\",
  rows_trim,
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "honestdid_trimmed_pre.tex"))

# ---------------------------------------------------------------------
# Combined side-by-side table: baseline vs. trimmed at each Mbar.
# ---------------------------------------------------------------------
stopifnot(all(sens$Mbar == sens_t$Mbar))
combined_rows <- vapply(seq_len(nrow(sens)), function(i) {
  mark_b <- if (sens$excl_zero[i])   "Yes" else "No"
  mark_t <- if (sens_t$excl_zero[i]) "Yes" else "No"
  sprintf("%.2f & [%.3f, %.3f] & %s & [%.3f, %.3f] & %s \\\\",
          sens$Mbar[i],
          sens$lb[i],   sens$ub[i],   mark_b,
          sens_t$lb[i], sens_t$ub[i], mark_t)
}, character(1))

writeLines(c(
  "% Auto-generated by code/R/15_honestdid_sensitivity.R",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  paste("$\\bar{M}$ & Baseline (e$\\in\\{-3,-2\\}$) CI &",
        "Excludes 0? & Trimmed (e$=-2$) CI & Excludes 0? \\\\"),
  "\\midrule",
  sprintf("Original & [%.3f, %.3f] & %s & [%.3f, %.3f] & %s \\\\",
          orig$lb, orig$ub, orig_excl,
          orig_t$lb, orig_t$ub, orig_excl_t),
  "\\midrule",
  combined_rows,
  "\\bottomrule",
  "\\end{tabular}"
), file.path(tables_dir, "honestdid_combined.tex"))

# ---------------------------------------------------------------------
# 6. Joint sensitivity plot: baseline bounds + trimmed bounds on the
#    same axes, labelled by specification.
# ---------------------------------------------------------------------
joint_df <- bind_rows(
  sens |>
    transmute(Mbar = as.numeric(Mbar),
              lb   = as.numeric(lb),
              ub   = as.numeric(ub),
              spec = sprintf("Baseline pre-period (e in -3,-2)")),
  sens_t |>
    transmute(Mbar = as.numeric(Mbar),
              lb   = as.numeric(lb),
              ub   = as.numeric(ub),
              spec = sprintf("Trimmed pre-period (e = -2 only)"))
) |>
  mutate(spec = factor(spec,
                       levels = c("Baseline pre-period (e in -3,-2)",
                                  "Trimmed pre-period (e = -2 only)")))

p_joint <- ggplot(joint_df, aes(x = Mbar, ymin = lb, ymax = ub,
                                 colour = spec)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = unname(palette_oi["black"]), alpha = 0.5) +
  geom_errorbar(position = position_dodge(width = 0.08),
                width = 0.04, linewidth = 0.6) +
  scale_colour_manual(
    values = setNames(
      c(unname(palette_oi["vermilion"]), unname(palette_oi["blue"])),
      c("Baseline pre-period (e in -3,-2)",
        "Trimmed pre-period (e = -2 only)")
    )
  ) +
  scale_x_continuous(breaks = Mbarvec,
                     labels = sprintf("%.2f", Mbarvec)) +
  labs(x = expression("Relative-magnitudes parameter"~bar(M)),
       y = "Avg. ATT on health index over event-times 0--3 (95% CI)",
       colour = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(figures_dir, "honestdid_joint_sensitivity.pdf"),
       plot = p_joint, width = 6.5, height = 4.2, units = "in",
       device = cairo_pdf, bg = "white")

log_stage("stage15", "HonestDiD sensitivity figures and tables written.")
log_stage("stage15", "done")
