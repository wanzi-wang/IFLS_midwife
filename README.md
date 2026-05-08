# Early-Life Exposure to Indonesia's Village Midwife Program and Child Human Capital

**Author:** Wanzi Wang
**Data:** IFLS waves 1–5 (restricted-use; see `data/README.md`)

---

## Abstract

Can early-life maternal care generate lasting gains in adult human capital? In 1989 Indonesia assigned a trained midwife to each rural village to cut maternal mortality. I link the 1984–1999 birth cohorts of the Indonesia Family Life Survey to the program's staggered community-level rollout and estimate long-run effects with the Callaway–Sant'Anna (CS) estimator, using never-treated communities as controls. Exposure before age three is associated with increases in an adult physical-health index of 0.248 standard deviations (SE 0.145) and a cognition index of 0.272 standard deviations (SE 0.142), both marginally significant at the 10% level and concentrated in the in-utero and age 0–3 critical-period windows; estimates on mental health and socioemotional traits are small and statistically indistinguishable from zero. A mechanism decomposition on mother-reported pregnancy histories points to a single plausible channel: village-midwife-attended births are higher by roughly 4.5 percentage points on a baseline of 8.3 percent, while aggregate skilled attendance, prenatal care, breastfeeding, and birth weight do not move. The pattern is consistent with the program reallocating deliveries toward village midwives without raising the overall rate of skilled care.

---

## How to run

```bash
# 1. Point to raw IFLS data (W1–W5 Stata .dta files, restricted-use)
export IFLS_DATA_ROOT="/Users/wwz/Academic/Thesis_Data_IFLS/old/IFLS"

# 2. Run the data-construction pipeline (≈21s cold; idempotent when cached)
Rscript code/R/00_run_all.R

# 3. Build paper and slides (once regression stages + LaTeX are written)
cd paper  && latexmk -xelatex -outdir=../output/paper paper.tex
cd slides && latexmk -xelatex -outdir=../output/slides slides.tex
```

## Pipeline

| Stage | Script | Purpose |
|-------|--------|---------|
| 01 | `01_load_ifls.R` | Raw `.dta` → typed `.rds` cache under `data/intermediate/raw/` |
| 02 | `02_birth_roster.R` | Individual panel (cohorts 1984–1999), three `commid_birth` variants, mother covariates |
| 03 | `03_midwife_rollout.R` | SAR and PKK `start_year` per commid93 + rollout comparison diagnostic |
| 04 | `04_outcomes.R` | W5 outcomes (health, adult cognition, CESD-10, Big-5) |
| 04b | `04b_mediators.R` | Five mother-reported pregnancy-history mediators |
| 05 | `05_analysis_frame.R` | Merge + exposure-timing variables × 3 birthplace variants + audit vs legacy `df.rds` |
| 06 | `06_sample_weights_indices.R` | Anderson inverse-covariance composite indices and IPW |
| 07 | `07_balance_descriptives.R` | Balance table + descriptive statistics |
| 08 | `08_main_regressions.R` | TWFE + CS + de Chaisemartin headline estimates |
| 09 | `09_event_study.R` | CS dynamic event-study coefficients |
| 10 | `10_sibling_fe.R` | Within-mother (sibling-FE) specifications |
| 11 | `11_robustness.R` | Sensitivity (IPW, PKK, alternative birthplaces) + placebos |
| 12 | `12_heterogeneity_finalize.R` | Subgroup splits by sex × education × maternal age × region |
| 13 | `13_multiple_inference.R` | Romano–Wolf stepdown across the four indices |
| 14 | `14_mechanisms.R` | First-stage on five mediators (Ahsan-style) |
| 16 | `16_appendix_figures.R` | Province-level rollout choropleth |

Library helpers live in `code/R/lib/`. Produces `data/intermediate/stage05/analysis_frame.rds` and `stage06/analysis_sample.rds` (both gitignored, regenerable). Cold run ≈70 s.

## Repository layout

Key dirs:

- `code/R/` — pipeline scripts (01..05) + `lib/` helpers
- `data/intermediate/` — gitignored; regenerable
- `paper/` — manuscript LaTeX sources, sections, tables, figures
- `slides/` — presentation sources
- `old/` — frozen pre-rebuild project, read-only reference only

## Data access

IFLS microdata is restricted-use and is **not** committed to this repo. It is referenced by path from `$IFLS_DATA_ROOT`. Raw `.dta` files are never copied into the repo; only regenerable `.rds` intermediates under `data/intermediate/`, all gitignored. See `data/README.md`.

## License

TBD.
