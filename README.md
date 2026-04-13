# Early-Life Exposure to Indonesia's Village Midwife Program and Child Human Capital

**Author:** Wanzi Wang
**Status:** In rebuild — journal-submission-ready target. Data pipeline complete (2026-04-13); regression stages pending.
**Data:** IFLS waves 1–5 (restricted-use; see `data/README.md`)

---

## Abstract

This study examines the long-term effects of early-life public health interventions on children's human capital in Indonesia's village midwife program. Human capital is measured by cognitive, mathematical, health, and socioemotional outcomes for the birth cohort aged 15–30, with a focus on early exposure before age three. The findings suggest that the village midwife program had an impact on health performance, particularly among girls. No significant impact on cognition, math, and socioemotional development. The findings highlight the significance of the timing of early-life health intervention programs on shaping later-life outcomes.

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

## Data-construction pipeline

| Stage | Script | Purpose |
|-------|--------|---------|
| 01 | `01_load_ifls.R` | Raw `.dta` → typed `.rds` cache under `data/intermediate/raw/` |
| 02 | `02_birth_roster.R` | Individual panel (cohorts 1984–1999), three `commid_birth` variants, mother covariates, `dead`/`mover*` as columns |
| 03 | `03_midwife_rollout.R` | SAR and PKK `start_year` per commid93 + rollout comparison diagnostic + other-facilities |
| 04 | `04_outcomes.R` | W5 outcomes (health, adult cognition via `w_abil`, CESD-10, Big-5), raw items preserved wide |
| 05 | `05_analysis_frame.R` | Merge + exposure-timing variables × 3 birthplace variants + diagnostics + audit vs legacy `df.rds` |

Library helpers live in `code/R/lib/`: `paths.R` (file catalog), `io.R` (raw read + cache), `harmonize.R` (commid anchoring), `education.R` (IFLS level+grade → years of schooling), `birthplace.R` (three `commid_birth` variants), `audit.R` (legacy-vs-new reconciliation).


## Repository layout


- `code/R/` — pipeline scripts (01..05) + `lib/` helpers
- `data/intermediate/` — gitignored; regenerable
- `old/` — frozen pre-rebuild project, read-only reference only

## Data access

IFLS microdata is restricted-use and is **not** committed to this repo. It is referenced by path from `$IFLS_DATA_ROOT`. Raw `.dta` files are never copied into the repo; only regenerable `.rds` intermediates under `data/intermediate/`, all gitignored. See `data/README.md`.

## License

TBD.
