# Early-Life Exposure to Indonesia's Village Midwife Program and Child Human Capital

**Author:** Wanzi Wang
**Status:** In rebuild — journal-submission-ready target
**Data:** IFLS waves 1–5 (restricted-use; see `data/README.md`)

---

## Abstract

This study examines the long-term effects of early-life public health interventions on children's human capital in Indonesia's village midwife program. Human capital is measured by cognitive, mathematical, health, and socioemotional outcomes for the birth cohort aged 15–30, with a focus on early exposure before age three. The findings suggest that the village midwife program had an impact on health performance, particularly among girls. No significant impact on cognition, math, and socioemotional development. The findings highlight the significance of the timing of early-life health intervention programs on shaping later-life outcomes.

---

## How to run

```bash
# 1. Make sure IFLS raw data is accessible
export IFLS_DATA_ROOT="/Users/wwz/Academic/Thesis_Data_IFLS/old/IFLS"

# 2. Run the analysis pipeline (produces figures/tables under paper/)
Rscript code/R/00_run_all.R

# 3. Build the paper and slides
cd paper  && latexmk -xelatex -outdir=../output/paper paper.tex
cd slides && latexmk -xelatex -outdir=../output/slides slides.tex
```

## Repository layout


## Data access

IFLS microdata is restricted-use and is **not** committed to this repo. It is referenced by path from `$IFLS_DATA_ROOT`. See `data/README.md`.

## License

TBD.
