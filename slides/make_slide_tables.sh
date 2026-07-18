#!/usr/bin/env bash
# =====================================================================
# Derive slide-specific table fragments from the paper's auto-generated
# tables in paper/tables/. The paper keeps "Yes"/"No" spec rows; the deck
# shows a checkmark for Yes and a blank cell for No.
#
# Re-run this whenever the R pipeline regenerates paper/tables/:
#   Rscript code/R/00_run_all.R && slides/make_slide_tables.sh
#
# Never edit slides/tables/*.tex by hand -- edits are overwritten here.
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

SRC=../paper/tables
OUT=tables
mkdir -p "$OUT"

# Spec rows: "Yes" -> \cmark (defined in slides.tex), "No" -> empty cell.
for f in main_results main_results_by_sex sibling_fe; do
  {
    echo "% Derived from paper/tables/$f.tex by slides/make_slide_tables.sh -- do not edit."
    sed -e 's/& Yes/\& \\cmark/g' -e 's/& No/\& /g' "$SRC/$f.tex"
  } > "$OUT/$f.tex"
done

# Summary statistics ships as a full float with a caption; the deck needs
# the bare tabular so it can live inside a frame.
{
  echo "% Derived from paper/tables/table1_descriptives.tex by slides/make_slide_tables.sh -- do not edit."
  sed -n '/\\begin{tabular}/,/\\end{tabular}/p' "$SRC/table1_descriptives.tex"
} > "$OUT/table1_descriptives.tex"

# Fertility-selection check: the deck cites these two numbers in a bullet,
# where a tabular will not fit. Expose them as macros so no coefficient is
# ever typed by hand.
{
  echo "% Derived from paper/tables/sibling_fertility_check.tex by slides/make_slide_tables.sh -- do not edit."
  awk -F'&' '/Exposed during reproductive years/ {
    est = $2; se = $3
    gsub(/[^0-9.]/, "", est)
    gsub(/[^0-9.]/, "", se)
    printf "\\newcommand{\\CoefFert}{%.3f}\n", est
    printf "\\newcommand{\\SEFert}{%.3f}\n", se
  }' "$SRC/sibling_fertility_check.tex"
} > "$OUT/fertility_macros.tex"

echo "Wrote slide tables to slides/$OUT/:"
ls -1 "$OUT"
