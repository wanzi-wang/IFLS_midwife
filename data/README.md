# Data

## Raw IFLS (restricted-use, NOT committed)

Raw Indonesian Family Life Survey (IFLS) waves 1–5 live outside the repo.

- **Default path:** `/Users/wwz/Academic/Thesis_Data_IFLS/old/IFLS/W1..W5`
- **Override:** set the environment variable `IFLS_DATA_ROOT`.

R scripts read from `IFLS_DATA_ROOT`. Do not copy these files into the repo — IFLS is under a restricted-data agreement.

## Intermediate data (regenerable, gitignored)

`data/intermediate/` holds R-generated `.rds` caches (birth roster, treatment assignment, merged analysis frame). Everything here is reproducible from `Rscript code/R/00_run_all.R` — safe to delete and rebuild.

## Legacy reference

The pre-rebuild (messy) project lives at `/Users/wwz/Academic/Thesis_Data_IFLS/old/`. Treat it as read-only. Canonical legacy scripts to mine for variable-construction logic:

- `old/Code/Birth_roster_full.Rmd`
- `old/Code/midwife.Rmd`
- `old/Code/cognition_birthroster.Rmd`
- `old/Code/Health.Rmd`
- `old/Code/Table.Rmd`

Main legacy analysis dataset: `old/Intermediate/Reg/df.rds`.
