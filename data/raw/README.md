# Raw inputs

`AIG_neg.xlsx` — the parameter sets from the 10,000 simulations of
`analysis/04_sensitivity_simulations.R` whose AIG in area 1 is negative. It
holds every column of `df_simulations.csv`, so `analysis/06_aig_negative_cases.R`
can read the initial conditions (`x0_1`, `x0_2`, `z0_1`, `z0_2`, `r`) directly
and reproduce the exact trajectories rather than nearby ones.

The file is derived from `data/derived/df_simulations.csv` by keeping the rows
with `AIG_area1 < 0`. Regenerate it after re-running `04` if the design or the
seeds change.
