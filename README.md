# Assessing the impact of absence of coordination in malaria intervention strategies

Code and analysis scripts for __Assessing the impact of absence of coordination in malaria intervention strategies: a modelling study__  
Younes Iggidr, Nick Ruktanonchai, Bilal Benhana, Valérian Turbé, Billy Bauzile, Abigail Ward, Justin Cohen, Emilie Pothin, Clara Champagne.

The core model is a two-patch SIS metapopulation with Lagrangian mobility,
derived from Ruktanonchai et al. (2016). Interventions multiply the transmission
rate by `1 - omega_i`, alternating `t_I` years on and `t_I` years off, after one
full pre-intervention year at endemic equilibrium. In the synchronous scenario
both areas share a schedule; in the asynchronous one area 2 is shifted by one
period. Area 1 receives an identical schedule in both, so any change measured
there is attributable solely to the timing in area 2.

Three model variants are used, all sharing the same schedules, mobility
convention and metrics:

| Variant | Where it is used |
|---|---|
| SIS, deterministic (odin) | Sensitivity analysis, intervention-duration comparison, negative-AIG cases |
| SIS, stochastic CTMC (TiPS) | Illustrative example and case studies, where transmission interruption matters |
| Ross-Macdonald, deterministic (odin) | Check that the effect survives explicit vector-host dynamics |

## Metrics

A **cycle** is one "on" phase followed by one "off" phase, i.e. `2 * t_I` years.
The study window is 3 cycles, i.e. `6 * t_I` years, starting with the first year
under intervention:

| Metric | Definition |
|---|---|
| `AIG` | `sum_t ( I_async(t) - I_sync(t) )`, in cases per 1000 people |
| `AIG_year` | `AIG / (6 * t_I)` |
| `AIGR` | `AIG / sum_t I_sync(t)`, a relative excess |
| `AIGR_year` | `AIGR / (6 * t_I)` |

Each has an overall version and per-area versions (`_area1`, `_area2`). Most
analyses report area 1, the counterfactual area.

Transmission is considered interrupted in a stochastic replicate once its annual
incidence reaches zero and stays there, and in a scenario once at least 95% of
the 1000 replicates have done so.

## Parameter conventions

Every parameter in the code is on the paper's scale. There is no internal
complement to keep track of.

| Paper | Code | Range | Meaning |
|---|---|---|---|
| `X_i` | `x0`, `x[i]` | [0, 1] | Prevalence in area `i` |
| `Z_i` | `z0`, `z[i]` | — | Cumulative reported cases |
| `p_{1,2}` | `p_12` | (0, 0.5) | Share of area-1 residents' time spent in area 2 |
| `p_{2,1}` | `p_21` | (0, 0.5) | Share of area-2 residents' time spent in area 1 |
| `omega_i` | `omega_1`, `omega_2` | (0, 0.5) | Intervention effectiveness while on; 0 = no control |
| `R_{0,i}` | `R0_1`, `R0_2` | (0.9, 2.2) | Basic reproduction number |
| `R_{C,i}` | `RC_1`, `RC_2` | derived | `R0_i * (1 - omega_i)` |
| `r` | `r` | (1/200, 1/60) day⁻¹ | Recovery rate; the designs carry `rinv = 1/r` in days |
| `t_I` | `time_intervention` | {0.5, 1, ..., 5} years | Duration of one intervention |
| `kappa_j` | `K[j]` | — | Prevalence met in area `j`, mixing residents and visitors |

`mobility_matrix(p_12, p_21)` builds the matrix passed to the models, so the
diagonal (time spent at home) is derived rather than specified.

## Repository layout

```
R/                    functions, no side effects, never run directly
  setup.R             dependencies, output paths, sources everything below
  interventions.R     mobility matrix, on/off schedules, daily expansion
  models_ode.R        odin SIS and Ross-Macdonald models and their runners
  model_stochastic.R  TiPS simulator and its runner
  equilibrium.R       endemic equilibria and calibration
  simulate.R          synchronous vs asynchronous runs of each model
  metrics.R           AIG and AIGR, interruption curves
  batch.R             parallel batch loops over parameter designs, and the
                      loader every analysis script reads the database through
  plots.R             figures
  decision_tree.R     CART fitting, bootstrap stability, and their figures
  floquet.R           asymptotic stability under periodic schedules
  case_studies.R      the four case-study parameter sets, shared by 08 and 09
  tables.R            publication tables

analysis/             one script per figure or table, run in order
data/                 everything the scripts write; the derived CSVs needed by
                      05, 06 and 07 are committed
figures/              all figures
inst/hpc/             cluster submission template for 04
```

## Installation

Requires a C toolchain: Rtools on Windows, the Xcode command line tools on
macOS, `r-base-dev` on Linux. Without it, odin falls back to generating models
in R rather than C. The results are identical but the simulations are
substantially slower, and the message `Generating model in r` at startup is the
only sign.

```r
renv::restore()
```

`renv.lock` pins R 4.4.1 and 107 packages, the versions the published results
were produced with. To bootstrap an environment without renv:

```r
source("install_dependencies.R")
```

`sensitivity` pulls in a large dependency tree (`dtwclust`, `shiny`, `sf` and
their own dependencies), none of which this project calls directly. On a managed
HPC system, `sf` and `s2` may need `gdal-bin` and `libabsl-dev`: load the
corresponding modules if `renv::restore()` reports them as missing. On sciCORE
the Posit Package Manager binaries are incompatible, so set
`RENV_CONFIG_PPM_ENABLED=FALSE` and let renv compile from source.

All paths are resolved from the project root by `here`, so scripts run
identically from RStudio (open `asynchrony.Rproj`), from the command line
(`Rscript analysis/01_illustrative_example.R`), and from a cluster job. No
script changes the working directory.

## Reproducing the results

Scripts 05, 06 and 07 read what 04 writes; the rest are independent. The first
script of a session compiles the odin models, which adds a few seconds.

| Script | Paper element | Outputs | Runtime |
|---|---|---|---|
| `01_illustrative_example.R` | Figure 2, Table 2, Figure S1 | `figure2_illustrative_example.png`, `figureS1_interruption.png`, `table2_illustrative_example.csv` | ~10 min |
| `02_illustrative_example_vector_model.R` | Figure S2 | `figureS2_sis_vs_ross_macdonald.png` | seconds |
| `03_intervention_duration.R` | Figure S5 | `figureS5_intervention_duration.png` | ~1 min |
| `04_sensitivity_simulations.R` | Table S2; inputs to 05, 06 and 07 | `df_simulations.csv`, `sobol_*.csv` | hours, cluster |
| `05_sensitivity_figures.R` | Figure 3, Figures S3, S4, S6, Table S1 | `figure3_sensitivity_AIG.png`, `figureS3_sensitivity_AIGR.png`, `figureS4_near_elimination.png`, `figureS6_recovery_duration_heatmap.png`, `tableS1_metric_quantiles.csv` | ~1 min |
| `06_decision_trees.R` | Figure 4, Figures S7, S8, S9, S10 | `figure4_decision_tree.png`, `figureS7_decision_tree_thresholds.png`, `figureS8_bootstrap_stability.png`, `figureS9_bootstrap_tree_size.png`, `figureS10_decision_tree_AIGR.png` | ~5 min |
| `07_negative_aig_cases.R` | Figure S11, Table S3 | `figureS11_negative_aig_cases.png` | seconds |
| `08_case_studies.R` | Figure 5, Figure S12, Table S4 (metrics) | `figure5_case_studies.png`, `figureS12_case_studies_interruption.png`, `tableS4_case_studies.csv` | ~40 min |
| `09_floquet_exponents.R` | Table S4 (Floquet exponents) | `tableS4_floquet_exponents.csv` | seconds |

Figure 1 is a schematic and has no code. Table S2 is `sobol_AIG_S2.csv`,
transcribed as a lower-triangular matrix; the parameter order in the paper
differs from the column order of the CSV.

### The sensitivity analysis

`04` builds two things. The simulation database is 10,000 Latin hypercube
parameter sets, two ODE solves each. The Sobol design is derived from the same
10,000 draws but has `10,000 * (2k + 2) = 180,000` rows for `k = 8` parameters,
again two solves each: that is what needs a cluster, and why
`AIG_AIGR_computation()` checkpoints. See
`inst/hpc/sensitivity_simulations.sbatch`. The number of workers is detected
from `SLURM_CPUS_PER_TASK`, so setting `--cpus-per-task` is enough and nothing
in the R code has to change.

A parameter set is **viable** when both patches have a positive endemic
equilibrium. Where the system-level reproduction number is at most 1, the
disease-free state is the only equilibrium, both trajectories are identically
zero, and the metrics degenerate asymmetrically: `AIG = 0` is a finite number
while `AIGR = 0/0` is undefined. Two things follow.

- `04` restricts the Sobol design to draws whose **every** one of the `2k + 2`
  blocks is viable, not just the two source samples, and rebuilds the design on
  the surviving indices so the paired structure the estimator relies on is
  preserved.
- `05`, `06` and `07` read the database through `load_simulation_database()`,
  which drops the non-viable rows.

The Latin hypercube designs use fixed seeds (`SEED_LHS1`, `SEED_LHS2` at the top
of the script), so the designs and the database are reproducible. Changing them
changes every downstream figure.

### Committed data

To let `05`, `06` and `07` be re-run without repeating `04`, `data/` contains:

```
df_simulations.csv        one row per Latin hypercube parameter set
sobol_AIG.csv             first-order and total indices, AIG per year
sobol_AIGR.csv            the same for AIGR per year
sobol_AIGR_rank.csv       the same on ranks
sobol_AIG_S2.csv          second-order indices (Supplementary Table S2)
sobol_AIGR_S2.csv         the same for AIGR per year
sobol_AIGR_rank_S2.csv    the same on ranks
```

`LHS1.csv`, `LHS2.csv` and `sobol_design_Y.csv` are not committed: the first two
are regenerated deterministically from their seeds, and the third is the
180,000-row evaluation of the full Sobol design, needed only inside `04`.

## Licence

GNU General Public License v3.0 — see [LICENSE](LICENSE).