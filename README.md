# Assessing the impact of absence of coordination in malaria intervention strategies

Code and analysis scripts for *Assessing the impact of absence of coordination
in malaria intervention strategies: a modelling study*.

<!-- TODO: authors, journal, DOI, preprint link, and a CITATION.cff once the
     repository has a permanent home. -->

Malaria control is increasingly tailored at subnational scales, but neighbouring
areas stay connected by human mobility. When two connected areas run the same
intervention on different schedules, each re-seeds the other during its own
"off" period. This repository quantifies that penalty with the **Asynchrony
Induced Growth (AIG)**: the excess cases caused by asynchronous timing, relative
to deploying exactly the same intervention in phase.

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

Over a study window of 3 intervention cycles, i.e. `6 * t_I` years:

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
  batch.R             parallel batch loops over parameter designs
  plots.R             figures
  decision_tree.R     CART fitting, bootstrap stability, and their figures
  floquet.R           asymptotic stability under periodic schedules
  case_studies.R      the four case-study parameter sets, shared by 08 and 09
  tables.R            publication tables

analysis/             one script per figure or table, run in order
data/         everything the scripts write
figures/              all figures
inst/hpc/             cluster submission template for 04
```

## Installation

Requires R (>= 4.1) and a C/C++ toolchain, since `odin` and `TiPS` compile their
models: Rtools on Windows, the Xcode command line tools on macOS, `r-base-dev`
on Linux.

```r
# Preferred, once renv.lock exists:
renv::restore()

# Or, to bootstrap a fresh environment:
source("install_dependencies.R")
```

<!-- TODO: create the lockfile on your machine, then commit it:
     renv::init(); renv::snapshot()
     A hand-written lockfile would not reflect the versions the results were
     actually produced with, so it has to be generated rather than supplied. -->

All paths are resolved from the project root by `here`, so scripts run
identically from RStudio (open `asynchrony.Rproj`), from the command line
(`Rscript analysis/01_illustrative_example.R`), and from a cluster job. No
script changes the working directory.

## Reproducing the results

Scripts 05, 06 and 07 read what 04 writes; the rest are independent. The first
script of a session compiles both odin models, which adds a few seconds.

| Script | Paper element | Outputs | Runtime |
|---|---|---|---|
| `01_illustrative_example.R` | Figure 2, Table 2, Figure S1 | `figure2_illustrative_example.png`, `figureS1_interruption.png`, `table2_illustrative_example.csv` | ~10 min |
| `02_illustrative_example_vector_model.R` | Figure S2 | `figureS2_sis_vs_ross_macdonald.png` | seconds |
| `03_intervention_duration.R` | Figure S6 | `figureS6_intervention_duration.png` | ~1 min |
| `04_sensitivity_simulations.R` | Table S2; inputs to 05 and 07 | `df_simulations.csv`, `sobol_*.csv` | hours, cluster |
| `05_sensitivity_figures.R` | Figure 3, Figures S3, S5, S7, Table S1 | `figure3_sensitivity_AIG.png`, `figureS3_sensitivity_AIGR.png`, `figureS5_near_elimination.png`, `figureS7_recovery_duration_heatmap.png`, `tableS1_metric_quantiles.csv` | ~1 min |
| `06_decision_trees.R` | Figure 4, Figures S8, S9, S10, S11 | `figure4_decision_tree.png`, `figureS8_decision_tree_thresholds.png`, `figureS9_bootstrap_stability.png`, `figureS10_bootstrap_tree_size.png`, `figureS11_decision_tree_AIGR.png` | ~5 min |
| `07_negative_aig_cases.R` | Figure S12, Table S3 | `figureS12_negative_aig_cases.png` | seconds |
| `08_case_studies.R` | Figure 5, Figure S13, Table S4 (metrics) | `figure5_case_studies.png`, `figureS13_case_studies_interruption.png`, `tableS4_case_studies.csv` | ~40 min |
| `09_floquet_exponents.R` | Table S4 (Floquet exponents) | `tableS4_floquet_exponents.csv` | seconds |

Figure 1 is a schematic and has no code.

Every analysis in the paper has code here. One caveat: the Figure S7 heatmap in
`05` was reimplemented from `df_simulations.csv` rather than carried over, so
check it against the published version.

### The sensitivity analysis

`04` runs 10,000 parameter sets for the simulation database and a Sobol design
of the same size, two ODE solves each. It is the only script that needs a
cluster; see `inst/hpc/sensitivity_simulations.sbatch`. The number of workers is
detected from `SLURM_CPUS_PER_TASK`, so setting `--cpus-per-task` is enough and
nothing in the R code has to change.

Parameter sets whose equilibrium prevalence solves to zero, i.e. `R0 < 1` in
both areas or below 1 in one and barely above in the other, are excluded from
the Sobol indices; the paper reports 163 such sets out of 10,000.

The Latin hypercube designs and the Sobol bootstrap use fixed seeds
(`SEED_LHS1`, `SEED_LHS2`, `SEED_SOBOL` at the top of the script), so the whole
chain is reproducible. Changing them changes every downstream figure.

The derived CSVs are committed, so `05` and `07` can be re-run without repeating
`04`.

## Notes for readers of the code

- **Two different setups inside the sensitivity analysis.** The simulation
  database starts interventions on day 365 over a 30,000-day horizon; the Sobol
  design uses day 730 over 35,000 days. The difference is cosmetic and returns
  identical metrics, because the model sits at the endemic equilibrium until the
  intervention starts and `compute_metrics()` locates its window from those same
  two arguments. Both functions now take them as arguments if you want to align
  them anyway.
- **AIGR near elimination.** When a synchronous run has already eliminated
  transmission inside the metric window, the AIGR denominator is zero and the
  ratio is undefined, while the AIG stays well defined. The stochastic summaries
  count and exclude those replicates rather than letting them blank out a whole
  row. `04` additionally writes rank-based AIGR indices
  (`sobol_AIGR_rank.csv`), insensitive to that tail; `05` plots the raw-scale
  ones and can be pointed at the rank version in one line.
- **Metric window indexing.** `compute_metrics()` locates each area's block of
  years by integer position, not by label. Both the deterministic and the
  stochastic runners therefore emit the same yearly layout: one entry per
  complete year, preceded by a `t = 0` placeholder.
- **Case management is not modelled.** The `alpha` term present in earlier
  versions of the code has been removed: it was fixed at 0 in every analysis in
  the paper.
- **The Floquet window starts on a cycle boundary.** Multipliers are the same
  wherever the window sits inside the periodic regime, but it has to sit inside
  it, and the schedule only becomes periodic once the interventions begin: the
  simulations open with a pre-intervention year so the endemic equilibrium is
  visible on the figures. `09` drops those days before integrating, and
  `check_schedule_periodicity()` refuses any window that is not a period.
- **Randomness in the trees.** `06` depends on the RNG twice: `rpart`'s
  10-fold cross-validation assigns folds at random, and the stability analysis
  resamples the database 500 times. One `set.seed(123)` at the top of the
  script covers both.
- **Parallelism.** The batch loops fork with `mclapply()`, so workers inherit
  the already-compiled odin model. On Windows they fall back to a single core.

## References

- Ruktanonchai N.W. et al. (2016) Identifying malaria transmission foci for
  elimination using human mobility data. *PLOS Computational Biology* 12(4),
  e1004846.
- Cosner C. et al. (2009) The effects of human movement on the persistence of
  vector-borne diseases. *Journal of Theoretical Biology* 258(4), 550-560.
- Zhao X.-Q. & Jing Z.-J. (1996) Global asymptotic behavior in some cooperative
  systems of functional differential equations. *Canadian Applied Mathematics
  Quarterly* 4(4).
- Wang W. & Zhao X.-Q. (2008) Threshold dynamics for compartmental epidemic
  models in periodic environments. *Journal of Dynamics and Differential
  Equations* 20(3), 699-717.
- FitzJohn R. et al. (2025) odin: ODE generation and integration. R package.
- Danesh G. et al. (2023) TiPS: rapidly simulating trajectories and phylogenies
  from compartmental models. *Methods in Ecology and Evolution* 14(2), 487-495.

## Licence

GNU General Public License v3.0 — see [LICENSE](LICENSE).
