#' Project setup: dependencies, paths and function loading
#'
#' Sourced at the top of every script in `analysis/`. Loads the packages the
#' project actually uses, defines the output directories, and sources the rest
#' of `R/`. Nothing here has an effect beyond attaching packages, creating
#' output directories and defining functions.
#' pkgbuild and pkgload are not called anywhere in this project, so renv::status()
#' reports them as unused. They are recorded deliberately: odin loads them to
#' compile its models to C, and without them it falls back to interpreted R
#' without saying so.

library(here)

library(odin)
library(TiPS)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpattern)
library(patchwork)
library(cowplot)
library(scales)

library(nleqslv)
library(lhs)
library(sensitivity)
library(rpart)
library(rpart.plot)
library(parallel)

options(ggpattern_use_R4.1_features = FALSE)

# ---------------------------------------------------------------------------
# Output locations
#
# All paths are resolved from the project root by `here`, so scripts run
# identically from RStudio, from `Rscript analysis/01_illustrative_example.R`,
# and from a cluster job, with no working-directory assumptions.
# ---------------------------------------------------------------------------

dir_data <- here::here("data")
dir_figures <- here::here("figures")

dir.create(dir_data, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

#' Build a path inside `data`
#' @param ... Path components, passed to `file.path()`.
#' @return A character path.
data_path <- function(...) file.path(dir_data, ...)

#' Build a path inside `figures`
#' @param ... Path components, passed to `file.path()`.
#' @return A character path.
figure_path <- function(...) file.path(dir_figures, ...)

# ---------------------------------------------------------------------------
# Function definitions
#
# `models_ode.R` and `floquet.R` compile three odin models between them, which
# takes a few seconds on the first call of an R session. Order matters only in
# that `simulate.R`, `batch.R` and `floquet.R` call into the files listed before
# them.
# ---------------------------------------------------------------------------

source(here::here("R", "interventions.R"))
source(here::here("R", "models_ode.R"))
source(here::here("R", "model_stochastic.R"))
source(here::here("R", "equilibrium.R"))
source(here::here("R", "simulate.R"))
source(here::here("R", "metrics.R"))
source(here::here("R", "batch.R"))
source(here::here("R", "plots.R"))
source(here::here("R", "decision_tree.R"))
source(here::here("R", "floquet.R"))
source(here::here("R", "tables.R"))
source(here::here("R", "case_studies.R"))
