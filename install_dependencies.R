# Install the packages this project needs.
#
# Prefer `renv::restore()` if a renv.lock is present: it reinstalls the exact
# versions the results were produced with. Use this script only to bootstrap a
# fresh environment before creating that lockfile.
#
#   Rscript install_dependencies.R

packages <- c(
  # Models
  "odin",        # compiles the ODE systems
  "pkgbuild",    # lets odin compile to C rather than falling back to R
  "TiPS",        # exact stochastic simulation

  # Data handling and figures
  "dplyr", "tidyr", "tibble",
  "ggplot2", "ggpattern", "patchwork", "cowplot", "scales",

  # Numerics and experimental design
  "nleqslv",     # nonlinear solver, used for the equilibria
  "lhs",         # Latin hypercube designs
  "sensitivity", # Sobol indices
  "rpart", "rpart.plot",  # classification trees

  # Input/output
  "here"
)

missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) {
  install.packages(missing)
} else {
  message("All dependencies are already installed.")
}

# `parallel` ships with R and needs no installation.
# `odin` and `TiPS` compile C/C++ code: a working toolchain is required
# (Rtools on Windows, Xcode command line tools on macOS, r-base-dev on Linux).
