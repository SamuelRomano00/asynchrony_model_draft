#' Deterministic metapopulation models
#'
#' Two odin model definitions and their runners:
#'  - a connected SIS model, the workhorse of the analysis;
#'  - a multi-patch Ross-Macdonald model (Ruktanonchai et al. 2016), used to
#'    check that the asynchrony effect is not an artefact of the SIS structure.
#'
#' Both share the same mobility convention, the same intervention format and
#' the same output shape, so the metric and plotting code is model-agnostic.

# ---------------------------------------------------------------------------
# Connected SIS model
# ---------------------------------------------------------------------------

#' Compiled generator for the connected SIS model
#'
#' An odin model definition, not a function: instantiate it with
#' `connectedSIS_deterministic_generator$new(user = parameters)` and integrate
#' with `$run(t)`. [run_SIS_deterministic()] wraps both steps.
#'
#' State variables, with `n` patches inferred from `length(N)`:
#'   `x[i]` prevalence in patch i, in \[0, 1\];
#'   `z[i]` cumulative reported cases in patch i.
connectedSIS_deterministic_generator <- odin::odin({

  ## Equations
  # Patch-i residents acquire infection at a rate driven by the force of
  # infection they meet across every patch they visit, and recover at rate r
  # straight back into the susceptible pool.
  deriv(x[]) <- r * ((1 - x[i]) * sum(pwlk[i, ]) - x[i])
  # Cumulative reported cases: true new infections scaled by the reporting
  # rate. No recovery term, this is an accumulator.
  deriv(z[]) <- r * (1 - x[i]) * rho_t[i] * sum(pwlkN[i, ])

  ## Initial conditions
  initial(x[]) <- x0[i]
  initial(z[]) <- z0[i]

  ## Inputs, supplied through `$new(user = parameters)`
  N[] <- user()       # population size per patch
  omega[, ] <- user() # vector-control efficacy schedule, [time x patch], in [0, 0.5]
  rho[, ] <- user()   # reporting-rate schedule, [time x patch]
  tt[] <- user()      # daily time grid the schedules are defined on
  p[, ] <- user()     # mobility matrix, p[i,j] = time share of patch-i residents in patch j
  x0[] <- user()
  z0[] <- user()
  R0[] <- user()      # basic reproduction number per patch
  r <- user()         # recovery rate, shared across patches

  ## Time-varying inputs
  # The solver evaluates the right-hand side at arbitrary intermediate times,
  # not only on whole days, so the daily schedules are turned into continuous
  # functions of time. The inputs are already at daily resolution, so linear
  # interpolation is effectively exact for what is conceptually a step function.
  omega_t[] <- interpolate(tt, omega, "linear")
  rho_t[] <- interpolate(tt, rho, "linear")

  ## Force of infection, i.e. the mobility coupling between patches
  pN[, ] <- p[i, j] * N[i]        # patch-i residents present in patch j
  pxN[, ] <- pN[i, j] * x[i]      # of those, how many are infected
  # Prevalence actually met in patch j: infected people present there, from any
  # patch of origin, over everyone present there. Residents and visitors mix.
  K[] <- sum(pxN[, i]) / sum(pN[, i])
  # Contribution to patch i's force of infection from time spent in patch j,
  # reduced by the vector control in force in patch j at that moment.
  pwlk[, ] <- p[i, j] * (1 - omega_t[j]) * R0[j] * K[j]
  pwlkN[, ] <- pwlk[i, j] * N[i]  # same, in case counts rather than per capita

  ## Dimensions
  dim(N) <- user()
  n <- length(N)

  dim(omega_t) <- n
  dim(rho_t) <- n
  dim(tt) <- user()
  dim(omega) <- c(length(tt), length(omega_t))
  dim(rho) <- c(length(tt), length(rho_t))
  dim(x) <- n
  dim(z) <- n
  dim(x0) <- n
  dim(z0) <- n
  dim(p) <- c(n, n)
  dim(pN) <- c(n, n)
  dim(pxN) <- c(n, n)
  dim(K) <- n
  dim(R0) <- n
  dim(pwlk) <- c(n, n)
  dim(pwlkN) <- c(n, n)

})

# ---------------------------------------------------------------------------
# Multi-patch Ross-Macdonald model
# ---------------------------------------------------------------------------

#' Compiled generator for the connected Ross-Macdonald model
#'
#' Human states `x[]` (prevalence) and `z[]` (cumulative cases) plus a mosquito
#' state `v[]` (sporozoite prevalence), with survival through the extrinsic
#' incubation period given by `exp(-mu * tau)`. Vector control reduces vectorial
#' capacity, the same lever as in the SIS model.
connectedRM_deterministic_generator <- odin::odin({

  ## Human force of infection, summed over every patch visited
  hFoI[, ] <- p[i, j] * m[j] * a[j] * b[j] * v[j]
  foi_x[] <- sum(hFoI[i, ])

  ## Human equations
  deriv(x[]) <- foi_x[i] * (1 - x[i]) - r[i] * x[i]
  deriv(z[]) <- rho * foi_x[i] * (1 - x[i]) * N[i]

  ## Mosquito equation
  deriv(v[]) <- a[i] * c[i] * (1 - omega_t[i]) * K[i] * (exp(-mu[i] * tau[i]) - v[i]) - mu[i] * v[i]

  ## Initial conditions
  initial(x[]) <- x0[i]
  initial(z[]) <- z0[i]
  initial(v[]) <- v0[i]

  ## Inputs
  N[] <- user()
  m[] <- user()       # mosquito-to-human ratio per patch
  p[, ] <- user()
  omega[, ] <- user() # vector-control efficacy schedule, [time x patch]
  tt[] <- user()
  x0[] <- user()
  z0[] <- user()
  v0[] <- user()
  a[] <- user()       # biting rate
  b[] <- user()       # transmission probability, mosquito to human
  c[] <- user()       # transmission probability, human to mosquito
  r[] <- user()       # human recovery rate
  mu[] <- user()      # mosquito mortality
  tau[] <- user()     # extrinsic incubation period
  rho <- user()       # reporting rate

  omega_t[] <- interpolate(tt, omega, "linear")

  ## Human prevalence met by the mosquitoes of patch i
  kappa_num[, ] <- p[i, j] * N[i] * x[i]
  kappa_den[, ] <- p[i, j] * N[i]
  K[] <- sum(kappa_num[, i]) / sum(kappa_den[, i])

  ## Dimensions
  dim(N) <- user()
  n <- length(N)
  dim(m) <- n
  dim(tt) <- user()
  dim(omega) <- c(length(tt), n)
  dim(omega_t) <- n
  dim(x) <- n
  dim(v) <- n
  dim(z) <- n
  dim(x0) <- n
  dim(v0) <- n
  dim(z0) <- n
  dim(a) <- n
  dim(b) <- n
  dim(c) <- n
  dim(r) <- n
  dim(mu) <- n
  dim(tau) <- n
  dim(p) <- c(n, n)
  dim(hFoI) <- c(n, n)
  dim(foi_x) <- n
  dim(kappa_num) <- c(n, n)
  dim(kappa_den) <- c(n, n)
  dim(K) <- n

})

# ---------------------------------------------------------------------------
# Shared post-processing
# ---------------------------------------------------------------------------

#' Turn cumulative reported cases into annual incidence
#'
#' Both odin models accumulate reported cases in `z[]`. This differences that
#' accumulator on year boundaries to get cases per year, in the long format the
#' metric and plotting functions expect.
#'
#' @param df Wide odin output, one row per day, with columns `t` and `z[i]`.
#' @param N Population size per patch, in the same order as the `z[]` columns.
#' @return A data frame with columns `t`, `variable` (`"z[1]"`, `"z[2]"`, ...),
#'   `value` (reported cases that year), `year`, `pop` and `incidence` (cases
#'   per 1000 population).
annual_incidence_from_cumulative <- function(df, N) {
  out <- df %>%
    select(t, starts_with("z")) %>%
    tidyr::pivot_longer(-t, names_to = "variable", values_to = "value") %>%
    arrange(variable, t) %>%
    mutate(year = ceiling(t / 365)) %>%
    filter(t %% 365 == 0) %>%
    group_by(variable) %>%
    # Differencing loses one value, so the t = 0 row has no earlier year to
    # subtract from. It is filled with year 1's incidence as a stand-in, which
    # keeps one row per year boundary per patch: the metric functions locate
    # each patch's block by integer position and rely on that layout.
    mutate(value = c(NA_real_, diff(value)),
           value = replace(value, 1, value[2])) %>%
    ungroup()

  pop <- data.frame(pop = N, variable = unique(out$variable))

  out %>%
    left_join(pop, by = "variable") %>%
    mutate(incidence = value * 1000 / pop)
}

#' Broadcast a constant schedule to a full daily matrix
#'
#' @param value Either a scalar/short vector, or an already-built
#'   `[n_days + 2 x n_patches]` daily matrix, which is returned unchanged.
#' @param n_days Simulation horizon, in days.
#' @param n_foci Number of patches.
#' @return A `[n_days + 2 x n_patches]` matrix.
as_daily_schedule <- function(value, n_days, n_foci) {
  if (is.matrix(value)) {
    return(value)
  }
  matrix(rep(value, n_days + 2), nrow = n_days + 2, ncol = n_foci, byrow = TRUE)
}

# ---------------------------------------------------------------------------
# Runners
# ---------------------------------------------------------------------------

#' Run the connected SIS model
#'
#' @param parameters Named list with `N`, `omega`, `rho`, `p`, `R0`, `r`, `x0`
#'   and `z0`. `omega` and `rho` may be given either as a constant or as a full
#'   daily matrix; `tt` is added automatically.
#' @param n_days Number of days to simulate.
#' @return A list with
#'   \describe{
#'     \item{simulation}{Raw wide odin output, one row per day.}
#'     \item{annual_incidence}{Long-format annual reported cases per patch, see
#'       [annual_incidence_from_cumulative()].}
#'     \item{prevalence}{Long-format daily prevalence per patch.}
#'   }
run_SIS_deterministic <- function(parameters, n_days) {

  t <- 0:n_days
  parameters$tt <- 0:(n_days + 1)   # one extra day of padding for interpolate()
  n_foci <- length(parameters$N)

  parameters$omega <- as_daily_schedule(parameters$omega, n_days, n_foci)
  parameters$rho <- as_daily_schedule(parameters$rho, n_days, n_foci)

  model <- connectedSIS_deterministic_generator$new(user = parameters)
  df <- as.data.frame(model$run(t))

  prevalence <- df %>%
    select(t, starts_with("x")) %>%
    tidyr::pivot_longer(-t, names_to = "variable", values_to = "value")

  list(simulation = df,
       annual_incidence = annual_incidence_from_cumulative(df, parameters$N),
       prevalence = prevalence)
}

#' Run the connected Ross-Macdonald model
#'
#' @param parameters Named list with `N`, `m`, `p`, `omega`, `x0`, `z0`, `v0`,
#'   `a`, `b`, `c`, `r`, `mu`, `tau` and `rho`. Entomological parameters given
#'   as a scalar are recycled across patches; `omega` may be a constant or a
#'   daily matrix.
#' @param n_days Number of days to simulate.
#' @return The same shape as [run_SIS_deterministic()], plus
#'   `mosquito_prevalence`.
run_RM_deterministic <- function(parameters, n_days) {

  t <- 0:n_days
  parameters$tt <- 0:(n_days + 1)
  n_foci <- length(parameters$N)

  parameters$omega <- as_daily_schedule(parameters$omega, n_days, n_foci)
  for (pname in c("a", "b", "c", "r", "mu", "tau")) {
    if (length(parameters[[pname]]) == 1) {
      parameters[[pname]] <- rep(parameters[[pname]], n_foci)
    }
  }

  model <- connectedRM_deterministic_generator$new(user = parameters)
  df <- as.data.frame(model$run(t))

  prevalence <- df %>%
    select(t, starts_with("x")) %>%
    tidyr::pivot_longer(-t, names_to = "variable", values_to = "value")

  mosquito_prevalence <- df %>%
    select(t, starts_with("v")) %>%
    tidyr::pivot_longer(-t, names_to = "variable", values_to = "value")

  list(simulation = df,
       annual_incidence = annual_incidence_from_cumulative(df, parameters$N),
       prevalence = prevalence,
       mosquito_prevalence = mosquito_prevalence)
}
