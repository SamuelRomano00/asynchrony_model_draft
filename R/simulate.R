#' Synchronous versus asynchronous simulation
#'
#' The core experiment: run the same model twice under the same intervention,
#' once with both patches switching it on and off together, once with patch 2
#' phase-shifted by one period. Everything else is held fixed, so the
#' difference between the two runs isolates the effect of timing.

#' Simulate the SIS model under both intervention schedules
#'
#' @param n_days Simulation horizon, in days.
#' @param start_interv Day the first intervention period begins.
#' @param length_intervention Duration of each on/off period, in years.
#' @param N Numeric vector of population sizes per patch.
#' @param params Named list with
#'   \describe{
#'     \item{R0}{`c(R0_1, R0_2)`.}
#'     \item{p}{`c(p_12, p_21)`, mobility fractions in \[0, 0.5\].}
#'     \item{omega}{`c(omega_1, omega_2)`, vector-control efficacies in
#'       \[0, 0.5\], in force while the intervention is on.}
#'     \item{rho}{`c(rho_1, rho_2)`, reporting rates.}
#'     \item{x0}{`c(x0_1, x0_2)`, initial prevalence.}
#'     \item{z0}{`c(z0_1, z0_2)`, initial cumulative reported cases.}
#'     \item{r}{Recovery rate.}
#'   }
#' @return A list with `synchronous` and `asynchronous`, each the full output of
#'   [run_SIS_deterministic()].
simulate_sis <- function(n_days, start_interv, length_intervention, N, params) {

  breaks <- intervention_breakpoints(n_days, start_interv, length_intervention)
  schedules <- intervention_schedules(breaks$on_first, params$omega)

  base <- list(
    N = N,
    rho = params$rho,
    R0 = params$R0,
    p = mobility_matrix(params$p[1], params$p[2]),
    z0 = params$z0,
    x0 = params$x0,
    r = params$r
  )

  run_one <- function(schedule) {
    pars <- base
    pars$omega <- create_timevarying_interventions(schedule, breaks$times)
    run_SIS_deterministic(parameters = pars, n_days = n_days)
  }

  list(synchronous = run_one(schedules$synchronous),
       asynchronous = run_one(schedules$asynchronous))
}

#' Simulate the Ross-Macdonald model under both intervention schedules
#'
#' Same schedules and mobility convention as [simulate_sis()], with the
#' mosquito-to-human ratio derived from the target `R0` via [R0_to_m()].
#'
#' @inheritParams simulate_sis
#' @param params Named list with `R0`, `p`, `omega`, `x0`, `z0`, `v0`, `a`,
#'   `b`, `c_rate`, `r`, `mu`, `tau` and `rho` (a scalar here).
#' @return A list with `synchronous` and `asynchronous`, each the full output of
#'   [run_RM_deterministic()].
simulate_rm <- function(n_days, start_interv, length_intervention, N, params) {

  breaks <- intervention_breakpoints(n_days, start_interv, length_intervention)
  schedules <- intervention_schedules(breaks$on_first, params$omega)

  base <- list(
    N = N,
    m = R0_to_m(params$R0, params$a, params$b, params$c_rate,
                params$r, params$mu, params$tau),
    p = mobility_matrix(params$p[1], params$p[2]),
    z0 = params$z0,
    x0 = params$x0,
    v0 = params$v0,
    a = params$a,
    b = params$b,
    c = params$c_rate,
    r = params$r,
    mu = params$mu,
    tau = params$tau,
    rho = params$rho
  )

  run_one <- function(schedule) {
    pars <- base
    pars$omega <- create_timevarying_interventions(schedule, breaks$times)
    run_RM_deterministic(parameters = pars, n_days = n_days)
  }

  list(synchronous = run_one(schedules$synchronous),
       asynchronous = run_one(schedules$asynchronous))
}

#' Simulate the stochastic SIS model under both intervention schedules
#'
#' The stochastic analogue of [simulate_sis()]. The schedules are the same, but
#' are passed to TiPS as one value per period rather than expanded to a daily
#' matrix, since TiPS resolves the current period internally.
#'
#' @inheritParams simulate_sis
#' @param simulator Simulator from [build_SIS_stochastic_simulator()]. Build it
#'   once per session and pass it in; compilation is the slow part.
#' @param params Named list with `R0`, `p`, `omega`, `x0` and `r`. `z0` and
#'   `rho` are not used: the stochastic model counts infection events directly
#'   rather than accumulating a reporting-scaled compartment.
#' @param n_reps Number of replicates per scenario.
#' @param method,seed,verbose Passed through to [run_SIS_stochastic()].
#' @return A list with `synchronous` and `asynchronous`, each the full output of
#'   [run_SIS_stochastic()].
simulate_sis_stochastic <- function(simulator, n_days, start_interv, length_intervention,
                                    N, params, n_reps = 20,
                                    method = "exact", seed = 1, verbose = FALSE) {

  breaks <- intervention_breakpoints(n_days, start_interv, length_intervention)
  schedules <- intervention_schedules(breaks$on_first, params$omega)

  base <- list(
    N = N,
    R0 = params$R0,
    p = mobility_matrix(params$p[1], params$p[2]),
    r = params$r,
    x0 = params$x0
  )

  run_one <- function(schedule) {
    pars <- base
    pars$omega <- list(omega1 = schedule[, 1], omega2 = schedule[, 2])
    run_SIS_stochastic(simulator, pars, times = breaks$times,
                       n_reps = n_reps, method = method, seed = seed,
                       verbose = verbose)
  }

  list(synchronous = run_one(schedules$synchronous),
       asynchronous = run_one(schedules$asynchronous))
}
