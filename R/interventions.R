#' Two-patch mobility matrix
#'
#' Builds the mobility matrix used by every model in the project from the two
#' off-diagonal mobility fractions.
#'
#' @param p_12 Fraction of patch-1 residents' time spent in patch 2, in
#'   \[0, 0.5\].
#' @param p_21 Fraction of patch-2 residents' time spent in patch 1, in
#'   \[0, 0.5\].
#' @return A 2x2 matrix `p` with `p[i, j]` the fraction of patch-i residents'
#'   time spent in patch j. Rows sum to 1, and the diagonal is the time spent
#'   at home.
mobility_matrix <- function(p_12, p_21) {
  matrix(c(1 - p_12, p_12,
           p_21,     1 - p_21),
         nrow = 2, byrow = TRUE)
}

#' Breakpoints of an alternating on/off intervention schedule
#'
#' Returns the period boundaries of a schedule that alternates "intervention
#' on" and "intervention off" every `length_intervention` years, starting at
#' `start_interv`, and covering the whole simulation horizon.
#'
#' @param n_days Simulation horizon, in days.
#' @param start_interv Day the first intervention period begins.
#' @param length_intervention Duration of each on/off period, in years.
#' @return A list with
#'   \describe{
#'     \item{times}{Numeric vector of period boundaries, starting at 0 and
#'       ending at `n_days + 2`. Has one more element than there are periods.}
#'     \item{on_first}{Logical vector, one value per period *after* the initial
#'       pre-intervention period, `TRUE` where the intervention is active in
#'       the synchronous scenario.}
#'   }
intervention_breakpoints <- function(n_days, start_interv, length_intervention) {
  # The number of on/off pairs is derived from n_days rather than fixed, so the
  # alternation always spans exactly the simulated horizon: too few periods
  # leaves the end of the run with no intervention at all, too many pushes the
  # final breakpoint below its predecessor and makes `times` non-monotonic.
  n_pairs <- ceiling(ceiling((n_days - start_interv) / (365 * length_intervention)) / 2)
  bps <- start_interv + seq(0, 2 * n_pairs) * 365 * length_intervention
  bps <- bps[bps < n_days + 2]

  list(times = c(0, bps, n_days + 2),
       on_first = seq_along(bps) %% 2 == 1)
}

#' Synchronous and asynchronous intervention schedules
#'
#' Builds the two period-level vector-control schedules compared throughout the
#' paper: both patches switching the intervention on and off together
#' (synchronous), and patch 2 phase-shifted by exactly one period
#' (asynchronous). Everything else about the two runs is identical, so any
#' difference in outcome is attributable to timing alone.
#'
#' @param on_first Logical vector from [intervention_breakpoints()].
#' @param omega Numeric vector `c(omega_1, omega_2)` of vector-control
#'   efficacies, in \[0, 0.5\]; 0 means no intervention.
#' @return A list with matrices `synchronous` and `asynchronous`, each with one
#'   row per period and one column per patch, holding the efficacy in force
#'   during that period. Both start with a leading 0 row for the period before
#'   `start_interv`.
intervention_schedules <- function(on_first, omega) {
  synchronous <- cbind(c(0, ifelse(on_first, omega[1], 0)),
                       c(0, ifelse(on_first, omega[2], 0)))
  asynchronous <- cbind(c(0, ifelse(on_first, omega[1], 0)),
                        c(0, ifelse(on_first, 0, omega[2])))

  list(synchronous = synchronous, asynchronous = asynchronous)
}

#' Expand a period-level schedule into a daily matrix
#'
#' The odin models interpolate their intervention inputs on a daily grid, so a
#' schedule expressed as one row per period has to be expanded to one row per
#' day before it can be passed to them.
#'
#' @param interventions Matrix with one row per period and one column per
#'   patch, e.g. an element of [intervention_schedules()]'s output.
#' @param interv_times Numeric vector of period boundaries; must start at 0 and
#'   have `nrow(interventions) + 1` elements.
#' @return A matrix with one row per day up to the final breakpoint and one
#'   column per patch.
create_timevarying_interventions <- function(interventions, interv_times) {
  if (interv_times[1] != 0) {
    stop("`interv_times` must start at 0.")
  }
  # Round the breakpoints, not the period lengths. Differences of a rounded
  # monotone sequence are whole numbers of days that telescope to exactly the
  # horizon, so no day is lost or recycled when `length_intervention` is
  # fractional (1.5, 2.5, ...).
  interv_times <- round(interv_times)

  myvect <- c()
  for (i in seq_len(length(interv_times) - 1)) {
    myvect <- c(myvect, rep(interventions[i, ], interv_times[i + 1] - interv_times[i]))
  }

  matrix(myvect,
         nrow = interv_times[length(interv_times)],
         ncol = ncol(interventions), byrow = TRUE)
}
