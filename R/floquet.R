#' Floquet exponents of the disease-free equilibrium
#'
#' The AIG is a finite-horizon quantity: it says how many extra cases asynchrony
#' causes over a fixed number of cycles, not whether the epidemic ultimately
#' dies out. This complements it with the asymptotic answer. Linearising the SIS
#' model around the disease-free equilibrium under a periodic intervention
#' schedule gives a linear system with periodic coefficients, whose stability is
#' characterised by the dominant Floquet exponent: negative means convergence to
#' the disease-free state, positive means the epidemic persists.
#'
#' Following Appendix A, with `q[i,j] = p[i,j] (1 - omega_j) lambda_j / a_j`,
#' `b[i,j] = p[i,j] N_i`, `a_j = sum_i b[i,j]` and `lambda_i = r R0_i`, the
#' linearised system is `dX/dt = J(0,t) X(t)` with `J(0,t) = Q(t) B(t)^T - r I`.

#' Compiled generator for the fundamental solution matrix
#'
#' Integrates `dY/dt = J(0,t) Y(t)` from `Y = I`, so that `Y(T)` is the
#' monodromy matrix over one period.
#'
#' State variable `Y[i,j]`, with `n` patches inferred from `length(N)`.
floquet_generator <- odin::odin({

  ## Jacobian of the linearised system at the disease-free equilibrium.
  # Split into an explicit triple-indexed product because odin's sum() takes a
  # plain indexed array rather than a matrix expression.
  qb_prod[, , ] <- q[i, k] * b[j, k]
  J[, ] <- sum(qb_prod[i, j, ]) - r * Id[i, j]

  ## Fundamental solution matrix
  prod[, , ] <- J[i, k] * Y[k, j]
  deriv(Y[, ]) <- sum(prod[i, j, ])
  initial(Y[, ]) <- Id[i, j]

  ## Inputs
  N[] <- user()
  p[, ] <- user()
  R0[] <- user()
  r <- user()
  tt[] <- user()
  omega[, ] <- user()  # daily efficacy schedule, [time x patch], in [0, 0.5]
  Id[, ] <- user()     # identity matrix, used both in J and as Y(0)

  ## Intermediates
  omega_t[] <- interpolate(tt, omega, "linear")
  lambda[] <- r * R0[i]
  b[, ] <- p[i, j] * N[i]
  a[] <- sum(b[, i])
  q[, ] <- p[i, j] * (1 - omega_t[j]) * lambda[j] / a[j]

  ## Dimensions
  dim(N) <- user()
  n <- length(N)
  dim(tt) <- user()
  dim(p) <- c(n, n)
  dim(R0) <- n
  dim(omega) <- c(length(tt), n)
  dim(omega_t) <- n
  dim(Id) <- c(n, n)
  dim(lambda) <- n
  dim(b) <- c(n, n)
  dim(a) <- n
  dim(q) <- c(n, n)
  dim(J) <- c(n, n)
  dim(Y) <- c(n, n)
  dim(qb_prod) <- c(n, n, n)
  dim(prod) <- c(n, n, n)

})

#' Dominant Floquet exponent from a daily intervention schedule
#'
#' @param N Population sizes per patch.
#' @param p Mobility matrix, see [mobility_matrix()].
#' @param R0 Basic reproduction number per patch.
#' @param r Recovery rate.
#' @param omega_daily Daily efficacy schedule, one row per day, one column per
#'   patch. Row 1 must be the start of a full intervention cycle.
#' @param T_period Length of one full intervention cycle, in days.
#' @param n_out Number of reporting points used to integrate over the period.
#' @return A list with the `monodromy` matrix, its dominant multiplier `rho`,
#'   and the exponent `Lambda = log(rho) / T_period`.
compute_floquet_exponent <- function(N, p, R0, r, omega_daily, T_period, n_out = 2000) {

  n <- length(N)
  tt <- 0:(nrow(omega_daily) - 1)

  model <- floquet_generator$new(N = N, p = p, R0 = R0, r = r,
                                 tt = tt, omega = omega_daily, Id = diag(n))

  out <- model$run(seq(0, T_period, length.out = n_out))

  Y_cols <- grep("^Y", colnames(out))
  # odin flattens Y in column-major order, which is how matrix() refills it.
  monodromy <- matrix(out[nrow(out), Y_cols], nrow = n, ncol = n)

  rho <- max(Mod(eigen(monodromy, only.values = TRUE)$values))

  list(monodromy = monodromy, rho = rho, Lambda = log(rho) / T_period)
}

#' Floquet exponents under both intervention schedules
#'
#' The counterpart of [simulate_sis()] for asymptotic behaviour: same
#' parameters, same schedules, but the growth rate of the disease-free
#' equilibrium instead of a trajectory.
#'
#' @param length_intervention Years per on/off period.
#' @param N Population sizes per patch.
#' @param params Named list with `R0`, `p` (`c(p_12, p_21)`), `omega`
#'   (`c(omega_1, omega_2)`) and `r`. Initial conditions are not used: the
#'   system is linearised at the disease-free equilibrium.
#' @param start_interv Day the first intervention period begins.
#' @return A list with `synchronous` and `asynchronous`, each the output of
#'   [compute_floquet_exponent()].
#'
#' The integration window starts on a cycle boundary, not at day 0. Floquet
#' multipliers are the same wherever the window sits inside the periodic regime,
#' but the window has to lie inside it, and the schedule only becomes periodic
#' at `start_interv`: the simulations open with a pre-intervention year so that
#' the endemic equilibrium is visible on the figures. Those leading days are
#' dropped here. Keeping them would truncate area 2's asynchronous control,
#' whose first "on" phase begins half a cycle after `start_interv`, to
#' `T/2 - start_interv` days out of `T/2`, understating the intervention and
#' overstating its growth rate. [check_schedule_periodicity()] refuses any
#' window that is not a period, so that cannot pass silently.
floquet_exponents <- function(length_intervention, N, params, start_interv = 365) {

  T_period <- 2 * length_intervention * 365

  # Enough cycles for the periodicity check to have a second period to compare
  # against, with room to spare.
  n_days <- ceiling(start_interv + 4 * T_period)

  breaks <- intervention_breakpoints(n_days, start_interv, length_intervention)
  schedules <- intervention_schedules(breaks$on_first, params$omega)

  run_one <- function(schedule) {
    omega_daily <- create_timevarying_interventions(schedule, breaks$times)
    omega_daily <- omega_daily[-seq_len(round(start_interv)), , drop = FALSE]
    check_schedule_periodicity(omega_daily, T_period)
    compute_floquet_exponent(N = N,
                             p = mobility_matrix(params$p[1], params$p[2]),
                             R0 = params$R0, r = params$r,
                             omega_daily = omega_daily, T_period = T_period)
  }

  list(synchronous = run_one(schedules$synchronous),
       asynchronous = run_one(schedules$asynchronous))
}

#' Check that a daily schedule really repeats over one period
#'
#' Floquet theory needs the coefficients to be periodic over the integration
#' window. If they are not, the monodromy matrix is not a monodromy matrix and
#' the exponent is meaningless, so this fails loudly rather than returning a
#' plausible-looking number.
#'
#' @param omega_daily Daily efficacy schedule.
#' @param T_period Period, in days.
#' @param tolerance Fraction of days allowed to differ between the first two
#'   periods.
#' @return `TRUE` invisibly.
#'
#' A handful of differing days is expected and harmless: half-year cycle lengths
#' put the breakpoints on half-days, and rounding them to the daily grid shifts
#' alternate cycle boundaries by one day. Anything beyond that means the window
#' straddles the pre-intervention transient, or covers the wrong number of
#' cycles.
check_schedule_periodicity <- function(omega_daily, T_period, tolerance = 0.01) {

  T_int <- round(T_period)
  if (nrow(omega_daily) < 2 * T_int) {
    stop("check_schedule_periodicity(): the schedule is shorter than two periods.")
  }

  first <- omega_daily[seq_len(T_int), , drop = FALSE]
  second <- omega_daily[T_int + seq_len(T_int), , drop = FALSE]
  n_diff <- sum(first != second)
  frac <- n_diff / length(first)

  if (frac > tolerance) {
    stop(sprintf(paste0("check_schedule_periodicity(): the schedule differs on %d of %d ",
                        "day-patch values (%.1f%%) between the first two periods, so ",
                        "[0, T] is not a period of the intervention and the Floquet ",
                        "exponent would be meaningless."),
                 n_diff, length(first), 100 * frac))
  }
  if (n_diff > 0) {
    message(sprintf("check_schedule_periodicity(): %d day-patch value(s) of %d differ ",
                    n_diff, length(first)),
            "between the first two periods, the expected rounding artefact of a ",
            "half-year cycle length.")
  }

  invisible(TRUE)
}
