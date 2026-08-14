#' Stochastic version of the connected SIS model
#'
#' A continuous-time Markov chain with the same structure as the deterministic
#' SIS model, simulated exactly with TiPS. Used to quantify how often
#' asynchrony delays or prevents transmission interruption, which a
#' deterministic model cannot express.

#' Build the TiPS simulator
#'
#' Compiles the reaction set to C++. Compilation takes roughly ten seconds, so
#' build the simulator once per session and reuse the returned function across
#' parameter sets and replicates.
#'
#' Multiplying the deterministic `deriv(x[i])` equation by `N[i]` and
#' substituting `K[j]`, the population-weighted prevalence met in patch j,
#' gives two elementary reactions per patch:
#'   `S[i] -> I[i]` at rate `r * sum_j( p[i,j] * (1 - omega_j) * R0_j * K[j] ) * S[i]`
#'   `I[i] -> S[i]` at rate `r * I[i]`
#' With `p11 = a`, `p12 = b`, `p21 = c`, `p22 = d`, that makes
#'   `K1 = (a*I1 + c*I2) / (a*N1 + c*N2)`
#'   `K2 = (b*I1 + d*I2) / (b*N1 + d*N2)`
#' which are inlined into each rate expression rather than passed through TiPS's
#' `functions` argument, so every symbol is picked up by its parameter parser.
#'
#' @return A simulator function, to be passed to [run_SIS_stochastic()].
build_SIS_stochastic_simulator <- function() {
  reactions <- c(
    "S1 [r*(a*(1-om1)*R01*(a*I1+c*I2)/(a*N1+c*N2) + b*(1-om2)*R02*(b*I1+d*I2)/(b*N1+d*N2))*S1] -> I1",
    "I1 [r*I1] -> S1",
    "S2 [r*(c*(1-om1)*R01*(a*I1+c*I2)/(a*N1+c*N2) + d*(1-om2)*R02*(b*I1+d*I2)/(b*N1+d*N2))*S2] -> I2",
    "I2 [r*I2] -> S2"
  )
  build_simulator(reactions)
}

#' Run a TiPS trajectory, retrying only on an unusable return value
#'
#' @param f Simulator function.
#' @param ... Passed to `f`.
#' @param seed Base seed. Retries perturb it, because re-calling a stochastic
#'   simulator with an identical seed reproduces the same failure forever.
#' @param max_tries Maximum number of attempts before giving up.
#' @return Whatever `f` returns on its first usable attempt.
#'
#' Extinction is not a failure and is not retried: TiPS returns a non-empty,
#' truncated trajectory ending at the extinction time with `I1 = I2 = 0`, which
#' the downstream summaries handle correctly. Only a genuinely empty return
#' value triggers a retry, so there is no elimination-outcome selection bias.
safe_run <- function(f, ..., seed = 0, max_tries = 50) {
  for (k in seq_len(max_tries)) {
    out <- f(..., seed = seed + 1000L * (k - 1L))
    if (length(out)) {
      if (k > 1) {
        warning("safe_run(): simulator returned an empty object ", k - 1,
                " time(s) before succeeding (seed ", seed, ").")
      }
      return(out)
    }
  }
  stop("safe_run(): simulator returned nothing usable after ", max_tries,
       " attempts (seed ", seed, ").")
}

#' Run the stochastic SIS model
#'
#' @param simulator Simulator from [build_SIS_stochastic_simulator()].
#' @param parameters Named list with
#'   \describe{
#'     \item{N}{`c(N1, N2)` population sizes.}
#'     \item{R0}{`c(R0_1, R0_2)`.}
#'     \item{p}{2x2 mobility matrix, see [mobility_matrix()].}
#'     \item{r}{Recovery rate.}
#'     \item{omega}{`list(omega1 = , omega2 = )`, each a vector with one
#'       efficacy value per period, i.e. `length(times) - 1` values.}
#'     \item{x0}{`c(x0_1, x0_2)` initial prevalence, converted to counts.}
#'   }
#' @param times Strictly increasing period boundaries of the intervention
#'   schedule, directly reusable as TiPS's `times` argument.
#' @param n_reps Number of independent replicates.
#' @param method `"exact"` (Gillespie, default), `"approximate"` (tau-leap) or
#'   `"mixed"`.
#' @param tau Tau-leap step, used only by the approximate methods.
#' @param seed Base seed; replicate k uses `seed + k - 1`.
#' @param verbose Passed through to TiPS.
#' @return A list with `simulation` (the raw per-replicate TiPS objects),
#'   `prevalence` and `annual_incidence`, the latter two in the same long
#'   format as [run_SIS_deterministic()] plus a `rep` column.
run_SIS_stochastic <- function(simulator, parameters, times, n_reps = 1,
                               method = "exact", tau = 0.01, seed = NULL,
                               verbose = FALSE) {

  N1 <- unname(parameters$N[1]); N2 <- unname(parameters$N[2])
  # Local names avoid shadowing `c`; the TiPS parameter names below must stay
  # a/b/c/d because they appear verbatim in the reaction expressions.
  m11 <- parameters$p[1, 1]; m12 <- parameters$p[1, 2]
  m21 <- parameters$p[2, 1]; m22 <- parameters$p[2, 2]

  paramValues <- list(
    r = unname(parameters$r),
    R01 = unname(parameters$R0[1]), R02 = unname(parameters$R0[2]),
    om1 = unname(parameters$omega$omega1), om2 = unname(parameters$omega$omega2),
    a = m11, b = m12, c = m21, d = m22,
    N1 = N1, N2 = N2
  )

  # TiPS indexes the per-period parameter vectors without a bounds check: one
  # value short is a segfault that kills the session, not a recoverable error.
  n_periods <- length(times) - 1
  for (nm in c("om1", "om2")) {
    if (length(paramValues[[nm]]) < n_periods) {
      stop(sprintf(paste0("run_SIS_stochastic(): %s has %d value(s) but `times` defines ",
                          "%d period(s). TiPS does not recycle short parameter vectors. ",
                          "Supply one value per period."),
                   nm, length(paramValues[[nm]]), n_periods))
    }
  }
  if (is.unsorted(times, strictly = TRUE)) {
    stop("run_SIS_stochastic(): `times` must be strictly increasing; got ",
         paste(times, collapse = ", "))
  }

  I1_0 <- unname(round(parameters$x0[1] * N1))
  I2_0 <- unname(round(parameters$x0[2] * N2))
  initialStates <- c(S1 = N1 - I1_0, I1 = I1_0, S2 = N2 - I2_0, I2 = I2_0)

  n_days <- max(times)
  # The yearly grid must match the deterministic one exactly, because
  # compute_metrics() locates each patch's block by integer position: one entry
  # per complete year plus a t = 0 placeholder, a trailing partial year dropped.
  n_years_full <- floor(n_days / 365)
  year_grid <- seq(0, n_years_full) * 365

  reps <- vector("list", n_reps)
  prevalence_list <- vector("list", n_reps)
  incidence_list <- vector("list", n_reps)
  grid <- 0:n_days

  for (k in seq_len(n_reps)) {
    this_seed <- if (is.null(seed)) 0 else seed + k - 1
    sim <- safe_run(simulator,
                    paramValues = paramValues,
                    initialStates = initialStates,
                    times = times,
                    method = method,
                    tau = tau,
                    seed = this_seed,
                    verbose = verbose)
    reps[[k]] <- sim
    df <- sim$traj

    ## Daily prevalence, last value carried forward.
    # A trajectory truncated by extinction is handled by the same mechanism:
    # the final state (I = 0) is carried to the end of the grid.
    idx <- findInterval(grid, df$Time)
    idx[idx == 0] <- 1
    st <- df[idx, c("S1", "I1", "S2", "I2")]
    prevalence_list[[k]] <- data.frame(
      t = rep(grid, 2),
      variable = rep(c("x[1]", "x[2]"), each = length(grid)),
      value = c(st$I1 / N1, st$I2 / N2),
      rep = k
    )

    ## Annual incidence: count S -> I events per patch per year.
    count_incidence <- function(infection_suffix) {
      ev <- df[grepl(infection_suffix, df$Reaction, fixed = TRUE), ]
      yr <- ceiling(ev$Time / 365)
      keep <- yr >= 1 & yr <= n_years_full
      ev <- ev[keep, ]; yr <- yr[keep]
      counts <- as.numeric(tapply(ev$Nrep, factor(yr, levels = seq_len(n_years_full)),
                                  sum, default = 0))
      c(counts[1], counts)   # t = 0 placeholder, mirroring the deterministic grid
    }
    incidence_list[[k]] <- rbind(
      data.frame(t = year_grid, variable = "z[1]", value = count_incidence("-> I1"), rep = k),
      data.frame(t = year_grid, variable = "z[2]", value = count_incidence("-> I2"), rep = k)
    )
  }

  list(
    simulation = reps,
    prevalence = do.call(rbind, prevalence_list),
    annual_incidence = do.call(rbind, incidence_list)
  )
}
