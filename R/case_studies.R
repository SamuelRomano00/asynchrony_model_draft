#' Case-study definitions
#'
#' The four settings of Figure 5, kept here rather than in a script because two
#' analyses consume them: the stochastic simulations of `08_case_studies.R` and
#' the Floquet exponents of `09_floquet_exponents.R`. Defining them once
#' guarantees the asymptotic growth rates describe exactly the scenarios that
#' were simulated.
#'
#' Each entry holds the mobility fractions, the schedule, the efficacies and the
#' recovery rate, plus either a target equilibrium prevalence `x0` (from which
#' `R0` is recovered) or an explicit `R0`. Values are those of the Figure 5
#' caption, on the paper's \[0, 0.5\] scale.
CASE_STUDIES <- list(
  A = list(p_12 = 0.2, p_21 = 0.2, start_year = 1, length_intervention = 1.5,
           omega = c(0.4, 0.4), r = 1 / 200, x0 = c(0.3, 0.3), R0 = NULL),
  B = list(p_12 = 0.1, p_21 = 0.01, start_year = 1, length_intervention = 3,
           omega = c(0.5, 0.3), r = 1 / 150, x0 = c(0.2, 0.2), R0 = NULL),
  C = list(p_12 = 0.1, p_21 = 0.01, start_year = 1, length_intervention = 1.5,
           omega = c(0.45, 0.45), r = 1 / 60, x0 = c(0.2, 0.2), R0 = NULL),
  # p_12 and p_21 are the rounded 0.011 and 0.032 of the Figure 5 caption; the
  # exact values simulated are these complements of 0.9/0.91 and 0.9/0.93.
  D = list(p_12 = 1 - 0.9 / 0.91, p_21 = 1 - 0.9 / 0.93, start_year = 1,
           length_intervention = 5, omega = c(0.3, 0.3), r = 1 / 60,
           x0 = NULL, R0 = c(1.05, 1.05))
)

#' Resolve a case study's reproduction numbers and equilibrium prevalence
#'
#' Whichever of `R0` and `x0` the case study specifies, this returns both.
#'
#' @param case One element of [CASE_STUDIES].
#' @param N Population sizes per patch.
#' @return A list with `R0` and `x0`.
#'
#' Case D is specified through `R0 = 1.05`, which yields the 4.76% prevalence
#' quoted in the Figure 5 caption: with equal populations and equal
#' reproduction numbers the equilibrium condition collapses to
#' `R0 = 1 / (1 - X)`, independently of the mobility.
resolve_case_study <- function(case, N) {
  if (is.null(case$R0)) {
    list(R0 = compute_R0(case$x0, N, case$r, case$p_12, case$p_21),
         x0 = case$x0)
  } else {
    list(R0 = case$R0,
         x0 = compute_equilibrium_prevalence(case$R0, N, case$r, case$p_12, case$p_21))
  }
}
