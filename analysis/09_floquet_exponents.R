# Floquet exponents of the case studies: Supplementary Table S4
#
# The AIG measures excess cases over a fixed number of cycles. This measures
# whether the epidemic ultimately dies out: the dominant Floquet exponent of the
# disease-free equilibrium under each periodic intervention schedule, negative
# for convergence to elimination, positive for persistence. It is what supports
# the claim, for case study C, that asynchrony does not merely delay elimination
# but prevents it.
#
# The exponent is a property of the linearised system, so no initial prevalence
# and no stochasticity are involved; each run is a single ODE integration over
# one intervention cycle and takes under a second.
#
# Outputs: data/tableS4_floquet_exponents.csv

source(here::here("R", "setup.R"))

N <- c(1000, 1000)

#' Run one scenario and report both exponents
run_floquet <- function(label, length_intervention, p_12, p_21, omega, r, R0) {

  params <- list(R0 = R0, p = c(p_12, p_21), omega = omega, r = r)
  res <- floquet_exponents(length_intervention, N, params)

  data.frame(
    scenario = label,
    length_intervention = length_intervention,
    R0_1 = R0[1], R0_2 = R0[2],
    omega_1 = omega[1], omega_2 = omega[2],
    p_12 = p_12, p_21 = p_21,
    Lambda_sync = res$synchronous$Lambda,
    Lambda_async = res$asynchronous$Lambda
  )
}

# --------------------------------------------------------------------------
# Illustrative example
#
# Not part of Table S4, but a useful check: both schedules reach transmission
# interruption there, so both exponents should be negative.
# --------------------------------------------------------------------------

illustrative <- run_floquet("Illustrative example",
                            length_intervention = 5,
                            p_12 = 0.10, p_21 = 0.01,
                            omega = c(0.5, 0.5), r = 1 / 200,
                            R0 = c(1.1, 1.1))

# --------------------------------------------------------------------------
# Case studies A to D
#
# Parameters come from CASE_STUDIES, the same definitions 08_case_studies.R
# simulates, so the exponents describe exactly the simulated scenarios.
# --------------------------------------------------------------------------

case_results <- bind_rows(lapply(names(CASE_STUDIES), function(label) {
  case <- CASE_STUDIES[[label]]
  resolved <- resolve_case_study(case, N)
  cat("Case study ", label, ": R0 = ",
      paste(round(resolved$R0, 3), collapse = ", "), "\n", sep = "")
  run_floquet(label, case$length_intervention, case$p_12, case$p_21,
              case$omega, case$r, resolved$R0)
}))

results <- bind_rows(illustrative, case_results)

write.csv(results, data_path("tableS4_floquet_exponents.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# Interpretation
# --------------------------------------------------------------------------

print(results %>% select(scenario, Lambda_sync, Lambda_async))

cat("\nSign of the dominant Floquet exponent, by scenario:\n")
verdict <- function(lambda) if (lambda < 0) "elimination" else "persistence"
for (i in seq_len(nrow(results))) {
  row <- results[i, ]
  cat(sprintf("  %-21s synchronous -> %-11s (%.2e), asynchronous -> %-11s (%.2e)%s\n",
              row$scenario,
              verdict(row$Lambda_sync), row$Lambda_sync,
              verdict(row$Lambda_async), row$Lambda_async,
              if (row$Lambda_sync < 0 && row$Lambda_async > 0) {
                "   <- asynchrony prevents elimination"
              } else {
                ""
              }))
}
