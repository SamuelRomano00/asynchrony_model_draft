# Testing the human-to-human assumption: Supplementary Figure S2
#
# Reproduces the illustrative example with an explicit vector-host model, to
# check that the asynchrony effect is not an artefact of collapsing mosquito
# dynamics into a human-to-human transmission term. Both models share the same
# populations, mobility, R0 and intervention schedule; the Ross-Macdonald
# mosquito-to-human ratio is derived from the target R0. Deterministic, unlike
# the stochastic version used for Figure 2.
#
# Outputs: figures/figureS2_sis_vs_ross_macdonald.png

source(here::here("R", "setup.R"))

options(scipen = 7)

# --------------------------------------------------------------------------
# Shared parameters
# --------------------------------------------------------------------------

n_years <- 100
n_days <- 365 * n_years

N <- c(1000, 1000)

p_12 <- 0.10
p_21 <- 0.01

start_year <- 1
start_interv <- start_year * 365
length_intervention <- 5
nb_studied_cycles <- 3

omega <- c(0.5, 0.5)
R0 <- c(1.1, 1.1)
r <- 1 / 200
rho <- 1

# --------------------------------------------------------------------------
# SIS model
# --------------------------------------------------------------------------

x0_sis <- compute_equilibrium_prevalence(R0, N, r, p_12, p_21)
z0_sis <- x0_sis * r * 365 * 1000

params_sis <- list(R0 = R0, p = c(p_12, p_21), omega = omega, rho = rep(rho, 2),
                   x0 = x0_sis, z0 = z0_sis, r = r)

sim_sis <- simulate_sis(n_days, start_interv, length_intervention, N, params_sis)

# --------------------------------------------------------------------------
# Ross-Macdonald model
#
# Entomological values are those given in the caption of Supplementary
# Figure S2. The human recovery rate and R0 are shared with the SIS run, and m
# is derived from R0, so the two models start from the same transmission
# intensity.
# --------------------------------------------------------------------------

a_rm <- 0.15        # biting rate
b_rm <- 0.10        # transmission probability, mosquito to human
c_rm <- 0.107       # transmission probability, human to mosquito
mu_rm <- 0.15       # mosquito mortality
tau_rm <- 10        # extrinsic incubation period, in days

eq_rm <- compute_equilibrium_prevalence_RM(N, a_rm, b_rm, c_rm, r, mu_rm, tau_rm,
                                           R0, p_12, p_21)
x0_rm <- eq_rm$X    # human equilibrium prevalence
v0_rm <- eq_rm$Y    # mosquito equilibrium prevalence
z0_rm <- rho * r * x0_rm * 365 * 1000

cat("Mosquito-to-human ratio implied by R0:",
    R0_to_m(R0, a_rm, b_rm, c_rm, r, mu_rm, tau_rm), "\n")
cat("R0 recovered from the next-generation matrix:",
    compute_R0_RM(N, a_rm, b_rm, c_rm, r, mu_rm, tau_rm, R0, p_12, p_21), "\n")

params_rm <- list(R0 = R0, p = c(p_12, p_21), omega = omega,
                  x0 = x0_rm, z0 = z0_rm, v0 = v0_rm,
                  a = a_rm, b = b_rm, c_rate = c_rm,
                  r = r, mu = mu_rm, tau = tau_rm, rho = rho)

sim_rm <- simulate_rm(n_days, start_interv, length_intervention, N, params_rm)

# --------------------------------------------------------------------------
# Combined figure
# --------------------------------------------------------------------------

stack_scenarios <- function(simulations, model_name) {
  bind_rows(
    simulations$synchronous$annual_incidence %>% mutate(scenario = "Synchronous"),
    simulations$asynchronous$annual_incidence %>% mutate(scenario = "Asynchronous")
  ) %>%
    mutate(area = ifelse(variable == "z[1]", "Area1", "Area2"), model = model_name) %>%
    select(t, variable, value, year, scenario, area, model)
}

all_incidence <- bind_rows(stack_scenarios(sim_sis, "SIS"),
                           stack_scenarios(sim_rm, "Ross-Macdonald")) %>%
  mutate(model = factor(model, levels = c("SIS", "Ross-Macdonald")))

metrics_sis <- compute_metrics(sim_sis$asynchronous$annual_incidence$value,
                               sim_sis$synchronous$annual_incidence$value,
                               n_days, start_interv, length_intervention,
                               nb_studied_cycles)
metrics_rm <- compute_metrics(sim_rm$asynchronous$annual_incidence$value,
                              sim_rm$synchronous$annual_incidence$value,
                              n_days, start_interv, length_intervention,
                              nb_studied_cycles)

metrics_labels <- bind_rows(
  data.frame(AIG_val = metrics_sis[c("AIG_area1", "AIG_area2")],
             AIG_year_val = metrics_sis[c("AIG_year_area1", "AIG_year_area2")],
             area = c("Area1", "Area2"), model = "SIS"),
  data.frame(AIG_val = metrics_rm[c("AIG_area1", "AIG_area2")],
             AIG_year_val = metrics_rm[c("AIG_year_area1", "AIG_year_area2")],
             area = c("Area1", "Area2"), model = "Ross-Macdonald")
) %>%
  mutate(model = factor(model, levels = c("SIS", "Ross-Macdonald")))

visualize_incidence(all_incidence, length_intervention,
                    year_start_plot = 0, year_end_plot = 30,
                    max_y = 185, gap = 6,
                    metrics_labels = metrics_labels,
                    facet_rows = "model", AIG = TRUE)

ggsave(figure_path("figureS2_sis_vs_ross_macdonald.png"),
       width = 10, height = 10, dpi = 300)

# --------------------------------------------------------------------------
# Elimination times and cases accrued
#
# The numbers quoted in the Supplementary Figure S2 text: how much later
# elimination is reached under asynchrony, and how many cases accrue before it.
# --------------------------------------------------------------------------

#' First year in which annual incidence in an area falls below a threshold
elimination_year <- function(simulation, series, threshold = 1e-5) {
  hit <- simulation$annual_incidence %>%
    filter(.data$variable == series, .data$value < threshold)
  if (nrow(hit) == 0) return(Inf)
  min(hit$t) / 365
}

#' Cumulated cases across both areas before a given year
cases_before <- function(simulation, year) {
  sum(simulation$annual_incidence$value[simulation$annual_incidence$t / 365 < year])
}

report_model <- function(simulations, x0, metrics, model_name) {
  elim_sync <- max(elimination_year(simulations$synchronous, "z[1]"),
                   elimination_year(simulations$synchronous, "z[2]"))
  elim_async <- max(elimination_year(simulations$asynchronous, "z[1]"),
                    elimination_year(simulations$asynchronous, "z[2]"))

  cat("\n---", model_name, "---\n")
  cat("Equilibrium prevalence        :", x0, "\n")
  cat("Elimination year, synchronous :", elim_sync, "\n")
  cat("Elimination year, asynchronous:", elim_async, "\n")
  cat("Delay caused by asynchrony    :", elim_async - elim_sync, "years\n")
  cat("Cases before elimination, synchronous :",
      cases_before(simulations$synchronous, elim_sync), "\n")
  cat("Cases before elimination, asynchronous:",
      cases_before(simulations$asynchronous, elim_async), "\n")
  print(metrics)
}

report_model(sim_sis, x0_sis, metrics_sis, "SIS")
report_model(sim_rm, x0_rm, metrics_rm, "Ross-Macdonald")
