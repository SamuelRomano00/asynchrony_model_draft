# Illustrative example: Figure 2, Table 2 and Supplementary Figure S1
#
# The paper's illustrative example, simulated as a continuous-time Markov chain
# over 1000 replicates. Two connected areas with the same R0 and the same
# intervention, run once with both areas' vector control switched on and off
# together, and once with area 2 phase-shifted by one period. Area 1 receives
# an identical schedule in both scenarios, so any difference there comes purely
# from the timing in area 2.
#
# Outputs:
#   figures/figure2_illustrative_example.png
#   figures/figureS1_interruption.png
#   data/table2_illustrative_example.csv

source(here::here("R", "setup.R"))

options(scipen = 7)

n_years <- 100
n_days <- 365 * n_years
n_reps <- 1000

# Compiling the TiPS simulator takes a few seconds; build it once and reuse it.
simulator <- build_SIS_stochastic_simulator()

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

N <- c(1000, 1000)

p_12 <- 0.10   # area-1 residents spend 10% of their time in area 2
p_21 <- 0.01   # area 2 is nearly closed

start_year <- 1
start_interv <- start_year * 365
length_intervention <- 5
nb_studied_cycles <- 3

omega <- c(0.5, 0.5)       # vector-control efficacy while the intervention is on
R0 <- c(1.1, 1.1)          # both areas just above the epidemic threshold
r <- 1 / 200               # mean infectious period of 200 days

x0 <- compute_equilibrium_prevalence(R0, N, r, p_12, p_21)

params <- list(R0 = R0, p = c(p_12, p_21), omega = omega, x0 = x0, r = r)

# --------------------------------------------------------------------------
# Simulation
# --------------------------------------------------------------------------

simulations <- simulate_sis_stochastic(simulator, n_days, start_interv,
                                       length_intervention, N, params,
                                       n_reps = n_reps, seed = 1)

all_incidence <- bind_rows(
  simulations$synchronous$annual_incidence %>% mutate(scenario = "Synchronous"),
  simulations$asynchronous$annual_incidence %>% mutate(scenario = "Asynchronous")
) %>%
  mutate(area = ifelse(variable == "z[1]", "Area1", "Area2"))

metrics <- compute_metrics_stochastic(simulations$asynchronous, simulations$synchronous,
                                      n_days, start_interv, length_intervention,
                                      nb_studied_cycles)

# --------------------------------------------------------------------------
# Figures and table
# --------------------------------------------------------------------------

metrics_labels <- data.frame(
  area = c("Area1", "Area2"),
  AIG_val = c(metrics$summary$AIG_area1_mean, metrics$summary$AIG_area2_mean),
  AIG_year_val = c(metrics$summary$AIG_year_area1_mean, metrics$summary$AIG_year_area2_mean)
)

visualize_incidence_stochastic(all_incidence, length_intervention,
                               year_start_plot = 0, year_end_plot = 30, gap = 5,
                               metrics_labels = metrics_labels,
                               show_aig_area = TRUE, AIG = TRUE)
ggsave(figure_path("figure2_illustrative_example.png"),
       width = 10, height = 7, units = "in", dpi = 600)

visualize_interruption_stochastic(all_incidence, year_start_plot = 0,
                                  year_end_plot = 30, smooth = TRUE)
ggsave(figure_path("figureS1_interruption.png"),
       width = 9, height = 5, units = "in", dpi = 600)

write.csv(build_publication_table(metrics$summary, all_incidence, start_interv),
          data_path("table2_illustrative_example.csv"), row.names = FALSE)
