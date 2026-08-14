# Effect of the intervention cycle length: Supplementary Figure S5
#
# Same scenario as the illustrative example, deterministic, with the on/off
# period length varied from 1 to 5 years instead of held at 5. One stacked
# panel per duration, showing the AIG growing with the cycle length: a longer
# cycle is also a longer delay before area 2 starts, and a longer break
# afterwards.
#
# Outputs: figures/figureS5_intervention_duration.png

source(here::here("R", "setup.R"))

options(scipen = 7)

# --------------------------------------------------------------------------
# Parameters, identical to the illustrative example except for the cycle length
# --------------------------------------------------------------------------

n_years <- 100
n_days <- 365 * n_years

N <- c(1000, 1000)

p_12 <- 0.10
p_21 <- 0.01

start_year <- 1
start_interv <- start_year * 365
nb_studied_cycles <- 3

omega <- c(0.5, 0.5)
R0 <- c(1.1, 1.1)
r <- 1 / 200

length_intervention_values <- 1:5

x0 <- compute_equilibrium_prevalence(R0, N, r, p_12, p_21)
z0 <- x0 * r * 365 * 1000

params <- list(R0 = R0, p = c(p_12, p_21), omega = omega, rho = rep(1, 2),
               x0 = x0, z0 = z0, r = r)

# --------------------------------------------------------------------------
# One simulation per cycle length
# --------------------------------------------------------------------------

duration_labels <- sprintf("%d year%s", length_intervention_values,
                           ifelse(length_intervention_values > 1, "s", ""))

incidence_list <- vector("list", length(length_intervention_values))
labels_list <- vector("list", length(length_intervention_values))
metrics_list <- vector("list", length(length_intervention_values))

for (k in seq_along(length_intervention_values)) {

  li <- length_intervention_values[k]
  simulations <- simulate_sis(n_days, start_interv, li, N, params)

  incidence_list[[k]] <- bind_rows(
    simulations$synchronous$annual_incidence %>% mutate(scenario = "Synchronous"),
    simulations$asynchronous$annual_incidence %>% mutate(scenario = "Asynchronous")
  ) %>%
    mutate(area = ifelse(variable == "z[1]", "Area1", "Area2"),
           length_intervention = duration_labels[k]) %>%
    select(t, variable, value, year, scenario, area, length_intervention)

  metrics <- compute_metrics(simulations$asynchronous$annual_incidence$value,
                             simulations$synchronous$annual_incidence$value,
                             n_days, start_interv, li, nb_studied_cycles)
  metrics_list[[k]] <- metrics

  labels_list[[k]] <- data.frame(
    AIG_val = metrics[c("AIG_area1", "AIG_area2")],
    AIG_year_val = metrics[c("AIG_year_area1", "AIG_year_area2")],
    area = c("Area1", "Area2"),
    length_intervention = duration_labels[k]
  )
}

all_incidence <- bind_rows(incidence_list) %>%
  mutate(length_intervention = factor(length_intervention, levels = duration_labels))
metrics_labels <- bind_rows(labels_list) %>%
  mutate(length_intervention = factor(length_intervention, levels = duration_labels))

# --------------------------------------------------------------------------
# Figure
# --------------------------------------------------------------------------

visualize_incidence(all_incidence,
                    length_intervention = length_intervention_values,
                    year_start_plot = 0, year_end_plot = 30,
                    max_y = 220, gap = 15,
                    metrics_labels = metrics_labels,
                    facet_rows = "length_intervention", AIG = TRUE)

ggsave(figure_path("figureS5_intervention_duration.png"),
       width = 10, height = 18, dpi = 600)

# --------------------------------------------------------------------------
# Summary across cycle lengths
# --------------------------------------------------------------------------

recap <- data.frame(
  length_intervention = duration_labels,
  AIG = sapply(metrics_list, function(x) x["AIG"]),
  AIGR = sapply(metrics_list, function(x) x["AIGR"]),
  AIG_area1 = sapply(metrics_list, function(x) x["AIG_area1"]),
  AIG_area2 = sapply(metrics_list, function(x) x["AIG_area2"])
)
print(recap)
