# Parameter sets where asynchrony REDUCES incidence in area 1: Figure S11
#
# Each set is re-simulated and drawn with the same conventions as the case-study
# figures, sorted by increasing AIG per year, so the most negative case comes
# first.
#
# The simulation database only keeps aggregated metrics, not the incidence time
# series, so the trajectories have to be replayed. Two ODE solves per parameter
# set: a few seconds, no parallelism needed.
#
# Outputs:
#   figures/figureS11_negative_aig_cases.png       the assembled grid, two
#                                                  parameter sets per row
#   figures/figureS11_negative_aig_cases/case_*.png one panel per parameter set

source(here::here("R", "setup.R"))

# --------------------------------------------------------------------------
# Parameter sets
# --------------------------------------------------------------------------

# Parameter sets with no endemic equilibrium have AIG_area1 == 0 and so would
# never be selected here anyway, but the same loader is used as in 05 and 06 so
# that every count reported from the database refers to the same sample. See
# load_simulation_database() in R/batch.R.
params <- load_simulation_database() %>%
  filter(AIG_area1 < 0) %>%
  arrange(AIG_year_area1)

cat("Parameter sets with a negative AIG:", nrow(params), "\n")

# --------------------------------------------------------------------------
# Re-simulation
# --------------------------------------------------------------------------

# x0, z0 and r are taken from the file rather than recomputed, so the
# trajectories reproduced are exactly those that produced these negative AIGs
# and not merely nearby ones.
simulate_one_case <- function(row) {
  
  n_days <- 30000
  N <- c(1000, 1000)
  start_interv <- 365
  nb_studied_cycles <- 3
  length_intervention <- row$time_intervention
  
  params_row <- list(
    R0 = c(row$R0_1, row$R0_2),
    p = c(row$p_12, row$p_21),
    omega = c(row$omega_1, row$omega_2),
    rho = rep(1, 2),
    x0 = c(row$x0_1, row$x0_2),
    z0 = c(row$z0_1, row$z0_2),
    r = row$r
  )
  
  simulations <- simulate_sis(n_days, start_interv, length_intervention, N, params_row)
  
  metrics <- compute_metrics(simulations$asynchronous$annual_incidence$value,
                             simulations$synchronous$annual_incidence$value,
                             n_days, start_interv, length_intervention,
                             nb_studied_cycles)
  
  to_long <- function(simulation, scenario_name) {
    simulation$annual_incidence %>%
      mutate(area = ifelse(variable == "z[1]", "Area1", "Area2"),
             scenario = scenario_name) %>%
      select(t, value, area, scenario)
  }
  
  list(incidence = bind_rows(to_long(simulations$synchronous, "Synchronous"),
                             to_long(simulations$asynchronous, "Asynchronous")),
       metrics = metrics,
       length_intervention = length_intervention)
}

cases <- lapply(seq_len(nrow(params)), function(i) simulate_one_case(params[i, ]))

# Check that replaying the parameter sets reproduces the AIG per year recorded
# in the file.
check <- data.frame(
  case = seq_len(nrow(params)),
  AIG_year_area1_file = params$AIG_year_area1,
  AIG_year_area1_resim = sapply(cases, function(cs) unname(round(cs$metrics["AIG_year_area1"])))
)
check$match <- check$AIG_year_area1_file == check$AIG_year_area1_resim
print(check)
if (!all(check$match)) {
  warning("Some re-simulated AIG per year values differ from the file; see the table above.")
}

# --------------------------------------------------------------------------
# Panels
# --------------------------------------------------------------------------

# Heights are derived per area from the incidence actually reached in the
# plotting window, because these ten cases span very different scales.
build_panel <- function(case) {
  
  li <- case$length_intervention
  year_start_plot <- 0
  year_end_plot <- max(12, ceiling(1 + 6 * li + 2))
  
  max_by_area <- case$incidence %>%
    filter(t / 365 >= year_start_plot, t / 365 <= year_end_plot) %>%
    group_by(area) %>%
    summarise(m = max(value, na.rm = TRUE), .groups = "drop")
  max_vec <- setNames(max_by_area$m, max_by_area$area)
  
  metrics_labels <- data.frame(
    area = c("Area1", "Area2"),
    value = c(paste0("AIG per year: ", round(case$metrics["AIG_year_area1"])),
              paste0("AIG per year: ", round(case$metrics["AIG_year_area2"])))
  )
  
  plot_incidence_panel(incidence_i = case$incidence,
                       metrics_labels_i = metrics_labels,
                       length_intervention_i = li,
                       year_start_plot = year_start_plot,
                       year_end_plot = year_end_plot,
                       max_y = max_vec * 1.3,
                       gap = max_vec * 0.09,
                       show_aig_area = TRUE,
                       text_scale = 1.6,
                       label_y = max_vec * 1.2)
}

panels <- lapply(cases, build_panel)

# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

panels_no_legend <- lapply(panels, function(p) p + theme(legend.position = "none"))
shared_legend <- cowplot::get_legend(panels[[1]] + theme(legend.position = "bottom"))

n_cases <- length(panels_no_legend)
n_rows <- ceiling(n_cases / 2)

grid_with_lines <- cowplot::ggdraw(
  cowplot::plot_grid(plotlist = panels_no_legend, ncol = 2, nrow = n_rows)
) +
  cowplot::draw_line(x = c(0.5, 0.5), y = c(0, 1), colour = "black", linewidth = 0.8)

for (y_boundary in seq_len(n_rows - 1) / n_rows) {
  grid_with_lines <- grid_with_lines +
    cowplot::draw_line(x = c(0, 1), y = c(y_boundary, y_boundary),
                       colour = "black", linewidth = 0.8)
}

figure_neg <- cowplot::plot_grid(grid_with_lines, shared_legend,
                                 ncol = 1, rel_heights = c(0.94, 0.06))

ggsave(figure_path("figureS11_negative_aig_cases.png"), figure_neg,
       width = 24, height = 6.4 * n_rows, limitsize = FALSE)

case_dir <- figure_path("figureS11_negative_aig_cases")
dir.create(case_dir, showWarnings = FALSE)
for (i in seq_along(panels)) {
  ggsave(file.path(case_dir, paste0("case_", i, ".png")), panels[[i]],
         width = 10, height = 5)
}