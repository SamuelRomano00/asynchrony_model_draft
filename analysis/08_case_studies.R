# Case studies A to D: Figure 5, Supplementary Figure S12 and Table S4
#
# Four settings chosen among the high-AIG scenarios of the sensitivity analysis,
# with intervention timings and recovery times consistent with real malaria
# parasites and control programmes:
#   A  high prevalence, P. falciparum-like recovery, 1.5-year ITN cycles
#   B  the same but with 3-year cycles, i.e. campaigns spaced by resource limits
#   C  high prevalence, short infectiousness (P. vivax or strong case
#      management); synchrony interrupts transmission, asynchrony does not
#   D  low prevalence, short infectiousness, 5-year elimination campaigns
#
# Each is simulated over 1000 replicates under both schedules.
#
# Outputs:
#   figures/figure5_case_studies.png
#   figures/figureS12_case_studies_interruption.png
#   data/tableS4_case_studies.csv

source(here::here("R", "setup.R"))

options(scipen = 7)

n_years <- 100
n_days <- 365 * n_years
n_reps <- 1000

simulator <- build_SIS_stochastic_simulator()

N <- c(1000, 1000)
nb_studied_cycles <- 3

# Case-study parameters live in R/case_studies.R, because 09_floquet_exponents.R
# reads the same definitions.

#' Simulate one case study and build its figures and table row
#'
#' @param case A list from [CASE_STUDIES].
#' @param label The case-study letter, used in the table.
#' @return A list with `incidence`, `metrics`, `plot_incidence`,
#'   `plot_interruption` and `table`.
run_case <- function(case, label) {

  start_interv <- case$start_year * 365

  resolved <- resolve_case_study(case, N)
  R0 <- resolved$R0
  x0 <- resolved$x0
  cat("Case study ", label, ": R0 = ", paste(round(R0, 3), collapse = ", "),
      "\n", sep = "")

  params <- list(R0 = R0, p = c(case$p_12, case$p_21), omega = case$omega,
                 x0 = x0, r = case$r)

  simulations <- simulate_sis_stochastic(simulator, n_days, start_interv,
                                         case$length_intervention, N, params,
                                         n_reps = n_reps, seed = 1)

  incidence <- bind_rows(
    simulations$synchronous$annual_incidence %>% mutate(scenario = "Synchronous"),
    simulations$asynchronous$annual_incidence %>% mutate(scenario = "Asynchronous")
  ) %>%
    mutate(area = ifelse(variable == "z[1]", "Area1", "Area2"))

  metrics <- compute_metrics_stochastic(simulations$asynchronous, simulations$synchronous,
                                        n_days, start_interv, case$length_intervention,
                                        nb_studied_cycles)

  metrics_labels <- data.frame(
    area = c("Area1", "Area2"),
    AIG_val = c(metrics$summary$AIG_area1_mean, metrics$summary$AIG_area2_mean),
    AIG_year_val = c(metrics$summary$AIG_year_area1_mean, metrics$summary$AIG_year_area2_mean)
  )

  list(
    incidence = incidence,
    metrics = metrics,
    plot_incidence = visualize_incidence_stochastic(
      incidence, case$length_intervention,
      year_start_plot = 0, year_end_plot = 25, gap = 20,
      metrics_labels = metrics_labels, start_year = case$start_year, AIG = TRUE
    ),
    plot_interruption = visualize_interruption_stochastic(
      incidence, year_start_plot = 0, year_end_plot = 25, smooth = TRUE
    ),
    table = build_publication_table(metrics$summary, incidence, start_interv) %>%
      mutate(`Case study` = label, .before = 1)
  )
}

results <- Map(run_case, CASE_STUDIES, names(CASE_STUDIES))

# --------------------------------------------------------------------------
# Stacked figures
#
# Only the top panel keeps its facet strips and only the bottom one its x-axis
# title, so the stack reads as a single figure.
# --------------------------------------------------------------------------

stack_panels <- function(plots, x_title, y_title = NULL) {
  n <- length(plots)
  plots <- lapply(seq_along(plots), function(i) {
    p <- plots[[i]] + labs(x = if (i == n) x_title else NULL)
    if (!is.null(y_title)) p <- p + labs(y = y_title)
    if (i > 1) p <- p + theme(strip.text.x = element_blank())
    p
  })

  patchwork::wrap_plots(plots, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A", tag_suffix = ".") &
    theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))
}

stack_panels(lapply(results, `[[`, "plot_incidence"), "Years", "Annual incidence")
ggsave(figure_path("figure5_case_studies.png"),
       width = 10, height = 14, units = "in", dpi = 600)

stack_panels(lapply(results, `[[`, "plot_interruption"), "Years")
ggsave(figure_path("figureS12_case_studies_interruption.png"),
       width = 10, height = 17, units = "in", dpi = 600)

# --------------------------------------------------------------------------
# Combined table
# --------------------------------------------------------------------------

write.csv(bind_rows(lapply(results, `[[`, "table")),
          data_path("tableS4_case_studies.csv"), row.names = FALSE)
