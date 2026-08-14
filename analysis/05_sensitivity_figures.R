# Uncertainty and sensitivity analysis: Figure 3 and Supplementary S3, S4, S6, T1
#
# Reads the CSVs written by 04_sensitivity_simulations.R. No simulation happens
# here, so this script is cheap to re-run while adjusting a figure.
#
# Outputs:
#   figures/figure3_sensitivity_AIG.png    distribution, Sobol indices and
#                                          parameter profiles, AIG per year in area 1
#   figures/figureS3_sensitivity_AIGR.png  the same for the AIGR per year
#   figures/figureS4_near_elimination.png  AIG and AIGR restricted to
#                                          R_C,1 < 1 and R_C,1 in [0.9, 1.1]
#   figures/figureS6_recovery_duration_heatmap.png  mean AIG by recovery time
#                                          and intervention duration
#   data/tableS1_metric_quantiles.csv       quantiles of both metrics

source(here::here("R", "setup.R"))

# Parameter sets with no endemic equilibrium are excluded here rather than
# figure by figure: their AIG is a finite zero while their AIGR is undefined, so
# reading the CSV directly would put the AIG and AIGR panels on different
# samples. See load_simulation_database() in R/batch.R.
df <- load_simulation_database()

# The Sobol indices are computed on a design that was restricted the same way in
# 04_sensitivity_simulations.R, on the stricter condition that every one of the
# 2k+2 blocks of a draw be viable. Both sides of this script therefore describe
# viable parameter sets only.
sobol_AIG <- read.csv(data_path("sobol_AIG.csv"))
sobol_AIGR <- read.csv(data_path("sobol_AIGR.csv"))

# A rank-based version of the AIGR indices is also available as
# sobol_AIGR_rank.csv; substitute it here to check the raw-scale indices against
# a version insensitive to the near-elimination tail.

# --------------------------------------------------------------------------
# Reusable panel builders
# --------------------------------------------------------------------------

#' Histogram of a metric, trimmed to its 1st-99th percentile range
#'
#' @param data Simulation database.
#' @param metric Column name, as a string.
#' @param title,subtitle Panel annotations.
#' @param x_label x-axis label.
#' @param binwidth Fixed bin width, or `NULL` to use a thirtieth of the trimmed
#'   range.
#' @param percent Whether to format the x axis as a percentage.
#' @return A ggplot object.
#'
#' Non-finite values are dropped before any range, bin or quantile computation:
#' a single Inf, which AIGR produces whenever a synchronous run has already
#' eliminated transmission, would otherwise corrupt the binning.
make_metric_hist <- function(data, metric, title, subtitle = NULL, x_label,
                             binwidth = NULL, percent = FALSE) {
  
  values <- data[[metric]]
  first_centile <- quantile(values, 0.01, na.rm = TRUE)
  last_centile <- quantile(values, 0.99, na.rm = TRUE)
  quantile75 <- quantile(values, 0.75, na.rm = TRUE)
  quantile90 <- quantile(values, 0.90, na.rm = TRUE)
  
  data_filtered <- data %>%
    filter(is.finite(.data[[metric]]),
           .data[[metric]] >= first_centile,
           .data[[metric]] <= last_centile)
  filtered_values <- data_filtered[[metric]]
  
  # AIG has a natural fixed unit (cases per year per 1000), so a fixed bin width
  # is meaningful; AIGR is a dimensionless ratio, so its bins are derived from
  # the observed range instead.
  if (is.null(binwidth)) {
    binwidth <- diff(range(filtered_values)) / 30
  }
  breaks_seq <- seq(min(filtered_values), max(filtered_values) + binwidth, by = binwidth)
  peak <- max(graphics::hist(filtered_values, breaks = breaks_seq, plot = FALSE)$counts)
  
  p <- ggplot(data_filtered, aes(x = .data[[metric]])) +
    geom_histogram(binwidth = binwidth, boundary = 0, fill = "skyblue", color = "white") +
    labs(x = x_label, y = "Frequency", title = title, subtitle = subtitle) +
    geom_vline(aes(xintercept = quantile75)) +
    annotate("text", x = quantile75, y = peak * 0.9, label = "q75", hjust = -0.15) +
    geom_vline(aes(xintercept = quantile90)) +
    annotate("text", x = quantile90, y = peak * 0.78, label = "q90", hjust = -0.15) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          plot.title = element_text(face = "bold", size = 10),
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(size = 7.5),
          axis.text.y = element_text(size = 7.5))
  
  if (percent) {
    p <- p + scale_x_continuous(labels = scales::percent)
  } else {
    p <- p + scale_x_continuous(breaks = seq(0, 200, by = 10))
  }
  
  p
}

#' First-order and total Sobol indices, one row per parameter
#'
#' @param sobol_df Data frame read from one of the `sobol_*.csv` files. The
#'   confidence-interval columns are named `min..c.i.`/`max..c.i.` because
#'   `read.csv()` sanitises the spaces and dots of the original headers.
#' @return A ggplot object.
plot_sobol_indices <- function(sobol_df) {
  sobol_df %>%
    mutate(param = factor(X, levels = PARAM_ORDER)) %>%
    ggplot() +
    geom_point(aes(x = original, y = param, col = index, shape = index), size = 2) +
    geom_linerange(aes(xmin = min..c.i., xmax = max..c.i., y = param, col = index),
                   linewidth = 1) +
    labs(y = "", x = "Sobol indices", shape = "", col = "", title = "B.") +
    scale_color_manual(values = c("cyan3", "darkblue")) +
    scale_y_discrete(labels = PARAM_LABELS) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 10),
          axis.title.x = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.83, 1),
          legend.background = element_rect(colour = 1),
          legend.box = "horizontal",
          strip.background = element_rect(fill = "white", color = "white"))
}

#' Ten parameter profiles on a common y scale
#'
#' @param metric Column name of the metric, as a string.
#' @param y_label Shared y-axis label.
#' @param y_percent Whether to format the y axis as a percentage.
#' @return A patchwork object.
build_profile_panels <- function(metric, y_label, y_percent = FALSE) {
  
  specs <- list(
    list(param = "R0_1", bin = 0.1, label = "R0 in area 1", y_name = FALSE),
    list(param = "R0_2", bin = 0.1, label = "R0 in area 2", y_name = FALSE),
    list(param = "p_12", bin = 0.05, label = expression(p["1,2"]), y_name = FALSE),
    list(param = "p_21", bin = 0.05, label = expression(p["2,1"]), y_name = FALSE),
    list(param = "RC_1", bin = 0.05, label = expression(R[C] ~ " in area 1"), y_name = TRUE),
    list(param = "RC_2", bin = 0.05, label = expression(R[C] ~ " in area 2"), y_name = FALSE),
    list(param = "omega_1", bin = 0.05, label = "Intervention efficiency in area 1", y_name = FALSE),
    list(param = "omega_2", bin = 0.05, label = "Intervention efficiency in area 2", y_name = FALSE),
    list(param = "rinv", bin = 10, label = "Duration of recovery", y_name = FALSE),
    list(param = "time_intervention", bin = 0.5, label = "Duration of the intervention", y_name = FALSE)
  )
  
  panels <- lapply(specs, function(s) {
    plot_metric_by_parameter(df, s$param, metric, s$bin, s$label,
                             y_name = s$y_name, y_label = y_label,
                             y_percent = y_percent)
  })
  
  # A common y scale across all ten panels, so their slopes are comparable.
  y_ranges <- sapply(panels, function(p) ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  ylim_common <- c(min(y_ranges[1, ]), max(y_ranges[2, ]))
  
  ((panels[[1]] | panels[[2]] | panels[[3]] | panels[[4]]) /
      (panels[[5]] | panels[[6]]) /
      (panels[[7]] | panels[[8]] | panels[[9]] | panels[[10]])) +
    plot_layout(guides = "collect", axis_titles = "collect_y") +
    plot_annotation(title = "C.",
                    theme = theme(plot.title = element_text(face = "bold", size = 10))) &
    coord_cartesian(ylim = ylim_common)   # `&` applies to every sub-plot
}

# --------------------------------------------------------------------------
# AIG figure
# --------------------------------------------------------------------------

hist_AIG <- make_metric_hist(
  df, "AIG_year_area1", title = "A.",
  x_label = "AIG per year in Area 1\n(cases per year per 1000 individuals)",
  binwidth = 5
)

figure_AIG <- (hist_AIG + plot_sobol_indices(sobol_AIG)) /
  free(wrap_elements(build_profile_panels("AIG_year_area1", "AIG per year in Area 1")),
       side = "l") +
  plot_layout(heights = c(0.4, 1))

ggsave(figure_path("figure3_sensitivity_AIG.png"), figure_AIG, width = 11.4, height = 10.8)

# --------------------------------------------------------------------------
# AIGR figure
# --------------------------------------------------------------------------

hist_AIGR <- make_metric_hist(
  df, "AIGR_year_area1", title = "A.",
  x_label = "AIGR per year in Area 1",
  percent = TRUE
)

figure_AIGR <- (hist_AIGR + plot_sobol_indices(sobol_AIGR)) /
  free(wrap_elements(build_profile_panels("AIGR_year_area1", "AIGR per year in Area 1",
                                          y_percent = TRUE)),
       side = "l") +
  plot_layout(heights = c(0.4, 1))

ggsave(figure_path("figureS3_sensitivity_AIGR.png"), figure_AIGR, width = 11.4, height = 10.8)

# --------------------------------------------------------------------------
# Near-elimination regimes
# --------------------------------------------------------------------------

df_RC1_below1 <- df[df$RC_1 < 1, ]
df_RC1_around1 <- df[df$RC_1 >= 0.9 & df$RC_1 <= 1.1, ]

cat("Simulations with RC_1 < 1         :", nrow(df_RC1_below1), "\n")
cat("Simulations with RC_1 in [0.9,1.1]:", nrow(df_RC1_around1), "\n")

figure_S4 <-
  (make_metric_hist(df_RC1_below1, "AIG_year_area1",
                    title = "AIG per year in Area 1",
                    subtitle = expression(R[C] ~ "in area 1 < 1"),
                    x_label = "AIG per year in Area 1", binwidth = 5) |
     make_metric_hist(df_RC1_around1, "AIG_year_area1",
                      title = "AIG per year in Area 1",
                      subtitle = expression(R[C] ~ "in area 1 in [0.9, 1.1]"),
                      x_label = "AIG per year in Area 1", binwidth = 5)) /
  (make_metric_hist(df_RC1_below1, "AIGR_year_area1",
                    title = "AIGR per year in Area 1",
                    subtitle = expression(R[C] ~ "in area 1 < 1"),
                    x_label = "AIGR per year in Area 1", percent = TRUE) |
     make_metric_hist(df_RC1_around1, "AIGR_year_area1",
                      title = "AIGR per year in Area 1",
                      subtitle = expression(R[C] ~ "in area 1 in [0.9, 1.1]"),
                      x_label = "AIGR per year in Area 1", percent = TRUE))

ggsave(figure_path("figureS4_near_elimination.png"), figure_S4)

# --------------------------------------------------------------------------
# Recovery time against intervention duration
#
# Both parameters set the time scale of the model and interact strongly (they
# are the strongest pair in the second-order Sobol indices), so their joint
# effect is worth showing directly rather than through two marginal profiles.
# --------------------------------------------------------------------------

heatmap_df <- df %>%
  mutate(rinv_bin = cut(rinv, breaks = c(seq(60, 190, by = 10), 200 + 1e-10), right = FALSE,
                        include.lowest = TRUE)) %>%
  filter(!is.na(rinv_bin)) %>%
  group_by(rinv_bin, time_intervention) %>%
  summarise(mean_AIG = mean(AIG_year_area1, na.rm = TRUE), .groups = "drop")

levels(heatmap_df$rinv_bin)[length(levels(heatmap_df$rinv_bin))] <- "[190,200]"

ggplot(heatmap_df, aes(x = rinv_bin, y = factor(time_intervention), fill = mean_AIG)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Mean AIG") +
  labs(x = "Recovery time (days)", y = "Duration of intervention") +
  theme_bw() +
  theme(
    panel.border = element_rect(
      colour = "grey47",
      fill = NA,
      linewidth = 0.8
    ),
    axis.line = element_line(colour = "grey47"),
    axis.ticks = element_line(colour = "grey47"),
    axis.text = element_text(colour = "grey47"),
    axis.title = element_text(colour = "grey47"),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "grey47"
    )
  )

ggsave(figure_path("figureS6_recovery_duration_heatmap.png"),
       width = 9, height = 6, dpi = 300)

# --------------------------------------------------------------------------
# Summary quantiles
# --------------------------------------------------------------------------

# `n` is reported alongside the quantiles so that any remaining difference in
# sample size between the AIG and AIGR rows is visible in the table itself
# rather than hidden by the is.finite() filter. After
# load_simulation_database(), the four rows should share the same n.
stats <- df %>%
  select(AIG_year_area1, AIGR_area1, AIG_year_area2, AIGR_area2) %>%
  tidyr::pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  filter(is.finite(value)) %>%
  group_by(variable) %>%
  summarise(n = n(),
            min = min(value),
            max = max(value),
            q1 = quantile(value, 0.01),
            q99 = quantile(value, 0.99),
            q25 = quantile(value, 0.25),
            q75 = quantile(value, 0.75),
            median = median(value),
            mean = mean(value),
            .groups = "drop")

# Formatted for direct transcription into Supplementary Table S1
ratio_rows <- grepl("^AIGR", stats$variable)
# `n` is a count, not a metric value: it must not be scaled by 100 or suffixed
# with a percent sign on the AIGR rows.
value_cols <- setdiff(names(stats), c("variable", "n"))

stats_formatted <- stats
stats_formatted[value_cols] <- lapply(value_cols, function(col) {
  ifelse(ratio_rows,
         paste0(sprintf("%.1f", 100 * stats[[col]]), "%"),
         sprintf("%.1f", stats[[col]]))
})

write.csv(stats_formatted, data_path("tableS1_metric_quantiles.csv"), row.names = FALSE)