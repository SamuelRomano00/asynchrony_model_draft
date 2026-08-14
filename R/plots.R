#' Figures
#'
#' Shared plotting conventions across every figure in the paper: the
#' synchronous scenario in black, the asynchronous one in orange, and the AIG
#' shown as the hatched area between the two incidence curves.

#' Colours identifying the two intervention schedules
SCENARIO_COLOURS <- c("Synchronous" = "black", "Asynchronous" = "darkorange3")

#' Sampled parameters, in the order used on every axis that lists them
PARAM_ORDER <- c("rinv", "time_intervention", "R0_1", "R0_2",
                 "omega_1", "omega_2", "p_12", "p_21")

#' Axis labels for the sampled parameters
#'
#' Uses `expression()` for the mobility parameters, so `p_12` renders as a
#' subscripted symbol. Valid wherever ggplot2 accepts a labels vector; use
#' [PARAM_LABELS_PLAIN] where only character strings are allowed.
PARAM_LABELS <- c(
  "rinv" = "Recovery time",
  "time_intervention" = "Length of intervention",
  "R0_1" = "R0 in area 1",
  "R0_2" = "R0 in area 2",
  "omega_1" = "Intervention efficiency in area 1",
  "omega_2" = "Intervention efficiency in area 2",
  "p_12" = expression(p["1,2"]),
  "p_21" = expression(p["2,1"])
)

#' Plain-text version of [PARAM_LABELS]
#'
#' `rpart.plot()` draws node labels from a character column of the fitted
#' object, which cannot hold expressions.
PARAM_LABELS_PLAIN <- c(
  "rinv" = "Duration of recovery",
  "time_intervention" = "Duration of intervention",
  "R0_1" = "R0 in area 1",
  "R0_2" = "R0 in area 2",
  "omega_1" = "Intervention efficiency in area 1",
  "omega_2" = "Intervention efficiency in area 2",
  "p_12" = "p_1,2",
  "p_21" = "p_2,1"
)

#' Horizontal segments marking the on/off intervention schedule
#'
#' @param xmin,xmax Range to cover, in years.
#' @param y Height at which to draw the segments.
#' @param length,space Length of each "on" segment and of the gap after it, in
#'   years.
#' @param colours Value of the `col` column, used to map the segments to a
#'   scenario.
#' @return A data frame of segments for `geom_segment()`.
make_intervention_segments <- function(xmin, xmax, y, length, space, colours) {
  seqs <- seq(xmin, xmax, by = length + space)
  segs <- do.call(rbind, lapply(seqs, function(start) {
    ends <- start + length
    if (ends > xmax) return(NULL)
    data.frame(x = start, xend = ends, y = y, yend = y)
  }))
  segs$col <- rep(colours, length.out = nrow(segs))
  segs
}

#' One incidence panel: two areas, two scenarios, AIG highlighted
#'
#' @param incidence_i Long-format annual incidence with columns `t`, `value`,
#'   `area` (`"Area1"`/`"Area2"`) and `scenario`.
#' @param metrics_labels_i Data frame with `area` and a `value` column holding
#'   the already-formatted label text.
#' @param length_intervention_i Years per on/off period, used to draw the
#'   schedule segments.
#' @param year_start_plot,year_end_plot Plotting window, in years.
#' @param max_y,gap Height of the asynchronous schedule segments and vertical
#'   offset of the synchronous ones. Either a scalar applied to both areas, or a
#'   named vector indexed by `"Area1"`/`"Area2"`.
#' @param show_aig_area Whether to hatch the area between the two curves.
#' @param text_scale Multiplier applied to every text size, to compensate for
#'   panels being shrunk when assembled into a grid.
#' @param label_y Height of the metric label; defaults to 95% of `max_y`.
#' @param title_i Optional panel title.
#' @return A ggplot object.
plot_incidence_panel <- function(incidence_i, metrics_labels_i, length_intervention_i,
                                 year_start_plot, year_end_plot, max_y, gap,
                                 show_aig_area = TRUE, text_scale = 1,
                                 label_y = NULL, title_i = NULL, start_year = 1) {

  # A scalar applies to both areas; a named vector is indexed by area. Areas
  # plotted on free scales need their own heights, otherwise the schedule
  # segments and the label land outside the shorter of the two panels.
  get_area_val <- function(v, area) {
    if (length(v) == 1) unname(v) else unname(v[[area]])
  }
  max_y1 <- get_area_val(max_y, "Area1"); gap1 <- get_area_val(gap, "Area1")
  max_y2 <- get_area_val(max_y, "Area2"); gap2 <- get_area_val(gap, "Area2")

  # Area 1's asynchronous schedule is never phase-shifted; area 2's is shifted
  # by exactly one period, which is what "asynchronous" means here.
  segments_all <- bind_rows(
    make_intervention_segments(start_year, year_end_plot + 1, max_y1,
                               length_intervention_i, length_intervention_i,
                               "Asynchronous") %>% mutate(area = "Area1"),
    make_intervention_segments(start_year, year_end_plot + 1, max_y1 + gap1,
                               length_intervention_i, length_intervention_i,
                               "Synchronous") %>% mutate(area = "Area1"),
    make_intervention_segments(start_year + length_intervention_i,
                               length_intervention_i + year_end_plot + 1, max_y2,
                               length_intervention_i, length_intervention_i,
                               "Asynchronous") %>% mutate(area = "Area2"),
    make_intervention_segments(start_year, year_end_plot + 1, max_y2 + gap2,
                               length_intervention_i, length_intervention_i,
                               "Synchronous") %>% mutate(area = "Area2")
  )

  # Filter to the plotting window BEFORE pivoting: geom_ribbon_pattern() can
  # fail on points far outside the final xlim(), unlike a plain geom_ribbon
  # which clips cleanly afterwards.
  wide_incidence <- incidence_i %>%
    filter(t / 365 >= year_start_plot, t / 365 <= year_end_plot) %>%
    tidyr::pivot_wider(names_from = scenario, values_from = value) %>%
    filter(!is.na(Synchronous), !is.na(Asynchronous))

  metrics_labels_i <- metrics_labels_i %>%
    mutate(.label_y = sapply(area, function(a) {
      if (is.null(label_y)) get_area_val(max_y, a) * 0.95 else get_area_val(label_y, a)
    }))

  incidence_i %>%
    ggplot(aes(x = t / 365, y = value)) +
    # fill = NA so the hatching never masks what sits underneath.
    {if (show_aig_area) geom_ribbon_pattern(
      data = wide_incidence,
      aes(x = t / 365,
          ymin = pmin(Synchronous, Asynchronous),
          ymax = pmax(Synchronous, Asynchronous)),
      inherit.aes = FALSE,
      pattern = "stripe", pattern_density = 0.05, pattern_spacing = 0.025,
      pattern_fill = "grey55", pattern_colour = NA, fill = NA, colour = NA
    )} +
    geom_segment(data = segments_all,
                 aes(x = x, xend = xend, y = y, yend = yend, color = col),
                 inherit.aes = FALSE, linewidth = 2) +
    geom_line(aes(color = scenario), linewidth = 1.5) +
    facet_wrap(area ~ ., scales = "free",
               labeller = labeller(area = c("Area1" = "Area 1", "Area2" = "Area 2"))) +
    labs(y = "Annual incidence", x = "Years", color = "Scenario", title = title_i) +
    xlim(year_start_plot, year_end_plot) +
    # Generous headroom above the schedule segments, so they, the hatched area
    # and the facet strip do not end up touching.
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.12))) +
    scale_color_manual(values = SCENARIO_COLOURS) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12 * text_scale),
      strip.text = element_text(size = 14 * text_scale,
                                margin = ggplot2::margin(t = 5, b = 5)),
      axis.title = element_text(size = 14 * text_scale),
      legend.text = element_text(size = 18 * text_scale),
      legend.title = element_text(size = 18 * text_scale),
      plot.title = element_text(size = 16 * text_scale, face = "bold", hjust = 0.5,
                                margin = ggplot2::margin(t = 2, b = 12)),
      plot.margin = ggplot2::margin(t = 10, r = 12, b = 10, l = 6)
    ) +
    geom_label(
      # vjust = 1 anchors the TOP of the label: with AIG = TRUE the text grows
      # from one line to two, and a centred label would grow upwards into the
      # schedule segments.
      data = metrics_labels_i,
      aes(x = year_end_plot - 5, y = .label_y, label = value, group = area),
      inherit.aes = FALSE,
      size = 5 * text_scale, hjust = 1, vjust = 1,
      color = "grey30", fill = "white", label.size = 0.5, fontface = "bold"
    )
}

#' Deterministic incidence figure, optionally stacked over several panels
#'
#' @param incidence Long-format annual incidence with columns `t`, `value`,
#'   `area` and `scenario`, plus the `facet_rows` column when one is given.
#' @param length_intervention Years per on/off period: a scalar, or one value
#'   per level of `facet_rows` when they differ across panels.
#' @param year_start_plot,year_end_plot Plotting window, in years.
#' @param max_y,gap Passed to [plot_incidence_panel()].
#' @param metrics_labels Data frame with `area`, `AIG_val`, `AIG_year_val` and
#'   the `facet_rows` column when one is given. The label text is built here.
#' @param facet_rows Name of the column identifying each stacked panel, e.g.
#'   `"model"` or `"length_intervention"`. `NULL` for a single panel.
#' @param show_aig_area Whether to hatch the area between the two curves.
#' @param AIG If `TRUE`, the label shows the total AIG above the AIG per year.
#' @param panel_heights Relative heights passed to `patchwork::wrap_plots()`.
#' @param panel_titles Whether to title each panel with its `facet_rows` level.
#' @return A ggplot or patchwork object.
#'
#' Panels are assembled with patchwork rather than a single
#' `facet_grid(facet_rows ~ area)` because a 2D facetting with free scales on
#' both axes makes `ggpattern::geom_ribbon_pattern()` fail on its hatching grid.
visualize_incidence <- function(incidence, length_intervention,
                                year_start_plot, year_end_plot, max_y, gap,
                                metrics_labels, facet_rows = NULL,
                                show_aig_area = TRUE, AIG = FALSE,
                                panel_heights = NULL, panel_titles = TRUE, 
                                start_year = 1) {

  metrics_labels <- metrics_labels %>%
    mutate(value = if (AIG) {
      paste0("AIG: ", round(AIG_val), "\nAIG per year: ", round(AIG_year_val))
    } else {
      paste0("AIG per year: ", round(AIG_year_val))
    })

  if (is.null(facet_rows)) {
    return(plot_incidence_panel(incidence, metrics_labels,
                                length_intervention_i = length_intervention,
                                year_start_plot = year_start_plot,
                                year_end_plot = year_end_plot,
                                max_y = max_y, gap = gap,
                                start_year = 1,
                                show_aig_area = show_aig_area))
  }

  levels_i <- if (is.factor(incidence[[facet_rows]])) {
    levels(droplevels(incidence[[facet_rows]]))
  } else {
    unique(incidence[[facet_rows]])
  }

  if (length(length_intervention) == 1) {
    length_intervention <- rep(length_intervention, length(levels_i))
  }
  if (length(length_intervention) != length(levels_i)) {
    stop("`length_intervention` must be a scalar or have one value per level of `",
         facet_rows, "` (", length(levels_i), " expected, ",
         length(length_intervention), " given).")
  }

  plots <- Map(function(lv, li) {
    plot_incidence_panel(incidence[incidence[[facet_rows]] == lv, ],
                         metrics_labels[metrics_labels[[facet_rows]] == lv, ],
                         length_intervention_i = li,
                         year_start_plot = year_start_plot,
                         year_end_plot = year_end_plot,
                         max_y = max_y, gap = gap,
                         start_year = 1,
                         show_aig_area = show_aig_area,
                         title_i = if (panel_titles) as.character(lv) else NULL)
  }, levels_i, length_intervention)

  patchwork::wrap_plots(plots, ncol = 1, guides = "collect", heights = panel_heights) &
    theme(legend.position = "bottom")
}

#' Stochastic incidence figure with replicate ribbons
#'
#' Same conventions as [plot_incidence_panel()], with the single deterministic
#' trajectory replaced by the mean across replicates and two uncertainty bands
#' (interquartile and decile-to-decile).
#'
#' @param incidence Long-format annual incidence with columns `t`, `value`,
#'   `rep`, `area` and `scenario`.
#' @param length_intervention Years per on/off period.
#' @param year_start_plot,year_end_plot Plotting window, in years.
#' @param gap Vertical offset between the two schedules' segments.
#' @param metrics_labels Data frame with `area`, `AIG_val` and `AIG_year_val`.
#' @param show_aig_area Whether to hatch the area between the two mean curves.
#' @param AIG If `TRUE`, the label shows the total AIG above the AIG per year.
#' @return A ggplot object.
visualize_incidence_stochastic <- function(incidence, length_intervention,
                                           year_start_plot, year_end_plot, gap,
                                           metrics_labels, show_aig_area = FALSE,
                                           AIG = FALSE, start_year = 1) {

  metrics_labels <- metrics_labels %>%
    mutate(value = if (AIG) {
      paste0("AIG: ", round(AIG_val), "\nAIG per year: ", round(AIG_year_val))
    } else {
      paste0("AIG per year: ", round(AIG_year_val))
    })

  max_y <- max(incidence$value)

  segments_all <- bind_rows(
    make_intervention_segments(start_year, year_end_plot + 1, max_y,
                               length_intervention, length_intervention,
                               "Asynchronous") %>% mutate(area = "Area1"),
    make_intervention_segments(start_year, year_end_plot + 1, max_y + gap,
                               length_intervention, length_intervention,
                               "Synchronous") %>% mutate(area = "Area1"),
    make_intervention_segments(start_year + length_intervention,
                               length_intervention + year_end_plot + 1, max_y,
                               length_intervention, length_intervention,
                               "Asynchronous") %>% mutate(area = "Area2"),
    make_intervention_segments(start_year, year_end_plot + 1, max_y + gap,
                               length_intervention, length_intervention,
                               "Synchronous") %>% mutate(area = "Area2")
  )

  summary_df <- summarise_replicates(incidence, c("t", "area", "scenario"))
  sync <- summary_df %>% filter(scenario == "Synchronous")
  async <- summary_df %>% filter(scenario == "Asynchronous")

  if (show_aig_area) {
    mean_wide <- summary_df %>%
      filter(t / 365 >= year_start_plot, t / 365 <= year_end_plot) %>%
      select(t, area, scenario, mean) %>%
      tidyr::pivot_wider(names_from = scenario, values_from = mean)
  }

  ggplot(summary_df, aes(x = t / 365)) +
    geom_ribbon(data = sync, aes(ymin = q10, ymax = q90), fill = "grey85", alpha = 0.8) +
    geom_ribbon(data = sync, aes(ymin = q25, ymax = q75), fill = "grey55", alpha = 0.9) +
    geom_ribbon(data = async, aes(ymin = q10, ymax = q90), fill = "navajowhite", alpha = 0.8) +
    geom_ribbon(data = async, aes(ymin = q25, ymax = q75), fill = "tan2", alpha = 0.85) +
    {if (show_aig_area) geom_ribbon_pattern(
      data = mean_wide,
      aes(x = t / 365,
          ymin = pmin(Synchronous, Asynchronous),
          ymax = pmax(Synchronous, Asynchronous)),
      inherit.aes = FALSE,
      pattern = "stripe", pattern_density = 0.05, pattern_spacing = 0.025,
      pattern_fill = "grey55", pattern_colour = NA, fill = NA, colour = NA
    )} +
    geom_line(aes(y = mean, color = scenario), linewidth = 0.9) +
    geom_segment(data = segments_all,
                 aes(x = x, xend = xend, y = y, yend = yend, color = col),
                 inherit.aes = FALSE, linewidth = 2) +
    scale_color_manual(values = SCENARIO_COLOURS) +
    facet_wrap(area ~ ., scales = "free",
               labeller = labeller(area = c("Area1" = "Area 1", "Area2" = "Area 2"))) +
    labs(y = "Annual incidence", x = "Years", color = "Scenario") +
    xlim(year_start_plot, year_end_plot) +
    ylim(0, NA) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(linewidth = 1.2))) +
    geom_label(data = metrics_labels,
               aes(x = year_end_plot - 5, y = max_y * 0.95, label = value, group = area),
               inherit.aes = FALSE,
               size = 4.5, hjust = 1, vjust = 1,
               color = "grey30", fill = "white", label.size = 0.5, fontface = "bold")
}

#' Monotone spline through an interruption curve
#'
#' @param df Data frame with `year` and `proportion`.
#' @param n_out Number of interpolated points.
#' @param x_min,x_max Optional window, clipped to the observed range.
#' @return A data frame with `year` and `proportion`, interpolated.
smooth_monotone_curve <- function(df, n_out = 200, x_min = NULL, x_max = NULL) {
  df <- df[order(df$year), ]
  f <- stats::splinefun(df$year, df$proportion, method = "monoH.FC")
  lo <- max(min(df$year), if (is.null(x_min)) min(df$year) else x_min)
  hi <- min(max(df$year), if (is.null(x_max)) max(df$year) else x_max)
  yg <- seq(lo, hi, length.out = n_out)
  data.frame(year = yg, proportion = pmin(pmax(f(yg), 0), 1))
}

#' Proportion of replicates with interrupted transmission, over time
#'
#' @param incidence Long-format annual incidence with columns `t`, `value`,
#'   `rep`, `area` and `scenario`.
#' @param year_start_plot,year_end_plot Plotting window, in years.
#' @param smooth Whether to draw a monotone spline instead of the raw steps.
#' @return A ggplot object.
visualize_interruption_stochastic <- function(incidence, year_start_plot, year_end_plot,
                                              smooth = FALSE) {

  curve_df <- compute_interruption_curve(incidence)

  p <- ggplot(curve_df, aes(x = year, y = proportion, color = scenario))

  if (smooth) {
    smoothed <- curve_df %>%
      group_by(area, scenario) %>%
      group_modify(~smooth_monotone_curve(.x, x_min = year_start_plot,
                                          x_max = year_end_plot)) %>%
      ungroup()
    p <- p + geom_line(data = smoothed,
                       aes(x = year, y = proportion, color = scenario),
                       linewidth = 0.9)
  } else {
    p <- p + geom_line(linewidth = 0.9) + geom_point(size = 1.5)
  }

  p +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    scale_color_manual(values = SCENARIO_COLOURS) +
    facet_wrap(~area, labeller = labeller(area = c("Area1" = "Area 1", "Area2" = "Area 2"))) +
    labs(x = "Years", y = "Prop. interrupted", color = "Scenario") +
    scale_x_continuous(limits = c(year_start_plot, year_end_plot),
                       breaks = seq(0, year_end_plot, by = 10),
                       minor_breaks = seq(0, year_end_plot, by = 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = sort(c(seq(0, 1, 0.25), 0.95))) +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
}

#' A metric against one input parameter, across the whole simulation database
#'
#' Bins the parameter and summarises the metric's distribution within each bin.
#'
#' @param data Data frame containing `param` and `metric`, e.g. the output of
#'   [metrics_computation()].
#' @param param Name of the x-axis column, as a string.
#' @param metric Name of the y-axis column, as a string.
#' @param bin_width Width of the bins the parameter is cut into.
#' @param name_param x-axis label; may be a character string or an
#'   `expression()`.
#' @param y_name Whether to label the y axis, useful when only one panel of a
#'   row should carry the shared label.
#' @param y_label Text of the y-axis label when `y_name` is `TRUE`.
#' @param y_percent Whether to format the y axis as a percentage. AIGR is a
#'   ratio and reads better in percent; AIG is a count and does not.
#' @return A ggplot object.
#'
#' The trend line is the MEDIAN, and the bands are quantiles, so a single
#' extreme simulation in a bin cannot drag the summary away from where the bulk
#' of that bin's values sit. That protection matters for heavy-tailed metrics
#' such as AIGR, whose ratio can spike near the elimination threshold.
plot_metric_by_parameter <- function(data, param, metric, bin_width, name_param,
                                     y_name = FALSE, y_label = "AIG per year in Area 1",
                                     y_percent = FALSE) {

  summary_df <- data %>%
    mutate(param_val = .data[[param]], metric_val = .data[[metric]]) %>%
    mutate(param_bin = cut(param_val,
                           breaks = seq(floor(min(param_val)),
                                        ceiling(max(param_val)), by = bin_width),
                           include.lowest = TRUE)) %>%
    group_by(param_bin) %>%
    summarise(param_mid = mean(param_val),
              I_median = median(metric_val, na.rm = TRUE),
              I_d1 = quantile(metric_val, 0.1, na.rm = TRUE),
              I_d9 = quantile(metric_val, 0.9, na.rm = TRUE),
              I_q1 = quantile(metric_val, 0.25, na.rm = TRUE),
              I_q3 = quantile(metric_val, 0.75, na.rm = TRUE),
              .groups = "drop")

  p <- ggplot(summary_df, aes(x = param_mid)) +
    geom_ribbon(aes(ymin = I_d1, ymax = I_d9), fill = "lightblue", alpha = 0.7) +
    geom_ribbon(aes(ymin = I_q1, ymax = I_q3), fill = "steelblue", alpha = 0.4) +
    geom_line(aes(y = I_median), color = "blue", linewidth = 1) +
    labs(x = name_param) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))

  p <- if (y_name) p + labs(y = y_label) else p + theme(axis.title.y = element_blank())
  if (y_percent) p <- p + scale_y_continuous(labels = scales::percent)

  p
}
