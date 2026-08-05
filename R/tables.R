#' Publication tables

#' Summary table for one stochastic scenario
#'
#' Assembles the numbers reported in the paper's result tables: AIG and AIGR
#' with their 95% replicate intervals, the cumulated incidence under each
#' schedule, and the time to transmission interruption.
#'
#' @param summary The `summary` element of [compute_metrics_stochastic()].
#' @param incidence Long-format annual incidence with columns `t`, `value`,
#'   `rep`, `area` and `scenario`, used for the interruption times.
#' @return A tibble with one block of two rows per area plus a combined block.
build_publication_table <- function(summary, incidence) {

  curve_df <- compute_interruption_curve(incidence)

  # "Interrupted" is read off the same 95% threshold drawn on the interruption
  # figures.
  time_to_95 <- function(area_name, scenario_name) {
    d <- curve_df %>%
      filter(area == area_name, scenario == scenario_name, proportion >= 0.95)
    if (nrow(d) == 0) return("No interruption")
    as.character(round(min(d$year)))
  }

  get <- function(col) summary[[col]][1]

  fmt_num <- function(mean, q_low, q_high, digits = 0) {
    paste0(round(mean, digits), " [", round(q_low, digits), "-", round(q_high, digits), "]")
  }
  fmt_pct <- function(mean, q_low, q_high, digits = 0) {
    paste0(round(mean * 100, digits), "% [", round(q_low * 100, digits), "%-",
           round(q_high * 100, digits), "%]")
  }
  fmt_metric <- function(prefix, fmt) {
    fmt(get(paste0(prefix, "_mean")),
        get(paste0(prefix, "_q2.5")),
        get(paste0(prefix, "_q97.5")))
  }

  tibble::tibble(
    Area = c("Area 1", "", "Area 2", "", "Both Areas", ""),
    `AIG` = c(fmt_metric("AIG_area1", fmt_num), "",
              fmt_metric("AIG_area2", fmt_num), "",
              fmt_metric("AIG", fmt_num), ""),
    `AIG per year` = c(fmt_metric("AIG_year_area1", fmt_num), "",
                       fmt_metric("AIG_year_area2", fmt_num), "",
                       fmt_metric("AIG_year", fmt_num), ""),
    `AIGR` = c(fmt_metric("AIGR_area1", fmt_pct), "",
               fmt_metric("AIGR_area2", fmt_pct), "",
               fmt_metric("AIGR", fmt_pct), ""),
    `AIGR per year` = c(fmt_metric("AIGR_year_area1", fmt_pct), "",
                        fmt_metric("AIGR_year_area2", fmt_pct), "",
                        fmt_metric("AIGR_year", fmt_pct), ""),
    Scenario = c("Synchronous", "Asynchronous", "Synchronous", "Asynchronous",
                 "Synchronous", "Asynchronous"),
    `Cumulated annual incidence per 1000 people before transmission interruption` = c(
      fmt_metric("SA1", fmt_num), fmt_metric("AsA1", fmt_num),
      fmt_metric("SA2", fmt_num), fmt_metric("AsA2", fmt_num),
      fmt_metric("S", fmt_num), fmt_metric("As", fmt_num)
    ),
    `Time to transmission interruption (in years)` = c(
      time_to_95("Area1", "Synchronous"), time_to_95("Area1", "Asynchronous"),
      time_to_95("Area2", "Synchronous"), time_to_95("Area2", "Asynchronous"),
      max(time_to_95("Area1", "Synchronous"), time_to_95("Area2", "Synchronous")),
      max(time_to_95("Area1", "Asynchronous"), time_to_95("Area2", "Asynchronous"))
    )
  )
}
