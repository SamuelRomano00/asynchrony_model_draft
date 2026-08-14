#' Publication tables

# ---------------------------------------------------------------------------
# Reading transmission interruption off the replicate curves
#
# Both helpers below are table-building concerns rather than metrics, so they
# live here rather than in metrics.R, next to the only function that uses them.
# ---------------------------------------------------------------------------

#' Year at which a scenario reaches transmission interruption in an area
#'
#' The same 95% threshold that is drawn on the interruption figures: the first
#' year by which at least 95% of replicates have interrupted for good.
#'
#' @param curve_df Output of [compute_interruption_curve()].
#' @param area_name,scenario_name Which curve to read.
#' @return The year, as a number, or `NA_real_` when the threshold is never
#'   reached within the simulated horizon.
#'
#' Returned as a number rather than a formatted string on purpose: the combined
#' "Both Areas" rows take a `max()` over the two areas, and on character values
#' that comparison is lexicographic, so `max("9", "16")` would be `"9"`.
interruption_year_95 <- function(curve_df, area_name, scenario_name) {
  d <- curve_df[curve_df$area == area_name &
                  curve_df$scenario == scenario_name &
                  curve_df$proportion >= 0.95, ]
  if (nrow(d) == 0) {
    return(NA_real_)
  }
  min(d$year)
}

#' Cumulated incidence per replicate, up to a given year
#'
#' @param incidence Long-format annual incidence with columns `t`, `value`,
#'   `rep`, `area` and `scenario`.
#' @param areas Areas to sum over: one of them, or both for a combined total.
#' @param scenario_name Which schedule to sum.
#' @param year Cut-off, in years, inclusive. `NA` returns `NULL`.
#' @param start_interv Day the first intervention period begins.
#' @return One total per replicate, or `NULL` when `year` is `NA`.
#'
#' The sum starts at the first year actually under intervention, exactly where
#' `compute_metrics()` starts its own window, so that this column and the AIG
#' describe the same period. Two things are excluded by the single `t >
#' start_interv` condition. The pre-intervention year, held at the endemic
#' equilibrium so that the starting level is visible on the figures, carries a
#' full year of untreated transmission and would dominate the total: in the
#' illustrative example it is on its own larger than everything accrued
#' afterwards. And the `t = 0` row, a placeholder carrying year 1's value again
#' so that each patch's block has one entry per year boundary for
#' `compute_metrics()` to index by position, would count that year twice.
cumulated_incidence_to_year <- function(incidence, areas, scenario_name, year,
                                        start_interv) {
  
  if (is.na(year)) {
    return(NULL)
  }
  
  incidence %>%
    filter(area %in% areas,
           scenario == scenario_name,
           t > start_interv,
           t / 365 <= year) %>%
    group_by(rep) %>%
    summarise(total = sum(value), .groups = "drop") %>%
    pull(total)
}

# ---------------------------------------------------------------------------
# Publication table
# ---------------------------------------------------------------------------

#' Summary table for one stochastic scenario
#'
#' Assembles the numbers reported in the paper's result tables: AIG and AIGR
#' with their 95% replicate intervals, the cumulated incidence accrued before
#' transmission interruption, and the time to that interruption.
#'
#' @param summary The `summary` element of [compute_metrics_stochastic()].
#' @param incidence Long-format annual incidence with columns `t`, `value`,
#'   `rep`, `area` and `scenario`.
#' @param start_interv Day the first intervention period begins, the same value
#'   the simulation was run with. Everything before it is pre-intervention time
#'   at the endemic equilibrium and is excluded from the totals.
#' @return A tibble with one block of two rows per area plus a combined block.
#'
#' The two halves of the table start at the same point, the first year under
#' intervention, but END differently, and the caption has to say so. AIG and
#' AIGR run to the close of the fixed metric window of three intervention
#' cycles, because that is what makes them comparable across scenarios. The
#' cumulated incidence runs to the year the scenario actually reaches
#' interruption, which is what "before transmission interruption" means and
#' which differs from one schedule to the next.
#'
#' Where a scenario never reaches the 95% threshold, there is no such year, and
#' the cumulated-incidence cell is left empty rather than falling back on the
#' metric window: a total labelled "before transmission interruption" for a
#' scenario that never interrupts would be meaningless.
build_publication_table <- function(summary, incidence, start_interv) {
  
  curve_df <- compute_interruption_curve(incidence)
  
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
  
  # Interruption years, per area and per schedule. A combined block counts as
  # interrupted only once BOTH areas are, and max() propagates the NA of an
  # area that never interrupts, which is the behaviour wanted here.
  year_of <- function(area_name, scenario_name) {
    interruption_year_95(curve_df, area_name, scenario_name)
  }
  years <- list(
    Area1 = c(Synchronous = year_of("Area1", "Synchronous"),
              Asynchronous = year_of("Area1", "Asynchronous")),
    Area2 = c(Synchronous = year_of("Area2", "Synchronous"),
              Asynchronous = year_of("Area2", "Asynchronous"))
  )
  years$Both <- c(
    Synchronous = max(years$Area1[["Synchronous"]], years$Area2[["Synchronous"]]),
    Asynchronous = max(years$Area1[["Asynchronous"]], years$Area2[["Asynchronous"]])
  )
  
  fmt_year <- function(year) {
    if (is.na(year)) "No interruption" else as.character(round(year))
  }
  
  #' Cumulated incidence accrued by the time the scenario interrupts
  #'
  #' Every replicate contributes, including the fewer than 5% that are still
  #' transmitting at that year: their cases were accrued before the scenario
  #' reached the threshold, and dropping them would bias the mean downwards.
  #' Replicates that interrupted earlier contribute nothing after their own
  #' last case, so they need no special handling.
  cumulated_cell <- function(areas, scenario_name, year) {
    totals <- cumulated_incidence_to_year(incidence, areas, scenario_name, year,
                                          start_interv)
    if (is.null(totals) || length(totals) == 0) {
      return("")
    }
    fmt_num(mean(totals),
            unname(quantile(totals, 0.025)),
            unname(quantile(totals, 0.975)))
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
      cumulated_cell("Area1", "Synchronous", years$Area1[["Synchronous"]]),
      cumulated_cell("Area1", "Asynchronous", years$Area1[["Asynchronous"]]),
      cumulated_cell("Area2", "Synchronous", years$Area2[["Synchronous"]]),
      cumulated_cell("Area2", "Asynchronous", years$Area2[["Asynchronous"]]),
      cumulated_cell(c("Area1", "Area2"), "Synchronous", years$Both[["Synchronous"]]),
      cumulated_cell(c("Area1", "Area2"), "Asynchronous", years$Both[["Asynchronous"]])
    ),
    `Time to transmission interruption (in years)` = c(
      fmt_year(years$Area1[["Synchronous"]]), fmt_year(years$Area1[["Asynchronous"]]),
      fmt_year(years$Area2[["Synchronous"]]), fmt_year(years$Area2[["Asynchronous"]]),
      fmt_year(years$Both[["Synchronous"]]), fmt_year(years$Both[["Asynchronous"]])
    )
  )
}