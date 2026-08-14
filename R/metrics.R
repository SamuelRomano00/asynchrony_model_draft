#' Asynchrony-induced growth metrics
#'
#' AIG is the absolute excess of cases caused by running the intervention
#' asynchronously rather than synchronously, over a study window of
#' `nb_studied_cycles` on/off cycles. AIGR is the same quantity relative to the
#' synchronous total. Both are also reported per year of that window.

#' Compare two annual-incidence series
#'
#' @param incidence_async,incidence_sync Numeric vectors of annual incidence,
#'   with patch 1's years followed by patch 2's. This is exactly
#'   `annual_incidence$value` ordered by `variable` then `t`. Each patch's block
#'   is located by integer position, not by label, so the ordering matters.
#' @param n_days Total number of days simulated.
#' @param start_interv Day the first intervention period began.
#' @param length_intervention Years per on/off period.
#' @param nb_studied_cycles Number of full intervention cycles defining the
#'   study window. One cycle is one "on" phase followed by one "off" phase, so
#'   `2 * length_intervention` years.
#' @return A named numeric vector with, per area and overall: the summed cases
#'   under each scenario (`As*`, `S*`), `AIGR*` (relative excess), `AIG*`
#'   (absolute excess, rounded), and the `*_year` variants of both.
compute_metrics <- function(incidence_async, incidence_sync, n_days,
                            start_interv, length_intervention,
                            nb_studied_cycles) {

  # Entries per patch, including the t = 0 placeholder; hence start_area2 is
  # the same offset shifted by one whole block.
  B <- floor(n_days / 365) + 1
  start_area1 <- floor(start_interv / 365) + 2   # first year actually under intervention
  start_area2 <- B + start_area1
  studied_duration_intervention <- ceiling(2 * length_intervention * nb_studied_cycles)
  end_area1 <- start_area1 + studied_duration_intervention - 1
  end_area2 <- start_area2 + studied_duration_intervention - 1

  sum_async_area1 <- sum(incidence_async[start_area1:end_area1])
  sum_async_area2 <- sum(incidence_async[start_area2:end_area2])
  sum_async <- sum_async_area1 + sum_async_area2

  sum_sync_area1 <- sum(incidence_sync[start_area1:end_area1])
  sum_sync_area2 <- sum(incidence_sync[start_area2:end_area2])
  sum_sync <- sum_sync_area1 + sum_sync_area2

  AIGR <- (sum_async - sum_sync) / sum_sync
  AIGR_area1 <- (sum_async_area1 - sum_sync_area1) / sum_sync_area1
  AIGR_area2 <- (sum_async_area2 - sum_sync_area2) / sum_sync_area2

  AIGR_year <- AIGR / studied_duration_intervention
  AIGR_year_area1 <- AIGR_area1 / studied_duration_intervention
  AIGR_year_area2 <- AIGR_area2 / studied_duration_intervention

  AIG <- round(sum_async - sum_sync)
  AIG_area1 <- round(sum_async_area1 - sum_sync_area1)
  AIG_area2 <- round(sum_async_area2 - sum_sync_area2)

  AIG_year <- round(AIG / studied_duration_intervention)
  AIG_year_area1 <- round(AIG_area1 / studied_duration_intervention)
  AIG_year_area2 <- round(AIG_area2 / studied_duration_intervention)

  c("As" = sum_async, "AsA1" = sum_async_area1, "AsA2" = sum_async_area2,
    "S" = sum_sync, "SA1" = sum_sync_area1, "SA2" = sum_sync_area2,
    "AIGR" = AIGR, "AIGR_area1" = AIGR_area1, "AIGR_area2" = AIGR_area2,
    "AIGR_year" = AIGR_year, "AIGR_year_area1" = AIGR_year_area1,
    "AIGR_year_area2" = AIGR_year_area2,
    "AIG" = AIG, "AIG_area1" = AIG_area1, "AIG_area2" = AIG_area2,
    "AIG_year" = AIG_year, "AIG_year_area1" = AIG_year_area1,
    "AIG_year_area2" = AIG_year_area2)
}

#' Metrics across stochastic replicates
#'
#' Applies [compute_metrics()] to each replicate separately, then summarises
#' the resulting distribution.
#'
#' @param sim_async,sim_sync Outputs of [run_SIS_stochastic()] for the
#'   asynchronous and synchronous scenarios.
#' @inheritParams compute_metrics
#' @return A list with
#'   \describe{
#'     \item{per_rep}{One row per replicate, all metrics plus `rep`.}
#'     \item{summary}{Mean and 2.5/97.5 percentiles of each metric.}
#'     \item{n_nan_AIGR}{Number of replicates with a non-finite AIGR.}
#'   }
#'
#' A replicate that has already gone extinct contributes zero synchronous cases
#' in the window, which makes AIGR undefined for that replicate while AIG, a
#' plain difference, stays well defined. Those replicates are excluded from the
#' AIGR summaries rather than allowed to blank out the whole row, and counted in
#' `n_nan_AIGR`.
compute_metrics_stochastic <- function(sim_async, sim_sync, n_days, start_interv,
                                       length_intervention, nb_studied_cycles) {

  n_reps <- length(unique(sim_sync$annual_incidence$rep))

  per_rep <- lapply(seq_len(n_reps), function(k) {
    inc_async <- sim_async$annual_incidence %>%
      filter(rep == k) %>% arrange(variable, t) %>% pull(value)
    inc_sync <- sim_sync$annual_incidence %>%
      filter(rep == k) %>% arrange(variable, t) %>% pull(value)
    m <- compute_metrics(inc_async, inc_sync, n_days, start_interv,
                         length_intervention, nb_studied_cycles)
    as.data.frame(t(m)) %>% mutate(rep = k)
  })
  per_rep <- do.call(rbind, per_rep)

  n_nan <- sum(!is.finite(per_rep$AIGR))
  if (n_nan > 0) {
    message(n_nan, " / ", n_reps, " replicate(s) had a non-finite AIGR ",
            "(no synchronous case in the metric window, usually because that ",
            "replicate had already gone extinct). AIG is unaffected.")
  }

  summary <- per_rep %>%
    mutate(across(-rep, ~replace(.x, !is.finite(.x), NA))) %>%
    select(-rep) %>%
    summarise(across(everything(),
                     list(mean = ~mean(.x, na.rm = TRUE),
                          q2.5 = ~unname(quantile(.x, 0.025, na.rm = TRUE)),
                          q97.5 = ~unname(quantile(.x, 0.975, na.rm = TRUE))),
                     .names = "{.col}_{.fn}"))

  list(per_rep = per_rep, summary = summary, n_nan_AIGR = n_nan)
}

#' Summarise a long-format series across replicates
#'
#' @param df Data frame with a `value` column, a `rep` column and the grouping
#'   columns named in `group_cols`.
#' @param group_cols Character vector of columns to group by.
#' @return One row per group, with `mean`, `median` and the 10th, 25th, 75th
#'   and 90th percentiles of `value` across replicates.
summarise_replicates <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(mean = mean(value),
              median = median(value),
              q10 = quantile(value, 0.10),
              q25 = quantile(value, 0.25),
              q75 = quantile(value, 0.75),
              q90 = quantile(value, 0.90),
              .groups = "drop")
}

#' Cumulative proportion of replicates with interrupted transmission
#'
#' Transmission is considered interrupted in a replicate from the first year
#' after its last recorded case, and never if its final observed year still has
#' cases.
#'
#' @param incidence Long-format annual incidence with columns `t`, `value`,
#'   `rep`, `area` and `scenario`.
#' @return A data frame with `area`, `scenario`, `year` and `proportion`.
compute_interruption_curve <- function(incidence) {

  definitive_year <- incidence %>%
    group_by(rep, area, scenario) %>%
    arrange(t) %>%
    summarise(
      interruption_year = {
        yr <- t / 365
        nz <- which(value != 0)
        if (length(nz) == 0) {
          min(yr)                # never any case: interrupted from the start
        } else if (max(nz) == length(value)) {
          Inf                    # cases in the last observed year: not confirmed
        } else {
          yr[max(nz)] + 1        # first year after the last case
        }
      },
      .groups = "drop"
    )

  years <- sort(unique(incidence$t)) / 365

  definitive_year %>%
    group_by(area, scenario) %>%
    reframe(year = years,
            proportion = sapply(years, function(y) mean(interruption_year <= y)))
}
