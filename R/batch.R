#' Batch simulation over parameter designs
#'
#' Each row of a design is an independent simulation, so these loops are run in
#' parallel with `mclapply()`. Forking is used deliberately: workers inherit the
#' already-compiled odin model, whereas a socket cluster would have to
#' recompile it in every worker. On Windows, `mclapply()` falls back to a single
#' core with a warning, so this stays cross-platform but is only parallel on
#' Linux and macOS.

#' Names of the metrics returned by [compute_metrics()], in order
METRIC_NAMES <- c("As", "AsA1", "AsA2", "S", "SA1", "SA2",
                  "AIGR", "AIGR_area1", "AIGR_area2",
                  "AIGR_year", "AIGR_year_area1", "AIGR_year_area2",
                  "AIG", "AIG_area1", "AIG_area2",
                  "AIG_year", "AIG_year_area1", "AIG_year_area2")

#' Pick a safe default number of workers
#'
#' Tries `SLURM_CPUS_PER_TASK`, then `SLURM_JOB_CPUS_PER_NODE`, then
#' `parallel::detectCores() - 1`.
#'
#' @return A positive integer.
#'
#' The last resort reports the node's total hardware cores, not what the job was
#' actually allocated, so it can badly oversubscribe a shared cluster. It warns
#' loudly for that reason: if you see that warning under SLURM,
#' `--cpus-per-task` was probably missing from the submission.
detect_ncores <- function() {
  slurm_cpus_per_task <- Sys.getenv("SLURM_CPUS_PER_TASK")
  if (nzchar(slurm_cpus_per_task)) {
    return(as.integer(slurm_cpus_per_task))
  }
  
  slurm_cpus_per_node <- Sys.getenv("SLURM_JOB_CPUS_PER_NODE")
  if (nzchar(slurm_cpus_per_node)) {
    # Can be formatted as e.g. "24(x2)" for a multi-node job; keep the leading
    # integer.
    parsed <- suppressWarnings(as.integer(sub("^([0-9]+).*", "\\1", slurm_cpus_per_node)))
    if (!is.na(parsed) && parsed > 0) {
      return(parsed)
    }
  }
  
  ncores <- max(1, parallel::detectCores() - 1)
  warning("detect_ncores(): no SLURM core count found, falling back to ",
          "parallel::detectCores() - 1 = ", ncores, ", which reflects the ",
          "node's total hardware cores rather than this job's allocation.",
          call. = FALSE)
  ncores
}

#' Set up and simulate one row of a parameter design
#'
#' @param row One-row data frame with columns `rinv`, `R0_1`, `R0_2`, `p_12`,
#'   `p_21`, `omega_1`, `omega_2` and `time_intervention`.
#' @param n_days Simulation horizon, in days.
#' @param start_interv Day the first intervention period begins.
#' @param nb_studied_cycles Number of on/off cycles in the metric window.
#' @param N Population sizes per patch.
#' @return A list with `x0`, `z0`, `r`, `metrics`, and `viable`, which is
#'   `FALSE` when the parameter set has no endemic equilibrium.
simulate_design_row <- function(row, n_days, start_interv, nb_studied_cycles, N,
                                skip_if_not_viable = FALSE) {
  
  # r is derived from the sampled infectious period: the designs carry `rinv`
  # (days) rather than a rate.
  r <- 1 / row$rinv
  R0 <- c(row$R0_1, row$R0_2)
  
  x0 <- compute_equilibrium_prevalence(R0, N, r, row$p_12, row$p_21)
  z0 <- x0 * r * 365 * 1000
  
  # Without an endemic equilibrium the trajectory is identically zero and the
  # metrics carry no information, so the two ODE integrations are pure waste.
  # metrics_computation() still wants them simulated, because that is what
  # produces the NaN AIGR its database records, hence the flag rather than an
  # unconditional skip.
  if (skip_if_not_viable && (x0[1] <= 0 || x0[2] <= 0)) {
    return(list(x0 = x0, z0 = z0, r = r, metrics = NULL, viable = FALSE))
  }
  
  params <- list(
    R0 = R0,
    p = c(row$p_12, row$p_21),
    omega = c(row$omega_1, row$omega_2),
    rho = rep(1, 2),
    x0 = x0,
    z0 = z0,
    r = r
  )
  
  simulations <- simulate_sis(n_days, start_interv, row$time_intervention, N, params)
  
  metrics <- compute_metrics(simulations$asynchronous$annual_incidence$value,
                             simulations$synchronous$annual_incidence$value,
                             n_days, start_interv, row$time_intervention,
                             nb_studied_cycles)
  
  list(x0 = x0, z0 = z0, r = r, metrics = metrics,
       viable = x0[1] > 0 && x0[2] > 0)
}

#' Full metric database over a parameter design
#'
#' @param myvars Data frame, one row per parameter combination, with columns
#'   `rinv`, `R0_1`, `R0_2`, `p_12`, `p_21`, `omega_1`, `omega_2` and
#'   `time_intervention`.
#' @param ncores Number of workers.
#' @param checkpoint_every Rows between checkpoint writes to
#'   `data/derived/checkpoint_metrics_computation.csv`.
#' @param n_days,start_interv Simulation horizon and first intervention day.
#'   See the note on [AIG_AIGR_computation()] before changing them.
#' @param nb_studied_cycles On/off cycles in the metric window.
#' @return A data frame with one row per design row: the input parameters, the
#'   resolved initial conditions, the controlled reproduction numbers `RC_i`,
#'   `viable`, and every [compute_metrics()] output.
#'
#' Every design row is kept, including the ones with no endemic equilibrium, so
#' that the database describes the whole design and the exclusion stays visible
#' rather than being silently baked into the file. Those rows are flagged
#' `viable = 0`; [load_simulation_database()] is what drops them, and every
#' analysis script goes through it.
metrics_computation <- function(myvars, ncores = detect_ncores(), checkpoint_every = 2000,
                                n_days = 30000, start_interv = 365,
                                nb_studied_cycles = 3) {
  
  # Rows are addressed by column name with `$`, which only works on a data
  # frame: `matrix[i, ]` returns a plain named vector. sobolSalt() hands back a
  # matrix, so coerce rather than trust the caller.
  myvars <- as.data.frame(myvars)
  
  n <- nrow(myvars)
  N <- c(1000, 1000)
  
  cat("metrics_computation(): running", n, "rows on", ncores, "core(s)\n")
  
  col_names <- c("x0_1", "x0_2", "z0_1", "z0_2", "p_12", "p_21", "r", "rinv",
                 "time_intervention", "omega_1", "omega_2", "R0_1", "R0_2",
                 "RC_1", "RC_2", "viable", METRIC_NAMES)
  
  process_one_row <- function(i) {
    row <- myvars[i, ]
    res <- simulate_design_row(row, n_days, start_interv, nb_studied_cycles, N)
    
    # `out` is assembled with c(), which coerces everything to numeric, so the
    # flag is written as 1/0 rather than TRUE/FALSE. That matches the `failed`
    # column of AIG_AIGR_computation(); load_simulation_database() converts it
    # back, because subsetting a data frame with a numeric column would index
    # by position instead of filtering.
    out <- c(res$x0[1], res$x0[2], res$z0[1], res$z0[2], row$p_12, row$p_21,
             res$r, row$rinv, row$time_intervention, row$omega_1, row$omega_2,
             row$R0_1, row$R0_2,
             row$R0_1 * (1 - row$omega_1), row$R0_2 * (1 - row$omega_2),
             as.numeric(res$viable),
             res$metrics)
    names(out) <- col_names
    out
  }
  
  all_chunks <- list()
  t0 <- Sys.time()
  
  for (chunk_start in seq(1, n, by = checkpoint_every)) {
    chunk_idx <- chunk_start:min(chunk_start + checkpoint_every - 1, n)
    
    # mc.preschedule = FALSE dispatches one row at a time instead of splitting
    # the chunk into fixed contiguous blocks. A single slow solve near an
    # elimination threshold would otherwise stall its whole block while the
    # other cores sit idle.
    chunk_results <- mclapply(chunk_idx, process_one_row,
                              mc.cores = ncores, mc.preschedule = FALSE)
    all_chunks[[length(all_chunks) + 1]] <- do.call(rbind, chunk_results)
    
    cat("Rows", min(chunk_idx), "-", max(chunk_idx), "/", n, "done -",
        format(Sys.time() - t0), "\n")
    
    write.csv(as.data.frame(do.call(rbind, all_chunks)),
              data_path("checkpoint_metrics_computation.csv"), row.names = FALSE)
  }
  
  as.data.frame(do.call(rbind, all_chunks))
}

#' AIG and AIGR for a Sobol design
#'
#' Same loop as [metrics_computation()], but keeps only the two metrics the
#' Sobol analysis needs. Both come out of a single [compute_metrics()] call, so
#' returning them together costs no extra simulation and avoids running the
#' whole design twice.
#'
#' @param myvars Data frame, same columns as in [metrics_computation()].
#' @param ncores Number of workers.
#' @param n_days,start_interv Simulation horizon and first intervention day.
#' @param nb_studied_cycles On/off cycles in the metric window.
#' @return A data frame with columns `AIG_year_area1` and `AIGR_year_area1`, one
#'   row per design row. Rows with no endemic equilibrium get zeros; the Sobol
#'   design is filtered to exclude them beforehand.
#'
#' The defaults here, interventions starting on day 730 over a 35,000-day
#' horizon, differ from those of [metrics_computation()], day 365 over 30,000
#' days. The difference is cosmetic: the model sits exactly at the endemic
#' equilibrium until the intervention starts, so the trajectory that follows is
#' the same one shifted by a year, and [compute_metrics()] locates its window
#' from `start_interv` and `n_days`, which shifts with it. Both setups return
#' identical metrics.
AIG_AIGR_computation <- function(myvars, ncores = detect_ncores(),
                                 n_days = 35000, start_interv = 730,
                                 nb_studied_cycles = 3,
                                 checkpoint_every = 10000,
                                 checkpoint_file = data_path("checkpoint_sobol_Y.csv")) {
  
  # See the note in metrics_computation(): sobolSalt()$X is a matrix, whose
  # rows do not support `$`.
  myvars <- as.data.frame(myvars)
  
  n <- nrow(myvars)
  N <- c(1000, 1000)
  
  # Resume from a checkpoint if one is present. The Sobol design is built from
  # seeded inputs, so row i means the same thing across runs and completed rows
  # can be trusted. Delete the checkpoint to force a fresh run.
  done <- list()
  n_done <- 0L
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    previous <- as.matrix(read.csv(checkpoint_file))
    if (nrow(previous) > n) {
      stop("AIG_AIGR_computation(): the checkpoint has more rows than the design. ",
           "It belongs to a different run; delete it before restarting.")
    }
    done <- list(previous)
    n_done <- nrow(previous)
    message(sprintf("Resuming from %s: %d of %d rows already done.",
                    basename(checkpoint_file), n_done, n))
  }
  
  cat("AIG_AIGR_computation(): running", n - n_done, "remaining rows on",
      ncores, "core(s)\n")
  t0 <- Sys.time()
  
  process_one_row <- function(i) {
    if (i %% 1000 == 0) {
      cat("Row", i, "/", n, "-", format(Sys.time() - t0), "\n")
    }
    res <- simulate_design_row(myvars[i, ], n_days, start_interv, nb_studied_cycles, N,
                               skip_if_not_viable = TRUE)
    if (!res$viable) {
      return(c(AIG_year_area1 = 0, AIGR_year_area1 = 0, failed = 1))
    }
    c(AIG_year_area1 = unname(res$metrics["AIG_year_area1"]),
      AIGR_year_area1 = unname(res$metrics["AIGR_year_area1"]),
      failed = 0)
  }
  
  if (n_done < n) {
    for (chunk_start in seq(n_done + 1L, n, by = checkpoint_every)) {
      chunk_idx <- chunk_start:min(chunk_start + checkpoint_every - 1L, n)
      
      raw <- mclapply(chunk_idx, process_one_row,
                      mc.cores = ncores, mc.preschedule = FALSE)
      
      # mclapply() reports failures as try-error elements rather than raising
      # them, so a silent failure would otherwise surface much later as an
      # unrelated subscript error on the rbind result. Surface the real message.
      failed <- vapply(raw, inherits, logical(1), "try-error")
      if (any(failed)) {
        stop(sprintf("AIG_AIGR_computation(): %d of %d rows in block %d-%d failed. First error: %s",
                     sum(failed), length(chunk_idx), min(chunk_idx), max(chunk_idx),
                     conditionMessage(attr(raw[[which(failed)[1]]], "condition"))))
      }
      
      done[[length(done) + 1L]] <- do.call(rbind, raw)
      
      cat("Rows", min(chunk_idx), "-", max(chunk_idx), "/", n, "done -",
          format(Sys.time() - t0), "\n")
      
      if (!is.null(checkpoint_file)) {
        write.csv(do.call(rbind, done), checkpoint_file, row.names = FALSE)
      }
    }
  }
  
  results <- do.call(rbind, done)
  stopifnot("The assembled results do not cover the whole design." = nrow(results) == n)
  
  cat("Rows without an endemic equilibrium:", sum(results[, "failed"]), "/", n, "\n")
  
  data.frame(AIG_year_area1 = results[, "AIG_year_area1"],
             AIGR_year_area1 = results[, "AIGR_year_area1"],
             failed = results[, "failed"])
}

# ---------------------------------------------------------------------------
# Reading the database back
# ---------------------------------------------------------------------------

#' Read the simulation database, keeping only viable parameter sets
#'
#' The single entry point for every analysis script that consumes
#' `df_simulations.csv`. Reading the file directly is a mistake: see below.
#'
#' @param path Path to the CSV written by [metrics_computation()].
#' @param quiet Whether to suppress the exclusion message.
#' @return A data frame with one row per viable parameter set, row names reset.
#'
#' A parameter set is viable when both patches have a positive endemic
#' equilibrium. Where they do not, `compute_equilibrium_prevalence()` returns
#' exactly `c(0, 0)`: the system-level reproduction number is at most 1, so the
#' disease-free state is the only equilibrium. Both trajectories are then
#' identically zero, and the metrics degenerate ASYMMETRICALLY:
#'
#'   AIG  = 0 - 0     = 0,    a finite number
#'   AIGR = (0 - 0)/0 = NaN,  undefined
#'
#' Left in place, that asymmetry silently gives the AIG and the AIGR analyses
#' different samples: the AIG histogram, its quantiles and its decision tree
#' would carry a spike of zeros that the AIGR versions drop through their own
#' `is.finite()` guards. It also drags the reported share of positive AIG down,
#' since those rows are zero rather than positive. Excluding them here, once,
#' keeps every figure and table on the same N, which is what the paper claims.
#'
#' The exclusion is not a numerical workaround: a parameter set with no
#' transmission carries no information about the consequences of intervention
#' asynchrony, because there is no epidemic to desynchronise.
#'
#' Databases written before `viable` existed are handled by recomputing the flag
#' from the equilibrium prevalences already stored in the file, so an existing
#' `df_simulations.csv` does not have to be regenerated.
load_simulation_database <- function(path = data_path("df_simulations.csv"),
                                     quiet = FALSE) {
  
  df <- read.csv(path)
  
  if (!"viable" %in% names(df)) {
    if (!all(c("x0_1", "x0_2") %in% names(df))) {
      stop("load_simulation_database(): ", basename(path), " has neither a ",
           "`viable` column nor the `x0_1`/`x0_2` columns needed to rebuild it.")
    }
    df$viable <- df$x0_1 > 0 & df$x0_2 > 0
  }
  
  # Written as 1/0 by metrics_computation(), as TRUE/FALSE when rebuilt above.
  # Without this, `df[df$viable, ]` on a numeric column would index by position.
  viable <- as.logical(df$viable)
  
  n_total <- nrow(df)
  df <- df[viable, , drop = FALSE]
  rownames(df) <- NULL
  
  if (!quiet) {
    message(sprintf(paste0("load_simulation_database(): %d of %d parameter sets ",
                           "excluded (no endemic equilibrium in at least one ",
                           "patch); %d retained."),
                    n_total - nrow(df), n_total, nrow(df)))
  }
  
  # After the exclusion the two metrics must have the same support. If this
  # fails, some viable parameter set produced no synchronous case inside the
  # metric window, which the analyses below are not written to handle.
  if (any(!is.finite(df$AIGR_area1))) {
    warning(sprintf(paste0("load_simulation_database(): %d viable parameter ",
                           "set(s) still have a non-finite AIGR_area1. The AIG ",
                           "and AIGR analyses will not be run on the same rows."),
                    sum(!is.finite(df$AIGR_area1))))
  }
  
  df
}