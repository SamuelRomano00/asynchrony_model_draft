# Uncertainty and sensitivity analysis: designs and simulations
#
# The expensive step of the project: 10,000 Latin hypercube parameter sets for
# the simulation database, plus a Sobol design of the same size. Expect hours
# on a cluster rather than minutes on a laptop; see inst/hpc/ for a submission
# template. Everything downstream reads the CSVs written here, so 05 and 06 can
# be re-run freely without touching this script.
#
# Outputs, all in data/:
#   LHS1.csv, LHS2.csv          the two independent Latin hypercube designs
#   df_simulations.csv          parameters and metrics for every design point
#   sobol_design_Y.csv          AIG and AIGR over the full Sobol design
#   sobol_AIG.csv, sobol_AIGR.csv, sobol_AIGR_rank.csv        first-order and total indices
#   sobol_AIG_S2.csv, sobol_AIGR_S2.csv, sobol_AIGR_rank_S2.csv  second-order
#                               indices; sobol_AIG_S2.csv is Supplementary Table S2

source(here::here("R", "setup.R"))

# Seeds are fixed so the designs, and therefore every CSV below, are
# reproducible. Changing them changes every downstream figure.
SEED_LHS1 <- 20250101
SEED_LHS2 <- 20250102
SEED_SOBOL <- 20250103

design_size <- 10000

myvars <- c("rinv", "time_intervention", "R0_1", "R0_2",
            "omega_1", "omega_2", "p_12", "p_21")

# --------------------------------------------------------------------------
# Latin hypercube designs
#
# Two independent samples of the same 8 parameters. X1 alone would be enough
# for the simulation database, but sobolSalt() needs two independent samples to
# build its column-swapped designs, which is how first-order, second-order and
# total indices are estimated (Saltelli 2010, scheme B).
#
# Each column is drawn on the unit cube and rescaled to its physical range from
# Table 1 of the paper. Efficacies and mobility fractions are sampled directly
# on the paper's [0, 0.5] scale, which is also the scale the models use.
# --------------------------------------------------------------------------

build_lhs_design <- function(n, seed) {
  set.seed(seed)
  X <- data.frame(as.matrix(maximinLHS(n = n, k = 8)))
  colnames(X) <- c("R0_1", "R0_2", "rinv", "time_intervention",
                   "omega_1", "omega_2", "p_12", "p_21")
  
  X$R0_1 <- 0.9 + X$R0_1 * (2.2 - 0.9)
  X$R0_2 <- 0.9 + X$R0_2 * (2.2 - 0.9)
  X$rinv <- floor(60 + X$rinv * (201 - 60))                      # infectious period, integer days in [60, 200]
  X$time_intervention <- floor(2 * (0.5 + X$time_intervention * 5)) / 2   # steps of 0.5 in [0.5, 5]
  X$omega_1 <- X$omega_1 * 0.5
  X$omega_2 <- X$omega_2 * 0.5
  X$p_12 <- X$p_12 * 0.5
  X$p_21 <- X$p_21 * 0.5
  
  X
}

X1 <- build_lhs_design(design_size, SEED_LHS1)
X2 <- build_lhs_design(design_size, SEED_LHS2)

write.csv(X1, data_path("LHS1.csv"), row.names = FALSE)
write.csv(X2, data_path("LHS2.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# Simulation database
# --------------------------------------------------------------------------

# Reuse an existing database rather than recomputing it. The designs are seeded,
# so a df_simulations.csv produced by an earlier run of this script describes the
# same 10,000 parameter sets; recomputing it would take hours to reproduce a file
# that is already correct. Delete the file to force a fresh run.
if (file.exists(data_path("df_simulations.csv"))) {
  message("df_simulations.csv already exists, reusing it. ",
          "Delete it to recompute the simulation database.")
  df <- read.csv(data_path("df_simulations.csv"))
  stopifnot("Existing df_simulations.csv does not match the current design size." =
              nrow(df) == nrow(X1))
} else {
  df <- metrics_computation(X1)
  write.csv(df, data_path("df_simulations.csv"), row.names = FALSE)
}

# --------------------------------------------------------------------------
# Sobol design
#
# The design is built and evaluated on the FULL sample first, then restricted:
# parameter sets with no endemic equilibrium contribute a placeholder 0 rather
# than a simulated metric, and would bias the indices. Evaluating the full
# design once and reusing the values avoids re-simulating the restricted one,
# and also supplies the viability flag the restriction is based on.
# --------------------------------------------------------------------------

set.seed(SEED_SOBOL)
sobol_design_full <- sobolSalt(model = NULL, X1[myvars], X2[myvars],
                               scheme = "B", nboot = 100)
# sobolSalt() drops the column names, which AIG_AIGR_computation() needs.
colnames(sobol_design_full$X) <- myvars

Y_full <- AIG_AIGR_computation(sobol_design_full$X)
write.csv(Y_full, data_path("sobol_design_Y.csv"), row.names = FALSE)

# A draw is kept only if ALL of its 2k+2 blocks admit an endemic equilibrium,
# not just the two source samples. Each block mixes columns from X1[j] and
# X2[j], so both sources can be viable while the combination is not; filtering
# on the sources alone leaves such rows in the design, contributing an AIG of
# zero that no simulation produced. Every block at index j derives from X1[j]
# and X2[j] alone, so rebuilding the design on the surviving indices reproduces
# exactly the rows checked here and preserves the paired structure the Sobol
# estimator relies on.
#
# Viability is read off the `failed` flag rather than recomputed: it is already
# one Newton solve per design row, and AIG_AIGR_computation() has done them.
n_lhs <- nrow(X1)
n_blocks <- 2 * length(myvars) + 2
stopifnot(
  "Unexpected Sobol design size." = nrow(sobol_design_full$X) == n_lhs * n_blocks,
  "Y_full carries no `failed` flag." = "failed" %in% names(Y_full)
)

# Blocks are stacked, so a column-major fill puts block b in column b.
viable <- matrix(Y_full$failed == 0, nrow = n_lhs)
valid_idx <- which(rowSums(!viable) == 0)

n_rows_dropped <- nrow(sobol_design_full$X) - length(valid_idx) * n_blocks
cat("Design rows with no endemic equilibrium:", sum(!viable), "/",
    nrow(sobol_design_full$X),
    sprintf("(%.2f%%)\n", 100 * mean(!viable)))
cat("LHS draws excluded (some block has no endemic equilibrium):",
    n_lhs - length(valid_idx), "/", n_lhs,
    sprintf("(%.2f%%)\n", 100 * (n_lhs - length(valid_idx)) / n_lhs))
cat("Design rows removed from the Sobol estimate:", n_rows_dropped, "/",
    nrow(sobol_design_full$X),
    sprintf("(%.2f%%)\n", 100 * n_rows_dropped / nrow(sobol_design_full$X)))

set.seed(SEED_SOBOL)
sobol_design <- sobolSalt(model = NULL, X1[valid_idx, myvars], X2[valid_idx, myvars],
                          scheme = "B", nboot = 100)
colnames(sobol_design$X) <- myvars

# Map each row of the restricted design back to its already-simulated value in
# the full design. Duplicate rows are collapsed first, so the join stays
# one-to-one.
full_design_df <- as.data.frame(sobol_design_full$X)
full_design_df$orig_row <- seq_len(nrow(full_design_df))
full_design_df <- full_design_df %>%
  distinct(across(all_of(myvars)), .keep_all = TRUE)

clean_design_df <- as.data.frame(sobol_design$X)
clean_design_df$new_row <- seq_len(nrow(clean_design_df))

matched <- left_join(clean_design_df, full_design_df, by = myvars)
stopifnot(
  "Some rows of the restricted design have no match in the full design." =
    all(!is.na(matched$orig_row)),
  "The join changed the number of rows of the restricted design." =
    nrow(matched) == nrow(clean_design_df)
)

Y <- Y_full[matched$orig_row[order(matched$new_row)], ]
rownames(Y) <- NULL

# --------------------------------------------------------------------------
# Sobol indices
#
# AIGR is also computed on RANKS. Near the elimination threshold the ratio can
# blow up, and a handful of such rows dominate the variance decomposition on
# the raw scale; ranks are insensitive to that.
# --------------------------------------------------------------------------

#' Collect first-order and total indices from a `sobolSalt` object
collect_indices <- function(obj) {
  first_order <- obj$S
  first_order$X <- myvars
  first_order$index <- "first order"
  
  total <- obj$T
  total$X <- myvars
  total$index <- "total"
  
  out <- rbind(first_order, total)
  rownames(out) <- NULL
  out
}

#' Collect second-order indices from a `sobolSalt` object
collect_indices_S2 <- function(obj) {
  out <- obj$S2
  out$X <- combn(myvars, 2, FUN = function(v) paste(v, collapse = "*"))
  out$index <- "second order"
  rownames(out) <- NULL
  out
}

sobol_AIG_obj <- sobol_design
tell(sobol_AIG_obj, Y$AIG_year_area1)

sobol_AIGR_obj <- sobol_design
tell(sobol_AIGR_obj, Y$AIGR_year_area1)

sobol_AIGR_rank_obj <- sobol_design
tell(sobol_AIGR_rank_obj, rank(Y$AIGR_year_area1, ties.method = "average"))

write.csv(collect_indices(sobol_AIG_obj), data_path("sobol_AIG.csv"), row.names = FALSE)
write.csv(collect_indices(sobol_AIGR_obj), data_path("sobol_AIGR.csv"), row.names = FALSE)
write.csv(collect_indices(sobol_AIGR_rank_obj), data_path("sobol_AIGR_rank.csv"), row.names = FALSE)

write.csv(collect_indices_S2(sobol_AIG_obj), data_path("sobol_AIG_S2.csv"), row.names = FALSE)
write.csv(collect_indices_S2(sobol_AIGR_obj), data_path("sobol_AIGR_S2.csv"), row.names = FALSE)
write.csv(collect_indices_S2(sobol_AIGR_rank_obj), data_path("sobol_AIGR_rank_S2.csv"), row.names = FALSE)