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

df <- metrics_computation(X1)
write.csv(df, data_path("df_simulations.csv"), row.names = FALSE)

# --------------------------------------------------------------------------
# Sobol design
#
# The design is built and evaluated on the FULL sample first, then restricted:
# parameter sets with no endemic equilibrium contribute a placeholder 0 rather
# than a simulated metric, and would bias the indices. Evaluating the full
# design once and reusing the values avoids re-simulating the restricted one.
# --------------------------------------------------------------------------

set.seed(SEED_SOBOL)
sobol_design_full <- sobolSalt(model = NULL, X1[myvars], X2[myvars],
                               scheme = "B", nboot = 100)
# sobolSalt() drops the column names, which AIG_AIGR_computation() needs.
colnames(sobol_design_full$X) <- myvars

Y_full <- AIG_AIGR_computation(sobol_design_full$X)
write.csv(Y_full, data_path("sobol_design_Y.csv"), row.names = FALSE)

#' Rows of a design whose parameters admit an endemic equilibrium
check_valid_prevalence <- function(Xdf) {
  sapply(seq_len(nrow(Xdf)), function(j) {
    x0 <- compute_equilibrium_prevalence(c(Xdf$R0_1[j], Xdf$R0_2[j]),
                                         c(1000, 1000), 1 / Xdf$rinv[j],
                                         Xdf$p_12[j], Xdf$p_21[j])
    x0[1] > 0 && x0[2] > 0
  })
}

valid_idx <- which(check_valid_prevalence(X1[myvars]) & check_valid_prevalence(X2[myvars]))
cat("LHS draws excluded from the Sobol indices (no endemic equilibrium):",
    nrow(X1) - length(valid_idx), "/", nrow(X1), "\n")

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
