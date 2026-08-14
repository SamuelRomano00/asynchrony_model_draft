# Decision trees for high-AIG scenarios: Figure 4 and Figures S7, S8, S9, S10
#
# A different question from the Sobol analysis: not which parameters drive the
# AIG, but under what explicit decision rules asynchrony produces an impact of
# clinical relevance. Simulations in the top decile of AIG per year in area 1
# are labelled High, the rest Low, and a pruned CART recovers the thresholds
# separating them.
#
# Reads the database written by 04_sensitivity_simulations.R; no simulation
# happens here. The 500 bootstrap refits are the slow part, a few minutes.
#
# Outputs:
#   figures/figure4_decision_tree.png            top-decile tree
#   figures/figureS7_decision_tree_thresholds.png       the same at five thresholds
#   figures/figureS7_decision_tree_thresholds/top_*.png one tree per threshold
#   figures/figureS8_bootstrap_stability.png     variable and first-split stability
#   figures/figureS9_bootstrap_tree_size.png    tree-size distribution
#   figures/figureS10_decision_tree_AIGR.png     top-decile tree, AIGR instead
#   data/decision_tree_summary.csv       one row per threshold definition
#   data/decision_tree_bootstrap.csv     first-split threshold stability

source(here::here("R", "setup.R"))

# The trees depend on the RNG twice over: rpart's 10-fold cross-validation
# assigns folds at random, and the stability analysis resamples the database.
set.seed(42)

# Parameter sets with no endemic equilibrium are excluded before anything else.
# Their AIG is a finite zero, so reading the CSV directly would feed the AIG
# tree a block of rows labelled Low that the AIGR tree drops through its own
# is.finite() guard, and would shift the top-decile threshold. See
# load_simulation_database() in R/batch.R.
df <- load_simulation_database()

AIG_variable <- "AIG_year_area1"
n_bootstrap <- 500

# Definitions of "high AIG", as quantiles of the metric. 10% is the definition
# used in the main text; the others test how much the conclusions depend on it.
high_percentiles <- c("5%" = 0.95, "10%" = 0.90, "20%" = 0.80,
                      "25%" = 0.75, "50%" = 0.50)

# --------------------------------------------------------------------------
# One tree per definition of "high AIG"
# --------------------------------------------------------------------------

trees <- lapply(high_percentiles, function(q) {
  build_aig_tree(df, quantile_high = q, response = AIG_variable)
})

tree_summary <- data.frame(
  definition = names(trees),
  AIG_threshold = sapply(trees, function(x) x$threshold),
  n_leaves = sapply(trees, function(x) tree_n_leaves(x$tree)),
  cp = sapply(trees, function(x) x$cp),
  variables = sapply(trees, function(x) {
    paste(sort(unique(tree_split_variables(x$tree))), collapse = ", ")
  }),
  row.names = NULL
)
print(tree_summary)
write.csv(tree_summary, data_path("decision_tree_summary.csv"), row.names = FALSE)

cat("\nVariable importance, top-decile tree:\n")
print(sort(trees[["10%"]]$tree$variable.importance, decreasing = TRUE))

# --------------------------------------------------------------------------
# Figure 4: the top-decile tree
# --------------------------------------------------------------------------

png(figure_path("figure4_decision_tree.png"), width = 14, height = 8,
    units = "in", res = 300)
plot_aig_tree(trees[["10%"]], "High AIG: top 10%")
dev.off()

# --------------------------------------------------------------------------
# Figure S7: the same at five thresholds
# --------------------------------------------------------------------------

png(figure_path("figureS7_decision_tree_thresholds.png"), width = 16, height = 18,
    units = "in", res = 300)
# base graphics, so the panels are arranged with par() rather than patchwork
par(mfrow = c(3, 2))
for (i in seq_along(trees)) {
  plot_aig_tree(trees[[i]], paste0("High AIG: top ", names(high_percentiles)[i]))
}
par(mfrow = c(1, 1))
dev.off()

# Same trees again, one file each. The assembled grid is unreadable at the
# thresholds where the tree grows large, so each panel is also written at full
# size for inspection.
threshold_dir <- figure_path("figureS7_decision_tree_thresholds")
dir.create(threshold_dir, showWarnings = FALSE)

for (i in seq_along(trees)) {
  # "5%" is not usable in a filename; keep the number only.
  label <- names(high_percentiles)[i]
  slug <- sub("%$", "", label)
  
  png(file.path(threshold_dir, paste0("top_", slug, "pc.png")),
      width = 14, height = 8, units = "in", res = 300)
  plot_aig_tree(trees[[i]], paste0("High AIG: top ", label))
  dev.off()
}

# --------------------------------------------------------------------------
# Figures S8 and S9: bootstrap stability
#
# CART is known to be sensitive to sample perturbations, so the question is not
# whether individual trees are identical but whether they are built around the
# same variables with the same split thresholds.
# --------------------------------------------------------------------------

boot_trees <- bootstrap_aig_trees(df, B = n_bootstrap, quantile_high = 0.90,
                                  response = AIG_variable)
boot <- summarise_bootstrap_trees(boot_trees)

cat("\nVariable selection frequency (%):\n")
print(boot$variable_frequency)
cat("\nFirst-split frequency (%):\n")
print(boot$first_split_frequency)
cat("\nFirst-split threshold stability:\n")
print(boot$first_split_thresholds)
cat("\nTree size (terminal nodes):\n")
print(summary(boot$n_leaves))

write.csv(boot$first_split_thresholds,
          data_path("decision_tree_bootstrap.csv"), row.names = FALSE)

ggsave(figure_path("figureS8_bootstrap_stability.png"),
       plot_bootstrap_stability(boot), width = 16, height = 12, dpi = 300)

ggsave(figure_path("figureS9_bootstrap_tree_size.png"),
       plot_tree_size_distribution(boot$n_leaves), width = 10, height = 6, dpi = 300)

# --------------------------------------------------------------------------
# Figure S10: the same tree built on the AIGR
#
# AIG and AIGR rank scenarios differently near the elimination threshold, where
# a small synchronous case count inflates the ratio, so the tree is refitted on
# the relative metric as a check.
# --------------------------------------------------------------------------

tree_AIGR <- build_aig_tree(df, quantile_high = 0.90, response = "AIGR_year_area1")

cat("\nAIGR tree, threshold:", tree_AIGR$threshold,
    "- terminal nodes:", tree_n_leaves(tree_AIGR$tree), "\n")
print(sort(tree_AIGR$tree$variable.importance, decreasing = TRUE))

png(figure_path("figureS10_decision_tree_AIGR.png"), width = 14, height = 8,
    units = "in", res = 300)
plot_aig_tree(tree_AIGR, "High AIGR: top 10%")
dev.off()