#' Classification trees for high-AIG scenarios
#'
#' The sensitivity analysis answers "which parameters drive the AIG". This
#' answers a different, more operational question: under what explicit decision
#' rules does asynchrony produce an impact worth acting on. Simulations are
#' labelled High or Low according to a quantile of the metric, and a CART
#' recovers the parameter thresholds separating the two.

#' Parameters offered to the tree as candidate splits
AIG_TREE_PREDICTORS <- c("R0_1", "R0_2", "omega_1", "omega_2",
                         "rinv", "p_12", "p_21", "time_intervention")

#' Fit and prune a classification tree for high versus low metric values
#'
#' @param data Simulation database, one row per parameter set.
#' @param quantile_high Quantile of `response` above which a simulation is
#'   labelled High, e.g. 0.90 for the top decile.
#' @param response Column name of the metric, as a string.
#' @param predictors Candidate split variables.
#' @return A list with the pruned `tree`, the unpruned `full_tree`, the
#'   `threshold` value of `response` used for the labelling, the selected
#'   complexity parameter `cp`, and the full `cptable`.
#'
#' The two classes are strongly imbalanced by construction (a top decile is one
#' High for nine Lows), so observations are reweighted to give each class the
#' same total weight; without that the tree would score well by predicting Low
#' everywhere. The tree is grown fully, then pruned by cost-complexity using
#' 10-fold cross-validation and the 1-SE rule: among the trees whose
#' cross-validated error is within one standard error of the minimum, the
#' simplest is kept. `rpart()`'s cptable is ordered from largest to smallest cp,
#' i.e. smallest to largest tree, so the first candidate row is that simplest
#' tree.
build_aig_tree <- function(data, quantile_high, response = "AIG_year_area1",
                           predictors = AIG_TREE_PREDICTORS) {

  # AIGR is a ratio whose denominator is the synchronous case count, so it is
  # undefined wherever the synchronous scenario eliminated transmission inside
  # the metric window. Those rows carry no information about what separates high
  # from low values and would poison the quantile, the labelling and the class
  # weights, so they are dropped rather than passed to rpart as NA.
  n_before <- nrow(data)
  data <- data[is.finite(data[[response]]), , drop = FALSE]
  if (nrow(data) < n_before) {
    message(sprintf("build_aig_tree(): dropped %d of %d rows with a non-finite %s.",
                    n_before - nrow(data), n_before, response))
  }
  
  threshold <- quantile(data[[response]], quantile_high)

  data$label_AIG_tree <- factor(
    ifelse(data[[response]] >= threshold, "High", "Low"),
    levels = c("Low", "High")
  )

  n_total <- nrow(data)
  n_high <- sum(data$label_AIG_tree == "High")
  n_low <- sum(data$label_AIG_tree == "Low")

  w <- ifelse(data$label_AIG_tree == "High",
              n_total / (2 * n_high),
              n_total / (2 * n_low))

  formula <- stats::as.formula(
    paste("label_AIG_tree ~", paste(predictors, collapse = " + "))
  )

  tree <- rpart(formula, data = data, weights = w, method = "class",
                control = rpart.control(xval = 10))

  cp_table <- tree$cptable
  min_index <- which.min(cp_table[, "xerror"])
  threshold_1se <- cp_table[min_index, "xerror"] + cp_table[min_index, "xstd"]
  best_cp <- cp_table[which(cp_table[, "xerror"] <= threshold_1se)[1], "CP"]

  list(tree = prune(tree, cp = best_cp),
       full_tree = tree,
       threshold = threshold,
       cp = best_cp,
       cptable = cp_table)
}

#' Variables used as splits in a fitted tree
#'
#' @param tree A pruned `rpart` object.
#' @return Character vector of split variables, in node order, with repeats.
tree_split_variables <- function(tree) {
  vars <- tree$frame$var
  as.character(vars[vars != "<leaf>"])
}

#' Number of terminal nodes of a fitted tree
#'
#' @param tree A pruned `rpart` object.
#' @return An integer.
tree_n_leaves <- function(tree) {
  sum(tree$frame$var == "<leaf>")
}

#' Variable and threshold of a tree's first split
#'
#' @param tree A pruned `rpart` object.
#' @return A one-row data frame with `variable` and `threshold`, both `NA` if
#'   the tree was pruned back to a single leaf.
#'
#' `rpart` stores splits in node order, primary split first at each internal
#' node, so row 1 of `$splits` is the root's primary split and its `index`
#' column holds the cut point.
tree_first_split <- function(tree) {

  if (tree_n_leaves(tree) <= 1 || is.null(tree$splits)) {
    return(data.frame(variable = NA_character_, threshold = NA_real_))
  }

  data.frame(variable = as.character(tree$frame$var[1]),
             threshold = unname(tree$splits[1, "index"]))
}

#' Refit the tree on bootstrap resamples
#'
#' @param data Simulation database.
#' @param B Number of resamples.
#' @param quantile_high,response Passed to [build_aig_tree()]. The High/Low
#'   threshold is recomputed within each resample, so the labelling is part of
#'   what is being resampled.
#' @param verbose Whether to report progress.
#' @return A list of [build_aig_tree()] results, one per resample.
bootstrap_aig_trees <- function(data, B = 500, quantile_high = 0.90,
                                response = "AIG_year_area1", verbose = TRUE) {
  lapply(seq_len(B), function(b) {
    if (verbose && b %% 100 == 0) cat("Bootstrap tree", b, "/", B, "\n")
    idx <- sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE)
    build_aig_tree(data[idx, ], quantile_high = quantile_high, response = response)
  })
}

#' Summarise a set of bootstrap trees
#'
#' @param boot_trees Output of [bootstrap_aig_trees()].
#' @return A list with
#'   \describe{
#'     \item{variable_frequency}{Percentage of trees selecting each variable at
#'       least once, counted once per tree.}
#'     \item{first_split}{One row per tree: its first-split `variable` and
#'       `threshold`.}
#'     \item{first_split_frequency}{Percentage of trees splitting first on each
#'       variable.}
#'     \item{first_split_thresholds}{Per first-split variable, the frequency and
#'       the spread of the threshold across trees.}
#'     \item{n_leaves}{Terminal-node count of each tree.}
#'   }
summarise_bootstrap_trees <- function(boot_trees) {

  B <- length(boot_trees)

  variable_frequency <- boot_trees %>%
    lapply(function(x) unique(tree_split_variables(x$tree))) %>%
    unlist() %>%
    table()
  variable_frequency <- sort(100 * variable_frequency / B, decreasing = TRUE)

  first_split <- bind_rows(lapply(boot_trees, function(x) tree_first_split(x$tree)))

  first_split_frequency <- sort(100 * table(first_split$variable) / B,
                                decreasing = TRUE)

  first_split_thresholds <- first_split %>%
    filter(!is.na(variable), !is.na(threshold)) %>%
    group_by(variable) %>%
    summarise(n = n(),
              frequency = 100 * n() / B,
              threshold_median = median(threshold),
              threshold_mean = mean(threshold),
              threshold_sd = sd(threshold),
              threshold_Q025 = quantile(threshold, 0.025),
              threshold_Q25 = quantile(threshold, 0.25),
              threshold_Q75 = quantile(threshold, 0.75),
              threshold_Q975 = quantile(threshold, 0.975),
              .groups = "drop") %>%
    arrange(desc(frequency))

  list(variable_frequency = variable_frequency,
       first_split = first_split,
       first_split_frequency = first_split_frequency,
       first_split_thresholds = first_split_thresholds,
       n_leaves = sapply(boot_trees, function(x) tree_n_leaves(x$tree)))
}

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

#' Draw a fitted tree
#'
#' @param fit A [build_aig_tree()] result.
#' @param title Plot title.
#' @return Called for its side effect; returns `NULL` invisibly.
#'
#' Split variables are renamed on a COPY of the fitted object: `rpart.plot()`
#' reads node labels from `$frame$var`, and overwriting that column in place
#' would leave the tree unusable by [tree_split_variables()] afterwards.
plot_aig_tree <- function(fit, title) {

  tree <- fit$tree
  labels <- c(PARAM_LABELS_PLAIN, "<leaf>" = "<leaf>")
  renamed <- labels[as.character(tree$frame$var)]
  # Any variable absent from the label table keeps its raw name rather than
  # becoming NA and blanking the node.
  renamed[is.na(renamed)] <- as.character(tree$frame$var)[is.na(renamed)]
  tree$frame$var <- unname(renamed)

  rpart.plot(tree, extra = 104, main = title, cex = 0.9, tweak=1.2,
             box.palette = colorRampPalette(c("gold", "green3"))(4))

  invisible(NULL)
}

#' Bar chart of a bootstrap frequency table
#'
#' @param frequency A named percentage vector.
#' @param x_label,title Axis label and title.
#' @return A ggplot object.
plot_bootstrap_frequency <- function(frequency, x_label, title) {

  df <- data.frame(variable = names(frequency),
                   frequency = as.numeric(frequency))

  ggplot(df, aes(x = frequency, y = reorder(variable, frequency))) +
    geom_col(fill = "steelblue") +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
    scale_y_discrete(labels = PARAM_LABELS) +
    labs(x = x_label, y = NULL, title = title) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold"),
          strip.background = element_rect(fill = "white", color = "white"))
}

#' Bootstrap distribution of the first-split threshold, per variable
#'
#' @param first_split The `first_split` element of
#'   [summarise_bootstrap_trees()].
#' @return A ggplot object.
plot_first_split_thresholds <- function(first_split) {

  first_split %>%
    filter(!is.na(variable), !is.na(threshold)) %>%
    ggplot(aes(x = threshold)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    facet_wrap(~variable, scales = "free",
               labeller = labeller(variable = PARAM_LABELS_PLAIN)) +
    labs(x = "First-split threshold", y = "Number of bootstrap trees",
         title = "Bootstrap distribution of first-split thresholds") +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold"))
}

#' Combined bootstrap stability figure
#'
#' @param summary Output of [summarise_bootstrap_trees()].
#' @return A patchwork object.
plot_bootstrap_stability <- function(summary) {

  p_variable <- plot_bootstrap_frequency(
    summary$variable_frequency,
    "Percentage of bootstrap trees in which variable was selected (%)",
    "CART variable-selection stability"
  )
  p_first <- plot_bootstrap_frequency(
    summary$first_split_frequency,
    "Frequency as first split (%)",
    "Stability of the first CART split"
  )

  (p_variable | p_first) / free(plot_first_split_thresholds(summary$first_split)) +
    plot_layout(heights = c(1, 1.2))
}

#' Bootstrap distribution of tree size
#'
#' @param n_leaves The `n_leaves` element of [summarise_bootstrap_trees()].
#' @return A ggplot object.
plot_tree_size_distribution <- function(n_leaves) {

  ggplot(data.frame(n_leaves = n_leaves), aes(x = n_leaves)) +
    geom_histogram(binwidth = 1, boundary = 0, fill = "steelblue", color = "white") +
    labs(x = "Number of terminal nodes", y = "Number of bootstrap trees",
         title = "Bootstrap distribution of tree size") +
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold"))
}
