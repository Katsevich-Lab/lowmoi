#' Permutation test
#'
#' Runs a marginal permutation test using a Possion GLM test statistic.
#'
#' @inherit abstract_interface
#' @param n_rep number of permutation replicates
#' @param progress print progress messages?
#' @param return_permuted_test_stats return the permuted test statistics (alongside the corresponding p-values)?
#' @param test_stat string indicating the test statistic to use (either "log_fold_change" or "mann_whit")
#' @export
permutation_test <- function(response_odm, grna_odm, response_grna_group_pairs, n_rep = 1000, progress = TRUE, return_permuted_test_stats = FALSE, test_stat = "log_fold_change") {
  # convert n_rep to integer type (if necessary)
  if (is.character(n_rep)) n_rep <- as.integer(n_rep)
  if (is.character(progress)) progress <- as.logical(progress)

  # obtain the library sizes
  lib_sizes <- get_library_sizes(response_odm)

  # define the permutation test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    target_cell_sizes <- lib_sizes[target_cell_indices]
    control_cell_sizes <- lib_sizes[control_cell_indices]

    # create data frame for permutation
    df <- data.frame(exp = c(target_cells, control_cells),
                     pert_indicator = c(rep(1, length(target_cells)), rep(0, length(control_cells))),
                     lib_size = (c(target_cell_sizes, control_cell_sizes)))

    # compute beta on ground truth data
    z_star <- compute_test_stat(df, test_stat)

    # compute permuted test statistics
    z_null <- replicate(n = n_rep, expr = {
      df$pert_indicator <- sample(df$pert_indicator)
      compute_test_stat(df, test_stat)
    })

    # compute the p-value via call to fit skew t (in sceptre)
    out <- run_skew_t_prepare_output(z_null, z_star, return_permuted_test_stats)
    return(out)
  }

  # run the permutation test
  cbind_res <- if (return_permuted_test_stats) TRUE else FALSE
  res <- abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress, cbind_res)
  return(res)
}


compute_test_stat <- function(df, test_stat) {
  target_cell_idxs <- df$pert_indicator == 1
  control_cell_idxs <- df$pert_indicator == 0

  target_cells <- df$exp[target_cell_idxs]
  target_cell_sizes <- df$lib_size[target_cell_idxs]
  control_cells <- df$exp[control_cell_idxs]
  control_cell_sizes <- df$lib_size[control_cell_idxs]

  if (test_stat == "log_fold_change") {
    test_stat <- log(sum(target_cells)/sum(target_cell_sizes)) - log(sum(control_cells)/sum(control_cell_sizes))
  }
  if (test_stat == "mann_whit") {
    x <- log(target_cells/target_cell_sizes + 1)
    y <- log(control_cells/control_cell_sizes + 1)
    test_stat <- stats::wilcox.test(x = x, y = y, exact = FALSE)$statistic[[1]]
  }
  return(test_stat)
}
