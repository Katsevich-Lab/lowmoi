#' Permutation test
#'
#' Runs a marginal permutation test using a Possion GLM test statistic.
#'
#' @inherit abstract_interface
#' @param n_rep number of permutation replicates
#' @export
permutation_test <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, n_rep = 10000) {
  # convert n_rep to integer type (if necessary)
  if (is.character(n_rep)) n_rep <- as.integer(n_rep)

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
    beta_star <- compute_log_fold_change(df)

    # compute permuted test statistics
    beta_null <- replicate(n = n_rep, expr = {
      df$pert_indicator <- sample(df$pert_indicator)
      compute_log_fold_change(df)
    })
    p_val <- mean(c(Inf, abs(beta_null)) >= abs(beta_star))
    return(p_val)
  }

  # run the permutation test
  res <- abstract_two_sample_test(response_odm, gRNA_odm, response_gRNA_group_pairs, two_sample_test)
  return(res)
}


compute_log_fold_change <- function(df) {
  target_cell_idxs <- df$pert_indicator == 1
  control_cell_idxs <- df$pert_indicator == 0

  target_cells <- df$exp[target_cell_idxs]
  target_cell_sizes <- df$lib_size[target_cell_idxs]
  control_cells <- df$exp[control_cell_idxs]
  control_cell_sizes <- df$lib_size[control_cell_idxs]

  beta_1 <- log(sum(target_cells)/sum(target_cell_sizes)) - log(sum(control_cells)/sum(control_cell_sizes))
  return(beta_1)
}


#' Abstract two sample test
#'
#' An abstract a two-sample test. Pass function `two_sample_test` to carry out a given two-sample test (e.g., a t-test or a permutation test).
#'
#' @inherit abstract_interface
#' @param two_sample_test a two-sample test; should take as arguments (i) vector of expressions of target cells, (ii) vector of expressions of control cells, (iii) the indices of cells receiving the targeting gRNA, and (iv) the indices of the cells receiving the NT gRNAs.
#' @export
#' @examples
#' two_sample_test <- function(target_cells, control_cells, response_id, target_cell_indices, control_cell_indices) t.test(target_cells, control_cells)$p.value
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' gRNA_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' response_gRNA_group_pairs <-
#'  expand.grid(gRNA_group = c("CCNE2-TSS", "HS2-enh"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 50))
#' abstract_two_sample_test(response_odm, gRNA_odm, response_gRNA_group_pairs, two_sample_test)
abstract_two_sample_test <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, two_sample_test) {
  # load response data
  response_mat <- load_whole_odm(response_odm, FALSE)

  # get gRNA assignments and target assignments; obtain indices of NT cells
  gRNA_targets <- get_target_assignments_via_max_op(gRNA_odm)
  control_cell_indices <- which(gRNA_targets == "non-targeting")

  # loop through the pairs, calculating a p-value for each
  p_vals <- apply(X = response_gRNA_group_pairs, MARGIN = 1, FUN = function(r) {
    gRNA_group <- as.character(r[["gRNA_group"]])
    target_cell_indices <- gRNA_targets == gRNA_group
    response_id <- as.character(r[["response_id"]])
    # get the target and control cells
    target_cells <- response_mat[response_id, target_cell_indices] |> unname()
    control_cells <- response_mat[response_id, control_cell_indices] |> unname()
    two_sample_test(target_cells, control_cells, target_cell_indices, control_cell_indices)
  })
  response_gRNA_group_pairs$p_value <- p_vals
  return(response_gRNA_group_pairs)
}
