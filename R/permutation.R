#' Permutation test
#'
#' Runs a marginal permutation test using a Possion GLM test statistic.
#'
#' @inherit abstract_interface
#' @param n_rep number of permutation replicates
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
    beta_star <- compute_test_stat(df, test_stat)

    # compute permuted test statistics
    beta_null <- replicate(n = n_rep, expr = {
      df$pert_indicator <- sample(df$pert_indicator)
      compute_test_stat(df, test_stat)
    })

    # compute the p-value via call to fit skew t (in sceptre)
    p_val <- sceptre:::fit_skew_t(beta_null, beta_star, "both")$out_p

    if (return_permuted_test_stats) {
      resamples_df <- as.data.frame(matrix(beta_null, nrow = 1, ncol = length(beta_null)))
      colnames(resamples_df) <- paste0("resample_", seq(1, length(beta_null)))
      out <- cbind(data.frame(p_value = p_val, test_stat = beta_star), resamples_df)
    } else {
      out <- p_val
    }
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
    test_stat <- wilcox.test(x = x, y = y, exact = FALSE)$statistic[[1]]
  }
  return(test_stat)
}


#' Abstract two sample test
#'
#' An abstract a two-sample test. Pass function `two_sample_test` to carry out a given two-sample test (e.g., a t-test or a permutation test).
#'
#' @inherit abstract_interface
#' @param two_sample_test a two-sample test; should take as arguments (i) vector of expressions of target cells, (ii) vector of expressions of control cells, (iii) the indices of cells receiving the targeting grna, and (iv) the indices of the cells receiving the NT grnas.
#' @export
#' @examples
#' \dontrun{
#' two_sample_test <- function(target_cells, control_cells, response_id,
#' target_cell_indices, control_cell_indices) {
#' t.test(target_cells, control_cells)$p.value
#' }
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' grna_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' response_grna_group_pairs <-
#'  expand.grid(grna_group = c("CCNE2-TSS", "HS2-enh"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 50))
#' abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test)
#' }
abstract_two_sample_test <- function(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress, cbind_res = FALSE) {
  set.seed(4)
  # get grna assignments and target assignments; obtain indices of NT cells
  grna_targets <- get_target_assignments_via_max_op(grna_odm)
  control_cell_indices <- which(grna_targets == "non-targeting")

  # loop through the pairs, calculating a p-value for each
  res <- apply(X = response_grna_group_pairs, MARGIN = 1, FUN = function(r) {
    grna_group <- as.character(r[["grna_group"]])
    target_cell_indices <- grna_targets == grna_group
    response_id <- as.character(r[["response_id"]])
    if (progress) print(paste0("Analyzing ", response_id, " and ", grna_group))
    # get the target and control cells
    target_cells <- response_odm[[response_id, target_cell_indices]] |> as.numeric()
    control_cells <- response_odm[[response_id, control_cell_indices]] |> as.numeric()
    two_sample_test(target_cells, control_cells, target_cell_indices, control_cell_indices)
  }, simplify = FALSE)
  if (cbind_res) {
    to_attach <- dplyr::bind_rows(res)
    response_grna_group_pairs <- cbind(response_grna_group_pairs, to_attach)
  } else {
    response_grna_group_pairs$p_value <- unlist(res)
  }
  return(response_grna_group_pairs)
}
