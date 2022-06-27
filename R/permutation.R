#' Permutation test
#'
#' Runs a marginal permutation test using a Possion GLM test statistic.
#'
#' @inherit abstract_interface
#' @param n_rep number of permutation replicates
#'
#' @return
#' @export
#'
#' @examples
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' gRNA_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' response_gRNA_group_pairs <-
#'  expand.grid(response_id = sample(ondisc::get_feature_ids(response_odm)),
#'             gRNA_group = sample(ondisc::get_feature_ids(gRNA_odm), 2))
#' res <- permutation_test(response_odm, gRNA_odm, response_gRNA_group_pairs)
permutation_test <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, n_rep = 100) {
  # convert n_rep to integer type (if necessary)
  if (is.character(n_rep)) n_rep <- as.integer(n_rep)

  # obtain the library sizes
  cell_covariates <- response_odm |> ondisc::get_cell_covariates()
  lib_size_cov <- if ("n_umis" %in% colnames(cell_covariates)) {
    "n_umis"
  } else if ("n_fragments" %in% colnames(cell_covariates)) {
    "n_fragments"
  } else {
    stop("Neither n_umis nor n_fragments are columns of the cell covariates matrix of response_odm.")
  }
  lib_sizes <- cell_covariates[[lib_size_cov]]

  # define the permutation test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    # create the data frame to pass to GLM
    target_cell_sizes <- lib_sizes[target_cell_indices]
    control_cell_sizes <- lib_sizes[control_cell_indices]
    df <- data.frame(exp = c(target_cells, control_cells),
                     lab = c(rep(1, length(target_cells)), rep(0, length(control_cells))),
                     lg_lib_size = log(c(target_cell_sizes, control_cell_sizes)))

    # ground truth test stat
    fit_star <- glm(formula = exp ~ lab + offset(lg_lib_size), family = poisson(), data = df)
    beta_star <- coef(fit_star)[["lab"]]

    # resampled test stats
    beta_null <- replicate(n = n_rep, expr = {
      df$lab <- sample(df$lab)
      tryCatch({
        fit_null <- glm(formula = exp ~ lab + offset(lg_lib_size), family = poisson(), data = df)
        coef(fit_null)[["lab"]]
      },
      error = function(cond) NULL,
      warning = function(cond) NULL)
    }, simplify = TRUE)

    p_val <- mean(c(Inf, abs(beta_null)) >= abs(beta_star))
    return(p_val)
  }

  # run the permutation test
  res <- abstract_two_sample_test(response_odm, gRNA_odm, response_gRNA_group_pairs, two_sample_test)
  return(res)
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
#'  expand.grid(response_id = sample(ondisc::get_feature_ids(response_odm)),
#'              gRNA_group = sample(ondisc::get_feature_ids(gRNA_odm), 2))
#' abstract_two_sample_test(response_odm, gRNA_odm, response_gRNA_group_pairs, two_sample_test)
abstract_two_sample_test <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, two_sample_test) {
  # load data
  response_mat <- load_whole_odm(response_odm, FALSE)
  response_ids <- row.names(response_mat)
  gRNA_mat <- load_whole_odm(gRNA_odm, FALSE)
  # assign gRNAs to cells
  gRNA_assignments <- apply(X = gRNA_mat, MARGIN = 2, FUN = function(col) names(which.max(col)))
  # obtain the indices of the NT cells
  neg_control_gRNAs <- row.names(dplyr::filter(gRNA_odm |> ondisc::get_feature_covariates(),
                                               target_type == "non-targeting"))
  control_cell_indices <- gRNA_assignments %in% neg_control_gRNAs

  # loop through the pairs, calculating a p-value for each
  p_vals <- apply(X = response_gRNA_group_pairs, MARGIN = 1, FUN = function(r) {
    gRNA_id <- as.character(r[["gRNA_group"]])
    target_cell_indices <- gRNA_assignments == gRNA_id
    response_id <- as.character(r[["response_id"]])
    print(paste0("Testing ", response_id, " against ", gRNA_id))
    # get the target and control cells
    target_cells <- response_mat[response_id, target_cell_indices]
    control_cells <- response_mat[response_id, control_cell_indices]
    two_sample_test(target_cells, control_cells, target_cell_indices, control_cell_indices)
  })
  response_gRNA_group_pairs$p_value <- p_vals
  return(response_gRNA_group_pairs)
}
