#' Mimosca
#' Implements the MIMOSCA method.
#'
#' @inherit abstract_interface
#' @param n_rep number of permutations to use
#' @export
mimosca <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, n_rep = 50) {
  if (is.character(n_rep)) n_rep <- as.integer(n_rep)
  # load the gene and gRNA data, transposing the matrices (as is required by mimosca)
  response_mat_t <- Matrix::t(load_whole_odm(response_odm))
  response_names <- colnames(response_mat_t)

  # load the GROUPED gRNA matrix into memory; transpose and convert to dgC
  all_targets <- gRNA_odm |> ondisc::get_feature_covariates() |> dplyr::pull(target) |> unique()
  gRNA_mat_t <- ondisc::load_thresholded_and_grouped_gRNA(covariate_odm = gRNA_odm,
                                                          gRNA_group = all_targets,
                                                          gRNA_group_name = "target") |> t()
  gRNA_mat_t <- methods::as(gRNA_mat_t, "dgCMatrix")

  # obtain the dataset-specific confounder model matrix, and construct the design matrix
  form <- response_odm@misc$mimosca_formula
  cell_cov_mat <- response_odm |> ondisc::get_cell_covariates()
  confounder_mat <- stats::model.matrix(object = form, data = cell_cov_mat)
  confounder_mat <- methods::as(confounder_mat, "dgCMatrix")
  design_mat <- cbind(gRNA_mat_t, confounder_mat)

  # apply mimosca, looping over unique gRNA groups
  unique_gRNA_groups <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(unique_gRNA_groups, function(curr_gRNA_group) {
      cov_ind <- which(curr_gRNA_group ==  colnames(design_mat)) - 1L
      p_vals <- run_mimosca(get_sparse_matrix_pieces(response_mat_t),
                            get_sparse_matrix_pieces(design_mat),
                            cov_ind,
                            n_rep)
      out_df <- data.frame(response_id = response_names,
                           gRNA_group = curr_gRNA_group,
                           p_value = p_vals)
  })
  ret <- dplyr::left_join(response_gRNA_group_pairs,
                          do.call(what = "rbind", args = res_list))
  return(ret)
}
