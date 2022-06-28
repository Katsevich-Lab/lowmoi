#' Mimosca
#' Implements the MIMOSCA method.
#'
#' @inherit abstract_interface
#' @param n_rep number of permutations to use
#' @export
mimosca <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, n_rep = 1000) {
  if (is.character(n_rep))  n_rep <- as.integer(n_rep)
  # load the gene and gRNA data, transposing the matrices (as is required by mimosca)
  response_mat_t <- Matrix::t(load_whole_odm(response_odm))
  response_names <- colnames(response_mat_t)
  gRNA_mat_t <- Matrix::t(load_whole_odm(gRNA_odm))
  gRNA_names <- colnames(gRNA_mat_t)

  # if gRNA_odm is count-based (i.e., non-logical), assign gRNA IDs to cells via the max operation
  if (!gRNA_odm@ondisc_matrix@logical_mat) {
    gRNA_mat_t <- apply(X = gRNA_mat_t, MARGIN = 1, FUN = function(col) {
      out <- integer(length = ncol(gRNA_mat_t))
      out[which.max(col)] <- 1L
      return(out)
    }) |> Matrix::t()
    gRNA_mat_t <- methods::as(gRNA_mat_t, "dgCMatrix") == 1
  }

  # obtain the dataset-specific confounder model matrix
  form <- response_odm@misc$mimosca_formula
  cell_cov_mat <- response_odm |> ondisc::get_cell_covariates()
  confounder_mat <- model.matrix(object = form, data = cell_cov_mat)
  confounder_mat <- as(confounder_mat, "dgCMatrix")
  design_mat <- cbind(gRNA_mat_t, confounder_mat)

  # apply mimosca, looping over unique gRNAs
  unique_gRNAs <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(unique_gRNAs, function(curr_gRNA) {
      cov_ind <- which(curr_gRNA == gRNA_names) - 1L
      p_vals <- run_mimosca(get_sparse_matrix_pieces(response_mat_t),
                            get_sparse_matrix_pieces(design_mat),
                            cov_ind,
                            n_rep)
      out_df <- data.frame(response_id = response_names,
                           gRNA_group = curr_gRNA,
                           p_value = p_vals)
  })
  ret <- dplyr::left_join(response_gRNA_group_pairs,
                          do.call(what = "rbind", args = res_list))
  return(ret)
}
