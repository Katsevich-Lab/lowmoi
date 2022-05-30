#' Mimosca
#' Implements the MIMOSCA method.
#'
#' @inherit abstract_interface
#' @export
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' gRNA_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna")
#' response_gRNA_group_pairs <- expand.grid(response_id = (response_odm |> ondisc::get_feature_ids()), gRNA_group = c("GATA1-C", "GATA1-D"))
#' result <- mimosca(response_odm, gRNA_odm, response_gRNA_group_pairs)
#' }
mimosca <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, n_rep = 1000) {
  # load the data, transposing the matrices (as is required by mimosca)
  response_mat_t <- Matrix::t(load_whole_odm(response_odm))
  response_names <- colnames(response_mat_t)
  gRNA_mat_t <- Matrix::t(load_whole_odm(gRNA_odm))
  # assign gRNA IDs to cells via the max operation
  gRNA_mat_t_bin <- apply(X = gRNA_mat_t, MARGIN = 1, FUN = function(col) {
    out <- integer(length = ncol(gRNA_mat_t))
    out[which.max(col)] <- 1L
    return(out)
  }) |> Matrix::t()
  gRNA_mat_t_bin <- methods::as(gRNA_mat_t_bin, "dgCMatrix") == 1
  gRNA_names <- colnames(gRNA_mat_t)
  # apply mimosca
  unique_gRNAs <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(unique_gRNAs, function(curr_gRNA) {
      cov_ind <- which(curr_gRNA == gRNA_names) - 1L
      p_vals <- run_mimosca(get_sparse_matrix_pieces(response_mat_t),
                            get_sparse_matrix_pieces(gRNA_mat_t),
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

#' Return dense matrix
#'
#' Returns a dense matrix representation of response_odm
#'
#' @param response_odm an ondisc object
#' @export
#'
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' dense_mat <- return_dense_mat(response_odm)
#' }
return_dense_mat <- function(response_odm) {
  response_mat <- load_whole_odm(response_odm)
  out <- get_dense_array(get_sparse_matrix_pieces(response_mat))
  return(out)
}

#' Get sparse matrix pieces
#'
#' @param csc_mat a sparse matrix in CSC format
#'
#' @return a list containing x, i, p, Dim(1), and Dim(2) (in that order)
get_sparse_matrix_pieces <- function(csc_mat) {
  list(csc_mat@x, csc_mat@i, csc_mat@p, csc_mat@Dim[1], csc_mat@Dim[2])
}
