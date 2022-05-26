#' Mimosca
#' Implements the MIMOSCA method.
#'
#' @inherit abstract_interface
#' @export
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("frangieh/control/gene")
#' response_odm <- response_odm[1:100,]
#' gRNA_odm <- load_dataset_modality("frangieh/control/grna")
#' response_gRNA_group_pairs <- data.frame(response_id = response_odm |> ondisc::get_feature_ids(),
#' gRNA_group = "TTGCGGCCTCGATACGATAT")
#' }
mimosca <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  if (FALSE) {
  # first, load the data
  response_mat <- load_whole_odm(response_odm)
  gRNA_mat <- load_whole_odm(gRNA_odm)
  gRNAs <- unique(response_gRNA_group_pairs$gRNA_group)

  use_python("/Library/Frameworks/Python.framework/Versions/3.10/bin/python3")
  # simple example: increment an integer
  setwd("/Users/timbarry/research_code/lowmoi/")
  ret <- import("mimosca")

  # define get_pieces function
  get_pieces <- function(csc_mat) {
    list(csc_mat@x, csc_mat@i, csc_mat@p, csc_mat@Dim[1], csc_mat@Dim[2])
  }

  # load example data
  response_odm <- lowmoi::load_dataset_modality("frangieh/control/gene")
  response_odm <- response_odm[1:100,]
  gRNA_odm <- lowmoi::load_dataset_modality("frangieh/control/grna")
  my_gRNA <- "TTGCGGCCTCGATACGATAT"
  response_gRNA_group_pairs <- data.frame(response_id = response_odm |> get_feature_ids(),
                                          gRNA_group = my_gRNA)
  response_mat_t <- Matrix::t(lowmoi::load_whole_odm(response_odm))
  gRNA_mat_t <- Matrix::t(lowmoi::load_whole_odm(gRNA_odm))
  # call mimosca
  cov_ind <- which(my_gRNA == colnames(gRNA_mat_t)) - 1L
  out <- ret$run_mimosca(get_pieces(response_mat_t), get_pieces(gRNA_mat_t), cov_ind, 30)
  }
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
  # load the mimosca function
  retic <- tryCatch({
    reticulate::py_run_file(system.file("python", "mimosca.py", package = "lowmoi"))
  }, error = function(cond) {
    reticulate::py_run_file(system.file("python", "mimosca.py", package = "lowmoi"))
  })
  response_mat <- load_whole_odm(response_odm)
  out <- retic$get_dense_array(get_sparse_matrix_pieces(response_mat))
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
