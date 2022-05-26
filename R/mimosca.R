#' Mimosca
#' Implements the MIMOSCA method.
#'
#' @inherit abstract_interface
#' @export
#' @examples
#' response_odm <- load_dataset_modality("frangieh/control/gene")
#' response_odm <- response_odm[1:100,]
#' gRNA_odm <- load_dataset_modality("frangieh/control/grna")
#' response_gRNA_group_pairs <- data.frame(response_id = response_odm |> get_feature_ids(),
#' gRNA_group = "TTGCGGCCTCGATACGATAT")
mimosca <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  # first, load the data
  response_mat <- load_whole_odm(response_odm)
  gRNA_mat <- load_whole_odm(gRNA_odm)
  gRNAs <- unique(response_gRNA_group_pairs$gRNA_group)

  # next, pass the data to Python
  if (FALSE) {
    library(ondisc)
    library(reticulate)
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
