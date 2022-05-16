#' Load whole odm
#'
#' Loads data from disk into memory by cell.
#'
#' @param odm an ondisc_matrix object
#'
#' @return an in-memory matrix
#' @export
load_whole_odm <- function(odm) {
  x <- odm@ondisc_matrix
  x_dim <- dim(x)
  index_on_cell <- TRUE
  subset_vector <- ondisc:::get_subset_vector(x, index_on_cell)
  if (identical(subset_vector, NA_integer_)) {
    subset_vector <- seq(1, if (index_on_cell) x_dim[2] else x_dim[1])
  }
  out <- ondisc:::return_spMatrix_from_index(x@h5_file, subset_vector,
                                             index_on_cell, x@logical_mat, x@underlying_dimension)
  second_subset <- ondisc:::get_subset_vector(x, !index_on_cell)
  if (!identical(second_subset, NA_integer_)) {
    out <- if (index_on_cell)
      out[second_subset, , drop = FALSE]
    else out[, second_subset, drop = FALSE]
  }
  row.names(out) <- ondisc::get_feature_ids(odm)
  colnames(out) <- ondisc::get_cell_barcodes(odm)
  return(out)
}
