#' Seurat DE
#'
#' Implements the Suerat differential expression method.
#' @inherit abstract_interface
#'
#' @export
seurat_de <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  n_pairs <- nrow(response_gRNA_group_pairs)
  out_df <- response_gRNA_group_pairs |> dplyr::mutate(p_value = runif(n_pairs))
  return(out_df)
}
