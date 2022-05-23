#' Liscovitch method
#'
#' Implements the differential expression method of Liscovitch
#' @inherit abstract_interface
#'
#' @export
#' @examples
#' response_odm <- load_dataset_modality("liscovitch/experiment_small/chromatin")
#' gRNA_odm <- load_dataset_modality("liscovitch/experiment_small/grna")
liscovitch_method <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  chrom_data <- as.matrix(load_whole_odm(response_odm))
  grna_data <- load_whole_odm(gRNA_odm)
  cell_lib_size <- response_odm |>
    ondisc::get_cell_covariates() |>
    dplyr::pull(n_fragments)

}
