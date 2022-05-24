#' Liscovitch method
#'
#' Implements the differential expression method of Liscovitch
#' @inherit abstract_interface
#'
#' @export
#' @examples
#' response_odm <- load_dataset_modality("liscovitch/experiment_small/chromatin")
#' gRNA_odm <- load_dataset_modality("liscovitch/experiment_small/grna")
#' response_gRNA_group_pairs <- expand.grid(response_id = response_odm |> ondisc::get_feature_ids(),
#' grna_group = (gRNA_odm |> ondisc::get_feature_ids())[1:2])
#' res <- liscovitch_method(response_odm, gRNA_odm, response_gRNA_group_pairs)
liscovitch_method <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  chrom_data <- as.matrix(load_whole_odm(response_odm))
  grna_data <- load_whole_odm(gRNA_odm)
  cell_lib_sizes <- response_odm |>
    ondisc::get_cell_covariates() |>
    dplyr::pull(n_fragments) # if chrom, this strategy; else, sum to get lib size.
  # divide by the number of fragments per cell
  chrom_norm <- t(t(chrom_data)/cell_lib_sizes)
  # convert into (row-wise) z-scores
  z_score_mat <- apply(X = chrom_norm, MARGIN = 1, FUN = function(curr_row) scale(curr_row)[,1]) |> t()
  # determine which of the gRNAs are negative controls
  neg_control_gRNAs <- row.names(dplyr::filter(gRNA_odm |> ondisc::get_feature_covariates(),
                                               target_type == "non-targeting"))
  # assign gRNA IDs to cells via a max operation
  grna_assignments <- apply(X = grna_data, MARGIN = 2, FUN = function(col) names(which.max(col)))
  # obtain the indexes of the negative control grnas
  neg_control_idxs <- grna_assignments %in% neg_control_gRNAs
  # cycle through response-gRNA table
  p_vals <- apply(X = response_gRNA_group_pairs, MARGIN = 1, FUN = function(r) {
    grna <- as.character(r[["grna_group"]])
    response <- as.character(r[["response_id"]])
    target_cells <- z_score_mat[response, grna == grna_assignments]
    control_cells <- z_score_mat[response, neg_control_idxs]
    fit <- t.test(target_cells, control_cells)
    fit$p.value
  })
  response_gRNA_group_pairs$p_value <- p_vals
  return(response_gRNA_group_pairs)
}
