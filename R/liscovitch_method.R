#' Liscovitch method
#'
#' Implements the differential expression method of Liscovitch
#' @inherit abstract_interface
#'
#' @export
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("liscovitch/experiment_small/chromatin")
#' gRNA_odm <- load_dataset_modality("liscovitch/experiment_small/grna_assignment")
#' response_gRNA_group_pairs <- expand.grid(response_id = response_odm |> ondisc::get_feature_ids(),
#' gRNA_group = c("ARID1A", "ATRX"))
#' res <- liscovitch_method(response_odm, gRNA_odm, response_gRNA_group_pairs)
#' }
liscovitch_method <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  chrom_data <- as.matrix(load_whole_odm(response_odm))
  # obtain the library sizes
  cell_covariates <- response_odm |> ondisc::get_cell_covariates()
  lib_size_col_name <- if ("n_fragments" %in% colnames(cell_covariates)) "n_fragments" else "n_umis"
  cell_lib_sizes <- cell_covariates[[lib_size_col_name]]
  # divide by the number of fragments per cell
  chrom_norm <- t(t(chrom_data)/cell_lib_sizes)
  # convert into (row-wise) z-scores
  z_score_mat <- apply(X = chrom_norm, MARGIN = 1, FUN = function(curr_row) scale(curr_row)[,1]) |> t()

  # assign gRNAs to cells; obtain the indexes of the negative control grnas
  gRNA_targets <- get_target_assignments_via_max_op(gRNA_odm)
  neg_control_idxs <- which(gRNA_targets == "non-targeting")

  # cycle through response-gRNA table
  p_vals <- apply(X = response_gRNA_group_pairs, MARGIN = 1, FUN = function(r) {
    grna <- as.character(r[["gRNA_group"]])
    response <- as.character(r[["response_id"]])
    target_cells <- z_score_mat[response, grna == gRNA_targets]
    control_cells <- z_score_mat[response, neg_control_idxs]
    fit <- stats::t.test(target_cells, control_cells)
    fit$p.value
  })
  response_gRNA_group_pairs$p_value <- p_vals
  return(response_gRNA_group_pairs)
}
