#' Liscovitch method
#'
#' Implements the differential expression method of Liscovitch
#' @inherit abstract_interface
#'
#' @export
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("liscovitch/experiment_small/chromatin")
#' grna_odm <- load_dataset_modality("liscovitch/experiment_small/grna_assignment")
#' response_grna_group_pairs <- expand.grid(response_id = response_odm |> ondisc::get_feature_ids(),
#' grna_group = c("ARID1A", "ATRX"))
#' res <- liscovitch_method(response_odm, grna_odm, response_grna_group_pairs)
#' }
liscovitch_method <- function(response_odm, grna_odm, response_grna_group_pairs) {
  exp_data <- load_whole_odm(response_odm)
  # obtain the library sizes
  cell_covariates <- response_odm |> ondisc::get_cell_covariates()
  lib_size_col_name <- if ("n_fragments" %in% colnames(cell_covariates)) "n_fragments" else "n_umis"
  cell_lib_sizes <- cell_covariates[[lib_size_col_name]]
  # divide by the number of fragments per cell
  exp_data_norm <- Matrix::t(Matrix::t(exp_data)/cell_lib_sizes)
  rm(exp_data)

  # assign grnas to cells; obtain the indexes of the negative control grnas
  grna_targets <- get_target_assignments_via_max_op(grna_odm)
  neg_control_idxs <- which(grna_targets == "non-targeting")

  # cycle through response-grna table
  p_vals <- apply(X = response_grna_group_pairs, MARGIN = 1, FUN = function(r) {
    grna <- as.character(r[["grna_group"]])
    response <- as.character(r[["response_id"]])
    response_row <- exp_data_norm[response,]
    response_row_scale <- scale(response_row)[,1]

    target_cells <- response_row_scale[grna == grna_targets]
    control_cells <- response_row_scale[neg_control_idxs]
    fit <- stats::t.test(target_cells, control_cells)
    fit$p.value
  })
  response_grna_group_pairs$p_value <- p_vals
  return(response_grna_group_pairs)
}
