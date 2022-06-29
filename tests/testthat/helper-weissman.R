#' cell_pop_to_odm
#'
#' @param cell_pop CellPopulation object
#'
#' @return A list of length two, containing the gene_odm and grna_odm
cell_pop_to_odm <- function(cell_pop){
  # temporary directory to store ODMs
  tmpdir <- tempdir()

  # create gene ODM
  cat(sprintf("Creating gene ODM...\n"))
  expression_matrix <- cell_pop$matrix
  gene_odm <- ondisc::create_ondisc_matrix_from_R_matrix(
    r_matrix = expression_matrix |> Matrix::t(),
    barcodes = rownames(expression_matrix),
    features_df = data.frame(gene_id = colnames(expression_matrix)),
    odm_fp = paste0(tmpdir, "/gene_matrix.odm"),
    metadata_fp = paste0(tmpdir, "/gene_metadata.rds")
  )

  # create gRNA ODM
  cat(sprintf("Creating gRNA ODM...\n"))
  cell_info <- cell_pop$cells |>
    dplyr::rowwise() |>
    dplyr::mutate(guide_identity = sprintf("%s_%s", experiment, guide_identity)) |>
    dplyr::ungroup()
  grna_mat <- Matrix::sparseMatrix(i = cell_info$guide_identity |>
                                     as.factor() |>
                                     as.integer(),
                                   j = 1:nrow(cell_pop$cells),
                                   dimnames = list(cell_info |>
                                                     dplyr::pull(guide_identity) |>
                                                     as.factor() |>
                                                     levels(),
                                                   cell_pop$cells |> rownames()))

  grna_odm <- ondisc::create_ondisc_matrix_from_R_matrix(
    r_matrix = as.matrix(grna_mat),
    barcodes = colnames(grna_mat),
    features_df = data.frame(grna_id = rownames(grna_mat)),
    odm_fp = paste0(tmpdir, "/grna_matrix.odm"),
    metadata_fp = paste0(tmpdir, "/grna_metadata.rds")
  )

  identity_to_perturbation <- cell_info |>
    dplyr::select(guide_identity, long_experiment) |>
    unique()

  id_to_pt_vec <- identity_to_perturbation$long_experiment
  names(id_to_pt_vec) <- identity_to_perturbation$guide_identity

  grna_odm <- grna_odm |>
    ondisc::mutate_feature_covariates(target = id_to_pt_vec[grna_odm |> ondisc::get_feature_ids()]) |>
    ondisc::mutate_feature_covariates(target = ifelse(target == "DMSO_control", "non-targeting", target))

  list(gene_odm = gene_odm, grna_odm = grna_odm)
}
