#' Weissman method
#' Implements the Weissman method.
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
weissman <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {

  # convert to CellPopulation format
  cell_pop <- odm_to_cell_pop(response_odm, gRNA_odm)

  # normalize data
  cell_pop <- normalize_to_gemgroup_control(pop = cell_pop,
                                            control_cells = 'guide_target == "non-targeting"')

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


#' Converse from ODM to CellPopulation format
#'
#' @param response_odm ODM for response variable
#' @param gRNA_odm  ODM for gRNAs
#'
#' @return An object of type `CellPopulation` that the `perturbseq` library expects
odm_to_cell_pop <- function(response_odm, gRNA_odm){
  # load the data, transposing the matrices
  # TODO: Make proper indexed pandas data frame out of this
  response_mat_t <- response_odm |>
    load_whole_odm() |>
    Matrix::t() |>
    as.data.frame()

  feature_ids <- ondisc::get_feature_ids(response_odm)
  feature_names <- ondisc::get_feature_names(response_odm)
  genes_df <- tibble::tibble(gene_id = feature_ids)
  if(!all(is.na(feature_names))){
    genes_df$gene_name <- feature_names
  }

  gRNA_mat <- load_whole_odm(gRNA_odm)
  gRNA_names <- rownames(gRNA_mat)
  guide_identities <- gRNA_names[apply(X = gRNA_mat, MARGIN = 2, which.max)]

  cells_df <- ondisc::get_cell_covariates(response_odm) |>
    tibble::rownames_to_column(var = "cell_barcode") |>
    dplyr::mutate(guide_identity = guide_identities)
  if(!("batch" %in% names(cells_df))){
    cells_df$batch <- 1
  }
  cells_df <- cells_df |>
    dplyr::left_join(gRNA_odm |>
                       ondisc::get_feature_covariates() |>
                       tibble::rownames_to_column(var = "guide_identity") |>
                       dplyr::select(guide_identity, target),
                     by = "guide_identity") |>
    dplyr::rename(UMI_count = n_umis, gem_group = batch, guide_target = target) |>
    dplyr::mutate(gem_group = as.integer(as.factor(gem_group))) |>
    dplyr::select(cell_barcode, UMI_count, gem_group, guide_identity, guide_target)

  CellPopulation(matrix = response_mat_t,
                 cell_list = cells_df,
                 gene_list = genes_df,
                 calculate_statistics = FALSE)
}
