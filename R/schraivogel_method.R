#' Run Schraivogel's MAST.cov method
#'
#' @inherit abstract_interface
#' @param gRNA_groups_table A table specifying which gRNAs are in which groups, as in \code{sceptre}.
#' This argument is optional, and the default assumption is that each gRNA is in its own group.
#' @param gRNA_threshold A threshold for gRNA expression. This argument is optional, and defaults to 8,
#' which was Schraivogel et al's choice.
#'
#' @export
schraivogel_method <- function(response_odm,
                               gRNA_odm,
                               response_gRNA_group_pairs,
                               gRNA_groups_table = NULL,
                               gRNA_threshold = 8) {

  gene_odm <- response_odm
  grna_odm <- gRNA_odm
  gene_gRNA_group_pairs <- response_gRNA_group_pairs
  # pull the entire gRNA ODM into memory
  grna_data <- grna_odm[[1:nrow(grna_odm),1:ncol(grna_odm)]]
  rownames(grna_data) <- grna_odm |> ondisc::get_feature_ids()
  colnames(grna_data) <- grna_odm |> ondisc::get_cell_barcodes()

  # pull the entire gene expression ODM into memory
  gene_data <- gene_odm[[1:nrow(gene_odm),1:ncol(gene_odm)]]
  rownames(gene_data) <- gene_odm |> ondisc::get_feature_ids()
  colnames(gene_data) <- gene_odm |> ondisc::get_cell_barcodes()

  # threshold the gRNA matrix, unless it is already binary
  if(max(grna_data) >= 2){
    perturbation_matrix <- sceptre::threshold_gRNA_matrix(grna_data, threshold = gRNA_threshold)
  } else{
    perturbation_matrix <- grna_data
  }

  # combine the perturbation indicators within each gRNA group, unless gRNA_groups_table is not given
  if (!is.null(gRNA_groups_table)) {
    combined_perturbation_matrix <- sceptre::combine_perturbations(
      perturbation_matrix = perturbation_matrix,
      gRNA_groups_table = gRNA_groups_table
    )
  } else {
    combined_perturbation_matrix <- perturbation_matrix
  }

  # identify which of the gRNAs are negative controls
  scramble.cols <- grna_odm |>
    ondisc::get_feature_covariates() |>
    dplyr::mutate(NTC = target_type == "non-targeting") |>
    dplyr::pull(NTC) |>
    which()

  # transpose the perturbation matrix, as this is what the runSeuratTest function requires
  combined_pert_mat_t = combined_perturbation_matrix |> Matrix::t()

  # compute the number of genes expressed in each cell, which is the covariate that
  # this method adjusts for
  ngenes <- apply(gene_data > 0, 2, sum)

  # loop over distinct gRNAs present in gene_gRNA_group_pairs
  lapply(
    gene_gRNA_group_pairs |> dplyr::pull(gRNA_group) |> unique(),
    function(gRNA_group) {
      # call the workhorse function, runSeuratTest
      # Note: runSeuratTest will compute the p-values for a given gRNA against
      #       ALL GENES; we then can extract the pairs we are interested in.
      #       This is unavoidable because the MAST.cov test depends on which
      #       genes are included, and in the original paper it was run on all genes.
      runSeuratTest(
        g = gRNA_group,
        DGE = gene_data,
        pert = combined_pert_mat_t,
        covariate = ngenes,
        scrcols = scramble.cols,
        normfun = function(x) sum(x[x < stats::quantile(x, probs = 0.9)]) + 1
      )
    }
  ) |>
    # combine the results, extract gene-gRNA group pairs of interest, and format for output
    dplyr::bind_rows() |>
    dplyr::rename(response_id = gene, gRNA_group = guide, p_value = p_val) |>
    dplyr::select(response_id, gRNA_group, p_value) |>
    dplyr::semi_join(gene_gRNA_group_pairs, by = c("response_id", "gRNA_group")) |>
    tibble::as_tibble()
}
