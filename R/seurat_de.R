#' Seurat DE
#'
#' Implements the Suerat differential expression method.
#' @inherit abstract_interface
#'
#' @export
#' @examples
#' \dontrun{
#' library(ondisc)
#' papalexi_offsite <- .get_config_path("LOCAL_PAPALEXI_2021_DATA_DIR")
#' papalexi_gene = c(response_odm_fp = paste0(papalexi_offsite,
#' "processed/gene/expression_matrix.odm"),
#'                  response_metadata_fp = paste0(papalexi_offsite,
#'                  "processed/gene/metadata.rds"),
#'                  gRNA_odm_fp = paste0(papalexi_offsite,
#'                  "processed/gRNA/count_matrix.odm"),
#'                  gRNA_metadata_fp = paste0(papalexi_offsite,
#'                  "processed/gRNA/metadata.rds"))
#' response_odm <- read_odm(papalexi_gene["response_odm_fp"], papalexi_gene["response_metadata_fp"])
#' gRNA_odm <- read_odm(papalexi_gene["gRNA_odm_fp"], papalexi_gene["gRNA_metadata_fp"])
#' response_gRNA_group_pairs <- data.frame(response_id = get_feature_ids(response_odm),
#' gRNA_group = "CUL3g1")
#' result <- seurat_de(response_odm, gRNA_odm, response_gRNA_group_pairs)
#' }
seurat_de <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  # Step 0: load response ODM and gRNA ODM into memory
  response_mat <- load_whole_odm(response_odm)
  gRNA_mat <- load_whole_odm(gRNA_odm)
  gRNA_feature_covariates <- gRNA_odm |> ondisc::get_feature_covariates()

  # Step 1: assign perturbation identities to cells using the procedure outline in the Methods section of Papalexi 2021
  # i. compute gRNA library sizes
  cell_lib_sizes <- Matrix::colSums(gRNA_mat)
  # ii. remove cells with a gRNA library size of less than 5
  ok_cells <- cell_lib_sizes >= 5
  if (!all(ok_cells)) {
    response_mat <- response_mat[,ok_cells]
    gRNA_mat <- gRNA_mat[,ok_cells]
  }
  # iii. find the gRNA with the highest count in each cell, and assign it to that cell (confirmed match to Papalexi)
  gRNA_assignments <- apply(X = gRNA_mat, MARGIN = 2, FUN = function(col) names(which.max(col)))
  cell_metadata <- data.frame(perturbation = gRNA_assignments) |>
    dplyr::mutate(ondisc::get_cell_covariates(response_odm) |> dplyr::select(-n_nonzero, -n_umis))

  # Step 2: Put the gene expressions into a Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = response_mat, assay = "RNA", meta.data = cell_metadata)
  rm(response_mat)

  # Step 3: Normalize the expression data (consider also scTransform)
  seurat_obj <- Seurat::NormalizeData(seurat_obj)

  # Step 4: Aggregate the negative controls (and ONLY negative controls)
  neg_control_gRNAs <- gRNA_feature_covariates |>
    dplyr::filter(target_type == "non-targeting") |>
    row.names()
  aggregated_gRNAs <- ifelse(cell_metadata$perturbation %in% neg_control_gRNAs, "NT", cell_metadata$perturbation)
  seurat_obj$aggregate_gRNA <- aggregated_gRNAs
  Seurat::Idents(seurat_obj) <- "aggregate_gRNA"

  # Step 5: Carry out the differential expression analysis gRNA group by gRNA group
  unique_gRNA_groups <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(X = unique_gRNA_groups, FUN = function(curr_gRNA) {
    curr_response_gRNA_group_pairs <- response_gRNA_group_pairs |>
      dplyr::filter(gRNA_group == curr_gRNA)
    # run FindMarkers
    markers_res <- Seurat::FindMarkers(object = seurat_obj, ident.1 = curr_gRNA,
                               ident.2 = "NT", only.pos = FALSE,
                               logfc.threshold = 0.25, test.use = "wilcox")
    data.frame(response_id = row.names(markers_res),
               gRNA_group = curr_gRNA,
               p_value = markers_res$p_val)
  })
  res <- data.table::rbindlist(l = res_list)
  return(res)
}


load_whole_odm <- function(odm) {
  out <- odm[[,seq(1, ncol(odm))]]
  row.names(out) <- ondisc::get_feature_ids(odm)
  colnames(out) <- ondisc::get_cell_barcodes(odm)
  return(out)
}
