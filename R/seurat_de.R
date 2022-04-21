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
  response_mat <- load_whole_odm(response_odm)
  gRNA_mat <- load_whole_odm(gRNA_odm)
  gRNA_feature_covariates <- ondisc::get_feature_covariates(gRNA_odm)
  cell_lib_sizes <- Matrix::colSums(gRNA_mat)
  ok_cells <- cell_lib_sizes >= 5
  cell_covariates <- ondisc::get_cell_covariates(response_odm)
  if (!all(ok_cells)) {
    response_mat <- response_mat[, ok_cells]
    gRNA_mat <- gRNA_mat[, ok_cells]
    cell_covariates <- cell_covariates[ok_cells, ]
  }
  gRNA_assignments <- apply(X = gRNA_mat, MARGIN = 2, FUN = function(col) names(which.max(col)))
  cell_metadata <- dplyr::mutate(data.frame(perturbation = gRNA_assignments),
                                 dplyr::select(cell_covariates, -n_nonzero, -n_umis))
  seurat_obj <- Seurat::CreateSeuratObject(counts = response_mat,
                                           assay = "RNA", meta.data = cell_metadata)
  rm(response_mat)
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  neg_control_gRNAs <- row.names(dplyr::filter(gRNA_feature_covariates,
                                               target_type == "non-targeting"))
  aggregated_gRNAs <- ifelse(cell_metadata$perturbation %in%
                               neg_control_gRNAs, "NT", cell_metadata$perturbation)
  seurat_obj$aggregate_gRNA <- aggregated_gRNAs
  Seurat::Idents(seurat_obj) <- "aggregate_gRNA"
  unique_gRNA_groups <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(X = unique_gRNA_groups, FUN = function(curr_gRNA) {
    curr_response_gRNA_group_pairs <- dplyr::filter(response_gRNA_group_pairs, gRNA_group == curr_gRNA)
    markers_res <- Seurat::FindMarkers(object = seurat_obj,
                                       ident.1 = curr_gRNA, ident.2 = "NT", only.pos = FALSE,
                                       logfc.threshold = 0.0, test.use = "wilcox")
    if (nrow(markers_res) == 0) {
      ret <- as.data.frame(matrix(nrow = 0, ncol = 3))
      colnames(ret) <- c("gRNA_group", "p_value", "response_id")
    } else {
      ret <- data.frame(response_id = row.names(markers_res), gRNA_group = curr_gRNA,
                        p_value = markers_res$p_val)
    }
    return(ret)
  })
  res <- dplyr::bind_rows(l = res_list)
  return(res)
}


load_whole_odm <- function(odm) {
  out <- odm[[,seq(1, ncol(odm))]]
  row.names(out) <- ondisc::get_feature_ids(odm)
  colnames(out) <- ondisc::get_cell_barcodes(odm)
  return(out)
}
