#' Seurat DE
#'
#' Implements the Seurat differential expression method.
#' @inherit abstract_interface
#'
#' @export
seurat_de <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  response_mat <- load_whole_odm(response_odm)
  gRNA_mat <- load_whole_odm(gRNA_odm)
  gRNA_feature_covariates <- ondisc::get_feature_covariates(gRNA_odm)
  cell_covariates <- ondisc::get_cell_covariates(response_odm)

  # get gRNA assignments and target assignments
  gRNA_assignments <- apply(X = gRNA_mat, MARGIN = 2,
                            FUN = function(col) names(which.max(col))) |> unname()
  gRNA_to_target_map <- setNames(row.names(gRNA_feature_covariates), gRNA_feature_covariates$target)
  gRNA_targets <- names(gRNA_to_target_map)[match(x = gRNA_assignments, table = gRNA_to_target_map)]

  # update the cell covariates with the assignments
  cell_metadata <- dplyr::mutate(data.frame(perturbation = gRNA_targets), cell_covariates)
  row.names(cell_metadata) <- colnames(response_mat)
  seurat_obj <- Seurat::CreateSeuratObject(counts = response_mat, assay = "RNA", meta.data = cell_metadata)
  rm(response_mat)
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  neg_control_gRNAs <- row.names(dplyr::filter(gRNA_feature_covariates, target_type == "non-targeting"))
  aggregated_gRNAs <- ifelse(gRNA_assignments %in% neg_control_gRNAs, "NT", cell_metadata$perturbation)
  seurat_obj$aggregate_gRNA <- aggregated_gRNAs
  Seurat::Idents(seurat_obj) <- "aggregate_gRNA"
  unique_gRNA_groups <- as.character(unique(response_gRNA_group_pairs$gRNA_group))

  res_list <- lapply(X = unique_gRNA_groups, FUN = function(curr_gRNA_group) {
    curr_response_gRNA_group_pairs <- dplyr::filter(response_gRNA_group_pairs, gRNA_group == curr_gRNA_group)
    markers_res <- Seurat::FindMarkers(object = seurat_obj[curr_response_gRNA_group_pairs$response_id,],
                                       ident.1 = curr_gRNA_group, ident.2 = "NT", only.pos = FALSE,
                                       logfc.threshold = 0.0, test.use = "wilcox", min.pct = 0.0)

    if (nrow(markers_res) == 0) {
      ret <- as.data.frame(matrix(nrow = 0, ncol = 3))
      colnames(ret) <- c("gRNA_group", "p_value", "response_id")
    } else {
      ret <- data.frame(response_id = row.names(markers_res), gRNA_group = curr_gRNA_group,
                        p_value = markers_res$p_val)
    }
    return(ret)
  })
  res <- dplyr::bind_rows(l = res_list)
  return(res)
}
