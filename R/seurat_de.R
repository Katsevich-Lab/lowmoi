#' Seurat DE
#'
#' Implements the Seurat differential expression method.
#' @inherit abstract_interface
#'
#' @export
seurat_de <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  # load response data
  response_mat <- load_whole_odm(response_odm)

  # get gRNA assignments and target assignments
  gRNA_targets <- get_target_assignments_via_max_op(gRNA_odm)

  # create a cell covariate matrix with assignments
  cell_metadata <- data.frame(perturbation = gRNA_targets)
  row.names(cell_metadata) <- colnames(response_mat)

  # create seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = response_mat,
                                           assay = "RNA",
                                           meta.data = cell_metadata)
  rm(response_mat)

  # normalize data
  seurat_obj <- Seurat::NormalizeData(seurat_obj)

  # set the "Idents" to perturbation
  Seurat::Idents(seurat_obj) <- "perturbation"

  # test for differential expression by looping over gRNA groups
  unique_gRNA_groups <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(X = unique_gRNA_groups, FUN = function(curr_gRNA_group) {
    curr_response_gRNA_group_pairs <- dplyr::filter(response_gRNA_group_pairs, gRNA_group == curr_gRNA_group)
    markers_res <- Seurat::FindMarkers(object = seurat_obj[curr_response_gRNA_group_pairs$response_id,],
                                       ident.1 = curr_gRNA_group, ident.2 = "non-targeting", only.pos = FALSE,
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
