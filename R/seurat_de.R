seurat_de_workhorse <- function(response_odm, grna_odm, response_grna_group_pairs, method, normalize) {
  # load response data
  response_mat <- load_whole_odm(response_odm)

  # get grna assignments and target assignments
  grna_targets <- get_target_assignments_via_max_op(grna_odm)

  # create a cell covariate matrix with assignments
  cell_metadata <- data.frame(perturbation = grna_targets)
  row.names(cell_metadata) <- colnames(response_mat)

  # create seurat object
  seurat_obj <- Seurat::CreateSeuratObject(counts = response_mat,
                                           assay = "RNA",
                                           meta.data = cell_metadata)
  if (normalize) seurat_obj <- Seurat::NormalizeData(seurat_obj)
  rm(response_mat)

  # set the "Idents" to perturbation
  Seurat::Idents(seurat_obj) <- "perturbation"

  # test for differential expression by looping over grna groups
  unique_grna_groups <- as.character(unique(response_grna_group_pairs$grna_group))
  res_list <- lapply(X = unique_grna_groups, FUN = function(curr_grna_group) {
    curr_response_grna_group_pairs <- dplyr::filter(response_grna_group_pairs, grna_group == curr_grna_group)
    markers_res <- suppressWarnings(Seurat::FindMarkers(object = seurat_obj[curr_response_grna_group_pairs$response_id,],
                                    ident.1 = curr_grna_group,
                                    ident.2 = "non-targeting", only.pos = FALSE,
                                    logfc.threshold = 0.0, test.use = method,
                                    min.pct = 0.0, min.cells.feature = 0,
                                    min.cells.group = 0))

    if (nrow(markers_res) == 0) {
      ret <- as.data.frame(matrix(nrow = 0, ncol = 3))
      colnames(ret) <- c("grna_group", "p_value", "response_id")
    } else {
      ret <- data.frame(response_id = row.names(markers_res), grna_group = curr_grna_group,
                        p_value = markers_res$p_val)
    }
    return(ret)
  })
  res <- dplyr::bind_rows(l = res_list)
  return(res)
}


#' Seurat DE
#'
#' Implements the NB regression-based Seurat differential expression method.
#' @inherit abstract_interface
#'
#' @export
seurat_de_nb <- function(response_odm, grna_odm, response_grna_group_pairs) {
  seurat_de_workhorse(response_odm, grna_odm, response_grna_group_pairs, "negbinom", FALSE)
}


#' Seurat DE
#'
#' Implements the Seurat differential expression method. The (default) wilcox test is used,
#' @inherit abstract_interface
#'
#' @export
seurat_de <- function(response_odm, grna_odm, response_grna_group_pairs) {
  seurat_de_workhorse(response_odm, grna_odm, response_grna_group_pairs, "wilcox", TRUE)
}
