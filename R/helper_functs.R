#' Load whole odm
#'
#' Loads data from disk into memory by cell.
#'
#' @param odm an ondisc_matrix object
#'
#' @return an in-memory matrix
#' @export
load_whole_odm <- function(odm) {
  x <- odm@ondisc_matrix
  x_dim <- dim(x)
  index_on_cell <- TRUE
  subset_vector <- ondisc:::get_subset_vector(x, index_on_cell)
  if (identical(subset_vector, NA_integer_)) {
    subset_vector <- seq(1, if (index_on_cell) x_dim[2] else x_dim[1])
  }
  out <- ondisc:::return_spMatrix_from_index(x@h5_file, subset_vector,
                                             index_on_cell, x@logical_mat, x@underlying_dimension)
  second_subset <- ondisc:::get_subset_vector(x, !index_on_cell)
  if (!identical(second_subset, NA_integer_)) {
    out <- if (index_on_cell)
      out[second_subset, , drop = FALSE]
    else out[, second_subset, drop = FALSE]
  }
  row.names(out) <- ondisc::get_feature_ids(odm)
  colnames(out) <- ondisc::get_cell_barcodes(odm)
  return(out)
}



#' Load dataset modality
#'
#' @param data_fp (relative) file path to dataset modality
#' @param offsite_dir the SCETPRE2 offsite directory
#'
#' @return the (QC'ed) ODM
#' @export
#'
#' @examples
#' schraivogel_enhancer_screen_chr11_gene <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' schraivogel_enhancer_screen_chr11_grna <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna")
load_dataset_modality <- function(data_fp, offsite_dir = .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")) {
  modality_dir <- paste0(offsite_dir, "data/", data_fp)
  odm <- ondisc::read_odm(odm_fp = paste0(modality_dir, "/matrix.odm"),
                          metadata_fp = paste0(modality_dir, "/metadata_qc.rds"))
 return(odm)
}


#' Get data method RAM matrix from small result
#'
#' @param res the result data frame
#' @param p_increase multiplicative factor by which to increase RAM
#'
#' @return NULL (prints matrix to console)
#' @export
get_data_method_ram_matrix_from_small_result <- function(res, p_increase = 1.08) {
  summary <- res |>
    dplyr::group_by(dataset, method) |>
    dplyr::summarize(max_ram = 1.2 * max_ram[1]) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = "method", values_from = "max_ram") |>
    dplyr::arrange(dataset)

  ram_df <- summary |> dplyr::select(-dataset)
  ram_df_boost <- (p_increase * ram_df) |> ceiling()
  ram_df_boost <- ram_df_boost[, sort(colnames(ram_df_boost))]

  out_str <- "["
  for (i in 1:nrow(ram_df_boost)) {
    out_str <- paste0(out_str, "[", paste0(ram_df_boost[i,], collapse = ", "), "]", if (i == nrow(ram_df_boost)) "]" else ",", "\n")
  }
  cat(out_str)
}

