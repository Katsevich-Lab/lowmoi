#' Load whole odm
#'
#' Loads data from disk into memory by cell.
#'
#' @param odm an ondisc_matrix object
#' @param csc_format load in cell-accessible, CSC format (TRUE) or feature-accessible, CSR format (FALSE)?
#'
#' @return an in-memory matrix
#' @export
load_whole_odm <- function(odm, csc_format = TRUE) {
  x <- odm@ondisc_matrix
  x_dim <- dim(x)
  index_on_cell <- csc_format
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
#' \dontrun{
#' schraivogel_enhancer_ground_truth_tapseq_gene <-
#' load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' }
load_dataset_modality <- function(data_fp, offsite_dir = .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")) {
  modality_dir <- paste0(offsite_dir, "data/", data_fp)
  odm <- ondisc::read_odm(odm_fp = paste0(modality_dir, "/matrix.odm"),
                          metadata_fp = paste0(modality_dir, "/metadata_qc.rds"))
 return(odm)
}


#' Load dataset multimodal
#'
#' @param paper_fp (relative) file path to paper modality
#' @param offsite_dir the SCETPRE2 offsite directory
#'
#' @return the (QC'ed) multimodal ODM
#' @export
#'
#' @examples
#' \dontrun{
#' schraivogel_enhancer_screen_chr11 <-
#' load_dataset_multimodal("schraivogel/enhancer_screen_chr11")
#' }
load_dataset_multimodal <- function(paper_fp, offsite_dir = .get_config_path("LOCAL_SCEPTRE2_DATA_DIR")) {
  paper_dir <- paste0(offsite_dir, "data/", paper_fp, "/")
  modalities <- list.files(paper_dir)
  modalities <- modalities[modalities != "multimodal_metadata.rds"]
  odm_fps <- sapply(X = modalities, FUN = function(modality) paste0(paper_dir, modality, "/matrix.odm")) |> unname()
  multimodal_metadata_fp <- paste0(paper_dir, "multimodal_metadata.rds")
  ondisc::read_multimodal_odm(odm_fps, multimodal_metadata_fp)
}

#' Get data method RAM matrix from small result
#'
#' @param res the result data frame
#' @param p_increase multiplicative factor by which to increase RAM
#'
#' @return NULL (prints matrix to console)
#' @export
get_data_method_ram_matrix_from_small_result <- function(res, p_increase = 1.1) {
  summary <- res |>
    dplyr::group_by(dataset, method) |>
    dplyr::summarize(max_ram = max_ram[1]) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = "method", values_from = "max_ram") |>
    dplyr::arrange(dataset)

  ram_df <- summary |> dplyr::select(-dataset)
  ram_df_boost <- (p_increase * ram_df) |> ceiling()
  ram_df_boost <- ram_df_boost[, sort(colnames(ram_df_boost)), drop = FALSE]

  out_str <- "["
  for (i in 1:nrow(ram_df_boost)) {
    out_str <- paste0(out_str, "[", paste0(ram_df_boost[i,], collapse = ", "), "]", if (i == nrow(ram_df_boost)) "]" else ",", "\n")
  }
  cat(out_str)
}


#' Get sparse matrix pieces
#'
#' @param csc_mat a sparse matrix in CSC format
#'
#' @return a list containing x, i, p, Dim(1), and Dim(2) (in that order)
get_sparse_matrix_pieces <- function(csc_mat) {
  list(csc_mat@x, csc_mat@i, csc_mat@p, csc_mat@Dim[1], csc_mat@Dim[2])
}


#' Get the name of a gRNA dataset
#'
#' @param dataset_name the name of a dataset in paper/experiment/modality format
#' @param gRNA_modality the name of the gRNA modality, one of "assignment" or "expression"
#'
#' @return the name of the gRNA dataset in paper/experiment/modality format
#' @export
get_gRNA_dataset_name <- function(dataset_name, gRNA_modality) {
  paste0(sub('/[^/]*$', '', dataset_name), "/grna_", gRNA_modality)
}


#' Get undercover groups
#'
#' Produce a list of undercover groups given a list of NTCs.
#'
#' @param ntc_names a vector of NTC names
#' @param group_size size of the undercover groups
#' @param partition_count number of undercover groups to generate
#'
#' @return a vector of comma-separated undercover groups
#' @export
get_undercover_groups <- function(ntc_names, group_size, partition_count) {
  set.seed(4)
  n_ntcs <- length(ntc_names)
  total_possible_paritions <- choose(n_ntcs, group_size)
  if (partition_count > total_possible_paritions) stop("partition_count exceeds the total number of possible partitions.")
  my_undercover_groups <- character()
  ntc_names_copy <- ntc_names
  repeat {
    while (length(ntc_names_copy) >= group_size) {
      curr_grp <- sample(x = ntc_names_copy,
                         size = group_size,
                         replace = FALSE)
      curr_grp_string <- curr_grp |> sort() |> paste0(collapse = ",")
      if (!(curr_grp_string %in% my_undercover_groups)) {
        my_undercover_groups <- c(my_undercover_groups, curr_grp_string)
      }
      ntc_names_copy <- ntc_names_copy[!(ntc_names_copy %in% curr_grp)]
    }
    if (length (my_undercover_groups) >= partition_count) break()
    ntc_names_copy <- ntc_names
  }
  return(my_undercover_groups[seq(1, partition_count)])
}


#' Get gRNA assignments via max operation
#'
#' @param gRNA_odm a gRNA ODM;
#'
#' @return a vector of gRNA IDs obtained by applying a column-wise max operation to the gRNA ODM.
get_gRNA_assignments_via_max_op <- function(gRNA_odm) {
  ret <- gRNA_odm |>
    load_whole_odm() |>
    apply(MARGIN = 2, FUN = function(col) names(which.max(col))) |>
    unname()
  return(ret)
}


#' Get target assignments via max operation
#'
#' @param gRNA_odm a gRNA ODM
#'
#' @return a vector of "targets" obtained by applying a column-wise max operation to the gRNA ODM;
#' the target is present as a column of the feature covariate matrix of the gRNA ODM.
#' @export
#'
#' @examples
#' \dontrun{
#' gRNA_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' get_target_assignments_via_max_op(gRNA_odm)
#' }
get_target_assignments_via_max_op <- function(gRNA_odm) {
  gRNA_feature_covariates <- gRNA_odm |> ondisc::get_feature_covariates()
  gRNA_to_target_map <- stats::setNames(row.names(gRNA_feature_covariates), gRNA_feature_covariates$target)
  gRNA_assignments <- get_gRNA_assignments_via_max_op(gRNA_odm)
  gRNA_targets <- names(gRNA_to_target_map)[match(x = gRNA_assignments, table = gRNA_to_target_map)]
  return(gRNA_targets)
}


#' Get library sizes
#'
#' Gets the cell-wise library sizes from a response ODM (either "n_umis" or "n_fragments", as appropriate)
#'
#' @param response_odm a response ODM
#'
#' @return a vector of cell library sizes
get_library_sizes <- function(response_odm) {
  cell_covariates <- response_odm |> ondisc::get_cell_covariates()
  lib_size_cov <- if ("n_umis" %in% colnames(cell_covariates)) {
    "n_umis"
  } else if ("n_fragments" %in% colnames(cell_covariates)) {
    "n_fragments"
  } else {
    stop("Neither n_umis nor n_fragments are columns of the cell covariates matrix of response_odm.")
  }
  cell_covariates[[lib_size_cov]]
}


#' Read all modalities
#'
#' @param paper name of the paper to load
#' @param dataset name of the dataset within the paper to load
#' @param sceptre2_data_dir location of the sceptre2 data directory
#'
#' @return a multimodal odm for the specified paper/dataset
#' @export
read_all_modalities <- function(paper, dataset, sceptre2_data_dir = paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_vect <- list.files(dataset_dir)
  odm_list <- list()
  for (modality in modality_vect) {
    odm_dir <- paste0(dataset_dir, "/", modality, "/")
    curr_odm <- read_odm(odm_fp = paste0(odm_dir, "matrix.odm"),
                         metadata_fp = paste0(odm_dir, "metadata_orig.rds"))
    odm_list <- c(odm_list, curr_odm)
  }
  names(odm_list) <- modality_vect
  ret <- multimodal_ondisc_matrix(odm_list)
  return(ret)
}
