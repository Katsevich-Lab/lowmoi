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
  modalities <- modalities[modalities %in% c("gene", "protein", "grna_assignment", "grna_expression", "chromatin")]
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


#' Get the name of a grna dataset
#'
#' @param dataset_name the name of a dataset in paper/experiment/modality format
#' @param grna_modality the name of the grna modality, one of "assignment" or "expression"
#'
#' @return the name of the grna dataset in paper/experiment/modality format
#' @export
get_grna_dataset_name <- function(dataset_name, grna_modality) {
  paste0(sub('/[^/]*$', '', dataset_name), "/grna_", grna_modality)
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


#' Get grna assignments via max operation
#'
#' @param grna_odm a grna ODM;
#'
#' @return a vector of grna IDs obtained by applying a column-wise max operation to the grna ODM.
get_grna_assignments_via_max_op <- function(grna_odm) {
  ret <- grna_odm |>
    load_whole_odm() |>
    apply(MARGIN = 2, FUN = function(col) names(which.max(col))) |>
    unname()
  return(ret)
}


#' Get target assignments via max operation
#'
#' @param grna_odm a grna ODM
#'
#' @return a vector of "targets" obtained by applying a column-wise max operation to the grna ODM;
#' the target is present as a column of the feature covariate matrix of the grna ODM.
#' @export
#'
#' @examples
#' \dontrun{
#' grna_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' get_target_assignments_via_max_op(grna_odm)
#' }
get_target_assignments_via_max_op <- function(grna_odm) {
  grna_feature_covariates <- grna_odm |> ondisc::get_feature_covariates()
  grna_feature_covariates$target[is.na(grna_feature_covariates$target)] <- "candidate"
  grna_to_target_map <- stats::setNames(row.names(grna_feature_covariates), grna_feature_covariates$target)
  grna_assignments <- get_grna_assignments_via_max_op(grna_odm)
  grna_targets <- names(grna_to_target_map)[match(x = grna_assignments, table = grna_to_target_map)]
  return(grna_targets)
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
  modality_vect <- modality_vect[modality_vect %in% c("gene", "grna_assignment", "grna_expression", "chromatin", "protein")]
  odm_list <- list()
  for (modality in modality_vect) {
    odm_dir <- paste0(dataset_dir, "/", modality, "/")
    curr_odm <- ondisc::read_odm(odm_fp = paste0(odm_dir, "matrix.odm"),
                         metadata_fp = paste0(odm_dir, "metadata_orig.rds"))
    odm_list <- c(odm_list, curr_odm)
  }
  names(odm_list) <- modality_vect
  ret <- ondisc::multimodal_ondisc_matrix(odm_list)
  return(ret)
}


#' Save all modalities
#'
#' @param multimodal_odm a multimodal ODM object
#' @param paper name of the paper
#' @param dataset name of the dataset within the paper
#' @param metadata_file_name name of the metadata file to write (e.g., "metadata_qc.rds")
#' @param sceptre2_data_dir location of the sceptre2 data directory
#'
#' @return NULL
#' @export
save_all_modalities <- function(multimodal_odm, paper, dataset, metadata_file_name, sceptre2_data_dir = paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/")) {
  dataset_dir <- paste0(sceptre2_data_dir, paper, "/", dataset)
  modality_list <- multimodal_odm@modalities |> names()
  for (modality in modality_list) {
    ondisc::save_odm(odm = ondisc::get_modality(multimodal_odm, modality),
             metadata_fp = paste0(dataset_dir, "/", modality, "/", metadata_file_name))
  }
}

#' Get dataset for perturbation propensity analysis
#'
#' Get dataset for perturbation propensity analysis.
#'
#' @param paper "frangieh", "papalexi", etc.
#' @param dataset "ifn_gamma", "eccite_screen", etc.
#'
#' @return Data frame that contains the grna assignments and technical covariates.
#' @export
get_data_for_pert_prop <- function(paper, dataset) {
  # get grna ODM
  grna_fp <- paste0(paper, "/", dataset, "/grna_assignment")
  grna_odm <- lowmoi::load_dataset_modality(data_fp = grna_fp)

  # get response ODM
  if (paper == "liscovitch") {
    response_fp <- paste0(paper, "/", dataset, "/chromatin")
  } else {
    response_fp <- paste0(paper, "/", dataset, "/gene")
  }
  response_odm <- lowmoi::load_dataset_modality(data_fp = response_fp)

  # get list of NTC grnas
  ntcs <- grna_odm |>
    ondisc::get_feature_covariates() |>
    dplyr::filter(target == "non-targeting") |>
    rownames() |>
    unique()

  # get the grna assigned to each cell
  assigned_grnas_df <- grna_odm |>
    ondisc::get_cell_covariates() |>
    dplyr::select(assigned_grna)

  # get a logical vector of whether a cell was assigned an NTC
  ntc_cells <- assigned_grnas_df |>
    dplyr::mutate(ntc = assigned_grna %in% ntcs) |>
    dplyr::pull(ntc)

  # restrict assigned grnas to cells with an NTC
  assigned_grnas <- assigned_grnas_df[ntc_cells, , drop = F] |> dplyr::pull(assigned_grna)

  # extract cell covariates from gene ODM, and restrict attention to cells with NTC
  cell_covariates <- response_odm[, ntc_cells] |> ondisc::get_cell_covariates()

  # append assigned_grnas column to cell_covariates
  df <- cell_covariates |>
    dplyr::mutate(assigned_grna = assigned_grnas)

  # return
  df
}


#' Process multimodal odm
#'
#' @param mm_odm a multimodal odm
#' @param remove_grna_assignment_n_nonzero remove the column `grna_assignment_n_nonzero` (as it coincides with `grna_expression_n_nonzero`)?
#'
#' @return a streamlined, simpler multimodal odm
#' @export
process_multimodal_odm <- function(mm_odm, remove_grna_assignment_n_nonzero = TRUE) {
  cell_covariate_m <- mm_odm |> ondisc::get_cell_covariates()
  if (remove_grna_assignment_n_nonzero) {
    cell_covariate_m <- cell_covariate_m |> dplyr::mutate(grna_assignment_n_nonzero = NULL,
                                                          grna_assignment_n_umis = NULL)
  } else {
    cell_covariate_m <- cell_covariate_m |> dplyr::mutate(grna_assignment_n_umis = NULL)
  }

  cell_covariate_colnames <- colnames(cell_covariate_m)
  shared_covariates <- c("batch", "p_mito", "phase", "bio_rep", "lane")
  for (shared_covariate in shared_covariates) {
    match_col_names <- grep(pattern = shared_covariate, x = cell_covariate_colnames, value = TRUE)
    if (length(match_col_names) >= 1) {
      cell_covariate_m[[shared_covariate]] <- cell_covariate_m[[match_col_names[1]]]
      for (match_col_name in match_col_names) cell_covariate_m[[match_col_name]] <- NULL
    }
  }
  mm_odm@global_cell_covariates <- cell_covariate_m
  return(mm_odm)
}


#' Get sceptre function args for pair
#'
#' For a given response ID and undercover gRNA, obtains the arguments passed (by do.call) to sceptre in the undercover pipeline.
#'
#' @param response_id ID of the response
#' @param undercover_grna ID of the gRNA to go undercover
#' @param dataset_name name of the dataset
#' @param output_amount sceptre output amount
#' @param B number of resamples B
#'
#' @return a list of (swapped) arguments to pass (via do.call) to sceptre;
#' @export
get_sceptre_function_args_for_pair <- function(response_id, undercover_grna, dataset_name, output_amount, B) {
  undercover_ntc_name_in <- undercover_grna
  response_odm <- lowmoi::load_dataset_modality(dataset_name)
  grna_dataset_name <- lowmoi::get_grna_dataset_name(dataset_name, "assignment")
  grna_odm <- lowmoi::load_dataset_modality(grna_dataset_name)

  undercover_ntc_name <- strsplit(x = undercover_ntc_name_in, split = ",", fixed = TRUE) |> unlist()
  grna_feature_covariates <- grna_odm |> ondisc::get_feature_covariates()
  grna_feature_covariates[undercover_ntc_name, "target"] <- "undercover"
  if (!("non-targeting" %in% grna_feature_covariates$target)) {
    stop("After performing label swap, `non-targeting` is no longer string in the `target` column.")
  }
  grna_odm_swapped <- grna_odm |> ondisc::mutate_feature_covariates(target = grna_feature_covariates$target)
  response_grna_group_pairs <- data.frame(response_id = response_id,
                                          grna_group = "undercover")
  response_odm <- response_odm; grna_odm <- grna_odm_swapped; response_grna_group_pairs <- response_grna_group_pairs
  mm_odm <- lowmoi::process_multimodal_odm(ondisc::multimodal_ondisc_matrix(list(response = response_odm, grna = grna_odm)))
  form <- mm_odm@modalities$response@misc$sceptre_formula

  ret <- list(mm_odm = mm_odm,
              response_grna_group_pairs = response_grna_group_pairs,
              form = form,
              response_modality_name = "response",
              grna_modality_name = "grna",
              grna_group_column_name = "target",
              B = B,
              side = "both",
              output_amount = output_amount)
  return(ret)
}



#' Generate all pairs
#'
#' Generates the entire set of response-gRNA group pairs from an input response ODM and gRNA ODM.
#'
#' @param response_odm a response ODM
#' @param grna_odm a gRNA ODM
#'
#' @return a data frame with columns "response_id" and "grna_group" containing all pairs
#' @export
generate_all_pairs <- function(response_odm, grna_odm) {
  response_ids <- ondisc::get_feature_ids(response_odm)
  grna_group_data_frame <- grna_odm |> get_feature_covariates()
  grna_group_data_frame <- data.frame(grna_id = rownames(grna_group_data_frame),
                                      grna_group = grna_group_data_frame$target)
  grna_groups <-  grna_group_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |>
    dplyr::pull(grna_group) |> unique() |> factor()
  expand.grid(response_id = response_ids, grna_group = grna_groups)
}
