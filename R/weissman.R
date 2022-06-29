#' Weissman method
#' Implements the Weissman method.
#'
#' @inherit abstract_interface
#' @export
weissman_method <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  # obtain the cell-to-target assignments
  if (gRNA_odm@ondisc_matrix@logical_mat) {
    guide_targets <- get_target_assignments_via_max_op(gRNA_odm)
  } else{
    stop("The Weissman gRNA assignment method is not currently implemented.")
  }

  # convert to CellPopulation format
  cell_pop <- odm_to_cell_pop(response_odm, guide_targets)

  # normalize data
  cell_pop$normalized_matrix <- normalize_to_gemgroup_control(
    pop = cell_pop,
    control_cells = 'guide_target == "non-targeting"'
  )

  # apply Weissman method
  unique_gRNA_groups <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(unique_gRNA_groups, function(curr_gRNA_group) {
    # find responses to test this gRNA against
    response_vars <- response_gRNA_group_pairs |>
      dplyr::filter(gRNA_group == curr_gRNA_group) |>
      dplyr::pull(response_id) |>
      as.character()

    # get normalized matrix for cells with this gRNA and genes to test against
    cell_pop_targeting <- cell_pop$where(
      cells = sprintf('guide_target == \"%s\"', curr_gRNA_group),
      genes = response_vars |> add_ENSG(),
      normalized = TRUE
    )

    # get normalized matrix for cells with control gRNAs and genes to test against
    cell_pop_control <- cell_pop$where(
      cells = 'guide_target == "non-targeting"',
      genes = response_vars |> add_ENSG(),
      normalized = TRUE
    )

    # apply the python function ks_compare_pops
    p_vals <- ks_compare_pops(cell_pop_targeting, cell_pop_control)[[2]]
    out_df <- data.frame(
      response_id = names(p_vals) |> sapply(function(str)(substr(str, 6, nchar(str)))) |> unname(),
      gRNA_group = curr_gRNA_group,
      p_value = unname(p_vals)
    )
    out_df
  })
  # concatenate and return results
  dplyr::left_join(
    response_gRNA_group_pairs,
    do.call(what = "rbind", args = res_list)
  )
}

#' Convert from ODM to CellPopulation format
#'
#' @param response_odm ODM for response variable
#' @param guide_targets Character vector of perturbation assignments, one per cell
#'
#' @return An object of type `CellPopulation` that the `perturbseq` library expects
odm_to_cell_pop <- function(response_odm, guide_targets) {
  # load the data, transposing the matrices
  response_mat_t <- response_odm |>
    load_whole_odm() |>
    Matrix::t() |>
    as.matrix() |>
    as.data.frame()

  # add ENSG to column names for compatibility with perturbseq software
  colnames(response_mat_t) <- colnames(response_mat_t) |> add_ENSG()
  feature_ids <- ondisc::get_feature_ids(response_odm)
  feature_names <- ondisc::get_feature_names(response_odm)
  # if feature names are not present, let feature ids be feature names
  if (all(is.na(feature_names))) {
    feature_names <- feature_ids
  }
  feature_ids <- feature_ids |> add_ENSG()

  # create genes_df for input to CellPopulation constructor
  genes_df <- data.frame(gene_name = feature_names, row.names = feature_ids)

  # create cells_df for input to CellPopulation constructor
  cells_df <- ondisc::get_cell_covariates(response_odm) |>
    tibble::rownames_to_column(var = "cell_barcode") |>
    dplyr::mutate(guide_target = guide_targets) # add gRNA assignments
  if (!("batch" %in% names(cells_df))) {
    cells_df$batch <- 1 # add batch information if it is not present
  }
  # # join with guide_to_target_map
  cells_df <- cells_df |>
    dplyr::rename(UMI_count = n_umis, gem_group = batch) |>
    dplyr::mutate(gem_group = as.integer(as.factor(gem_group))) |>
    dplyr::select(cell_barcode, UMI_count, gem_group, guide_target) |>
    tibble::column_to_rownames(var = "cell_barcode")

  # construct and return the CellPopulation object expected by perturbseq software
  CellPopulation(
    matrix = response_mat_t,
    cell_list = cells_df,
    gene_list = genes_df
  )
}

#' Add ENSG prefix to a list of feature ids if necessary
#'
#' This is somewhat of a hack needed for compatibility with perturbseq code.
#'
#' @param feature_id Vector of feature ids
#'
#' @return The same set of feature ids with ENSG prefixes added if necessary
add_ENSG <- function(feature_ids) {
  sapply(feature_ids, function(feature_id) {
    paste0("ENSG_", feature_id)
  }) |> unname()
}
