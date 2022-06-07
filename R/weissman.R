#' Weissman method
#' Implements the Weissman method.
#'
#' @inherit abstract_interface
#' @export
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' gRNA_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna")
#' response_gRNA_group_pairs <- expand.grid(response_id = (response_odm |> ondisc::get_feature_ids()), gRNA_group = c("GATA1-C", "GATA1-D"))
#' result <- weissman(response_odm, gRNA_odm, response_gRNA_group_pairs)
#' }
weissman <- function(response_odm, response_gRNA_group_pairs, gRNA_assignment_method = "original", gRNA_odm = NULL) {
  # get perturbation assignments
  pert_assignments <- switch(gRNA_assignment_method,
    original = {
      response_odm |>
        ondisc::get_cell_covariates() |>
        pull(perturbation)
    },
    {
      stop("Invalid specification of gRNA_assignment_method.")
    }
  )

  # convert to CellPopulation format
  cell_pop <- odm_to_cell_pop(response_odm, pert_assignments)

  # normalize data
  cell_pop$normalized_matrix <- normalize_to_gemgroup_control(
    pop = cell_pop,
    control_cells = 'guide_target == "non-targeting"'
  )

  # apply Weissman method
  unique_gRNAs <- as.character(unique(response_gRNA_group_pairs$gRNA_group))
  res_list <- lapply(unique_gRNAs, function(curr_gRNA) {
    # find responses to test this gRNA against
    response_vars <- response_gRNA_group_pairs |>
      dplyr::filter(gRNA_group == curr_gRNA) |>
      dplyr::pull(response_id) |>
      as.character()

    # get normalized matrix for cells with this gRNA and genes to test against
    cell_pop_targeting <- cell_pop$where(
      cells = sprintf('guide_identity == \"%s\"', curr_gRNA),
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
    stopifnot(all(names(p_vals) == response_vars |> add_ENSG()))
    out_df <- data.frame(
      response_id = response_vars,
      gRNA_group = curr_gRNA,
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
#' @param pert_assignments Character vector of perturbation assignments, one per cell
#'
#' @return An object of type `CellPopulation` that the `perturbseq` library expects
odm_to_cell_pop <- function(response_odm, pert_assignments) {
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
    dplyr::mutate(guide_identity = pert_assignments) # add gRNA assignments
  if (!("batch" %in% names(cells_df))) {
    cells_df$batch <- 1 # add batch information if it is not present
  }
  cells_df <- cells_df |>
    dplyr::left_join(gRNA_odm |> # join with gRNA metadata to get guide targets
                       ondisc::get_feature_covariates() |>
                       tibble::rownames_to_column(var = "guide_identity") |>
                       dplyr::select(guide_identity, target),
                     by = "guide_identity"
    ) |>
    dplyr::rename(UMI_count = n_umis, gem_group = batch, guide_target = target) |>
    dplyr::mutate(gem_group = as.integer(as.factor(gem_group))) |>
    dplyr::select(cell_barcode, UMI_count, gem_group, guide_identity, guide_target) |>
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
    if (substr(feature_id, start = 0, stop = 4) != "ENSG") {
      paste0("ENSG_", feature_id)
    } else {
      feature_id
    }
  }) |> unname()
}
