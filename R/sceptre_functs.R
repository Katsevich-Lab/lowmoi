#' Implements sceptre (v3)
#'
#' @inherit abstract_interface
#'
#' @export
sceptre <- function(response_odm, grna_odm, response_grna_group_pairs, with_covariates = TRUE, regression_method = "poisson_glm", discovery_test_stat = "approximate", print_progress = TRUE) {
  # load the gene and grna matrix into memory
  response_matrix <- load_whole_odm(response_odm)
  grna_matrix <- load_whole_odm(grna_odm)
  grna_feature_df <- grna_odm |> ondisc::get_feature_covariates()
  grna_group_data_frame <- data.frame(grna_id = rownames(grna_feature_df),
                                      grna_group = grna_feature_df$target)
  # construct the multimodal odm
  mm_odm <- ondisc::multimodal_ondisc_matrix(list(response = response_odm, grna = grna_odm)) |>
    lowmoi::process_multimodal_odm()

  covariate_data_frame <- mm_odm |> ondisc::get_cell_covariates()
  formula_object <- response_odm@misc$sceptre_formula

  # get the formula
  if (with_covariates) {
    form <- mm_odm@modalities$response@misc$sceptre_formula
  } else {
    form <- stats::formula(~log(response_n_umis))
  }

  # run sceptre
  res <- sceptre::run_sceptre_lowmoi(response_matrix = response_matrix,
                              grna_matrix = grna_matrix,
                              covariate_data_frame = covariate_data_frame,
                              grna_group_data_frame = grna_group_data_frame,
                              formula_object = formula_object,
                              response_grna_group_pairs = response_grna_group_pairs,
                              calibration_check = FALSE,
                              print_progress = print_progress,
                              regression_method = regression_method,
                              discovery_test_stat = discovery_test_stat) |>
    dplyr::select(response_id, grna_group, p_value)
  return(res)
}


sceptre_no_covariates <- function(response_odm, grna_odm, response_grna_group_pairs) {
  sceptre(response_odm, grna_odm, response_grna_group_pairs)
}


#' @export
#' @inherit abstract_interface
sceptre_approximate_poisson <- function(response_odm, grna_odm, response_grna_group_pairs) {
  sceptre(response_odm = response_odm,
          grna_odm = grna_odm,
          response_grna_group_pairs = response_grna_group_pairs,
          with_covariates = TRUE,
          regression_method = "poisson_glm",
          discovery_test_stat = "approximate",
          print_progress = TRUE)
}


#' @export
#' @inherit abstract_interface
sceptre_exact_poisson <- function(response_odm, grna_odm, response_grna_group_pairs) {
  sceptre(response_odm = response_odm,
          grna_odm = grna_odm,
          response_grna_group_pairs = response_grna_group_pairs,
          with_covariates = TRUE,
          regression_method = "poisson_glm",
          discovery_test_stat = "exact",
          print_progress = TRUE)
}

#' @export
#' @inherit abstract_interface
sceptre_approximate_nb <- function(response_odm, grna_odm, response_grna_group_pairs) {
  sceptre(response_odm = response_odm,
          grna_odm = grna_odm,
          response_grna_group_pairs = response_grna_group_pairs,
          with_covariates = TRUE,
          regression_method = "nb_glm",
          discovery_test_stat = "approximate",
          print_progress = TRUE)
}

#' @export
#' @inherit abstract_interface
sceptre_exact_nb <- function(response_odm, grna_odm, response_grna_group_pairs) {
  sceptre(response_odm = response_odm,
          grna_odm = grna_odm,
          response_grna_group_pairs = response_grna_group_pairs,
          with_covariates = TRUE,
          regression_method = "nb_glm",
          discovery_test_stat = "exact",
          print_progress = TRUE)
}
