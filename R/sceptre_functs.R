#' sceptre (v2)
#'
#' Implements sceptre2.
#'
#' @inherit abstract_interface
#' @export
sceptre <- function(response_odm, grna_odm, response_grna_group_pairs) {
  # construct the multimodal ODM
  mm_odm <- ondisc::multimodal_ondisc_matrix(list(response = response_odm,
                                                  grna = grna_odm)) |>
    lowmoi::process_multimodal_odm()

  # get the formula
  form <- mm_odm@modalities$response@misc$sceptre_formula

  # call function
  res <- sceptre2::run_sceptre_low_moi(mm_odm = mm_odm,
                                       response_grna_group_pairs = response_grna_group_pairs,
                                       form = form,
                                       response_modality_name = "response",
                                       grna_modality_name = "grna",
                                       grna_group_column_name = "target",
                                       B = 2500,
                                       side = "both",
                                       full_output = FALSE)
}
