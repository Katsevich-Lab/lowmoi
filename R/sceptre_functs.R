#' sceptre (v2)
#'
#' Implements sceptre2.
#'
#' @inherit abstract_interface
#' @param B number of resamples to draw
#' @param output_amount amount of output to return from function (ranging from 1-3)
#' @export
sceptre <- function(response_odm, grna_odm, response_grna_group_pairs, B = 2500, output_amount = 1, sn_approx = TRUE) {
  if (!is.numeric(B)) B <- as.integer(B)
  if (!is.numeric(output_amount)) output_amount <- as.integer(output_amount)
  if (!is.logical(sn_approx)) sn_approx <- as.logical(sn_approx)

  # construct the multimodal ODM
  mm_odm <- ondisc::multimodal_ondisc_matrix(list(response = response_odm, grna = grna_odm)) |>
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
                                       B = B,
                                       side = "both",
                                       output_amount = output_amount,
                                       sn_approx = sn_approx) |> as.data.frame()

  # select p_val, grna_grop, response_id
  if (output_amount == 1) {
    res <- res |> dplyr::select(response_id, grna_group, p_value)
  }

  return(res)
}
