#' Process undercover results
#'
#' Processes undercover results. This function (i) replashes slashes with underscores in the dataset names, (ii) combines the Schraivogel enhancer screens, and (iii) appends n treatment (i.e., number of nonzero cells containing the undercover gRNA), n control (i.e., number of nonzero cells containing an NT gRNA, excluding the undercover gRNA), and effective_sample_size (i.e., total number of nonzero cells containing an NT gRNA)
#'
#' @param undercover_res an output of the undercover pipeline
#' @param sample_size_df the sample sizes data frame
#'
#' @return a proccessed `undercover_res`
#' @export
process_undercover_result <- function(undercover_res, sample_size_df) {
  # basic processing on the undercover results
  x <- undercover_res |>
    replace_slash_w_underscore() |>
    combine_schraivogel_enhancer_screens() |>
    update_dataset_and_method_names()

  # wrange the sample size df
  sample_size_df_nt <- sample_size_df |>
    dplyr::filter(target == "non-targeting") |>
    dplyr::mutate(dataset = dataset_concat,
                  dataset_concat = NULL, paper = NULL, modality = NULL) |>
    dplyr::rename(grna_group = target, response_id = feature_id) |>
    replace_slash_w_underscore() |>
    combine_schraivogel_enhancer_screens()
  ntc_effective_samp_size <- sample_size_df_nt |>
    dplyr::group_by(response_id, dataset) |>
    dplyr::summarize(effective_samp_size = sum(n_nonzero_cells))
  sample_size_df_nt <- dplyr::left_join(x = sample_size_df_nt,
                                        y = ntc_effective_samp_size,
                                        by = c("response_id", "dataset")) |>
    dplyr::mutate(n_nonzero_treatment = n_nonzero_cells,
                  n_nonzero_control = effective_samp_size - n_nonzero_cells) |>
    dplyr::rename("undercover_grna" = "grna_id") |>
    dplyr::select(response_id, undercover_grna, dataset,
                  n_nonzero_treatment, n_nonzero_control, effective_samp_size)

  out <- dplyr::left_join(x = x,
                          y = sample_size_df_nt,
                          by = c("response_id", "undercover_grna", "dataset"))
  return(out)
}

replace_slash_w_underscore <- function(undercover_res) {
  undercover_res |> dplyr::mutate(dataset = gsub(pattern = "/",
                                                 replacement = "_",
                                                 fixed = TRUE,
                                                 x = dataset))
}

combine_schraivogel_enhancer_screens <- function(undercover_res) {
  undercover_res |> dplyr::mutate(dataset = dataset |> forcats::fct_recode(schraivogel_enhancer_screen = "schraivogel_enhancer_screen_chr11_gene",
                                                                           schraivogel_enhancer_screen = "schraivogel_enhancer_screen_chr8_gene"))
}

update_dataset_and_method_names <- function(undercover_res) {
  undercover_res |>
    dplyr::mutate(dataset_rename = stringr::str_to_title(gsub(pattern = "_", replacement = " ", x = dataset)),
                  Method = stringr::str_to_title(gsub(pattern = "_", replacement = " ", x = method)))
}
