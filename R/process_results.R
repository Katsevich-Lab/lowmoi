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
  # wrange the sample size df
  sample_size_df_nt <- sample_size_df |>
    dplyr::filter(grna_group == "non-targeting")
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
  out <- dplyr::left_join(x = undercover_res,
                          y = sample_size_df_nt,
                          by = c("response_id", "undercover_grna", "dataset")) |>
    replace_slash_w_underscore()

  if ("schraivogel_enhancer_screen_chr11_gene" %in% out$dataset &&
      "schraivogel_enhancer_screen_chr8_gene" %in% out$dataset) {
    out <- out |> combine_schraivogel_enhancer_screens()
  }

  out <- out |>
    update_dataset_and_method_names() |>
    dplyr::mutate(p_value = ifelse(p_value <= 1e-250, 1e-250, p_value))
  return(out)
}


#' Process positive control result
#'
#' @param pc_res positive control result data frame
#' @param sample_size_df sample size data frame
#'
#' @return a processed `pc_res`
#' @export
process_pc_result <- function(pc_res, sample_size_df) {
  sample_size_df_pc <- sample_size_df |>
    filter(response_id %in% unique(pc_res$response_id),
           dataset %in% unique(pc_res$dataset))
  control_sample_size_df <- sample_size_df |>
    filter(grna_group == "non-targeting") |>
    group_by(response_id, dataset) |>
    summarize(n_control = sum(n_nonzero_cells))
  to_join <- sample_size_df_pc |>
    group_by(grna_group, response_id, dataset) |>
    summarize(n_treatment = sum(n_nonzero_cells)) |>
    select(response_id, grna_group, n_treatment, dataset)
  pc_res_w_ss <- left_join(x = pc_res,
                           y = to_join,
                           by = c("grna_group", "response_id", "dataset")) |>
    left_join(y = control_sample_size_df, by = c("response_id", "dataset")) |>
    replace_slash_w_underscore()
  if ("schraivogel_enhancer_screen_chr8_gene" %in% pc_res_w_ss$dataset &&
      "schraivogel_enhancer_screen_chr11_gene" %in% pc_res_w_ss$dataset) {
    pc_res_w_ss <- pc_res_w_ss |> combine_schraivogel_enhancer_screens()
  }
  pc_res_w_ss <- pc_res_w_ss |>
    update_dataset_and_method_names() |>
    dplyr::mutate(p_value = ifelse(p_value <= 1e-250, 1e-250, p_value))
  return(pc_res_w_ss)
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
    mutate(Method = forcats::fct_recode(method,
                                        "SCEPTRE" = "sceptre",
                                        "SCEPTRE (no covariates)" = "sceptre_no_covariates",
                                        "NB regression (w/ covariates)" = "nb_regression_w_covariates",
                                        "NB regression (no covariates)" = "nb_regression_no_covariates",
                                        "Seurat-Wilcox" = "seurat_de",
                                        "Seurat-NB" = "seurat_de_nb",
                                        "t-test" = "liscovitch_method",
                                        "KS test" = "weissman_method",
                                        "MAST" = "schraivogel_method",
                                        "MIMOSCA" = "mimosca") |>
             forcats::fct_relevel(c("SCEPTRE",
                                    "SCEPTRE (no covariates)",
                                    "NB regression (no covariates)",
                                    "NB regression (w/ covariates)",
                                    "Seurat-Wilcox",
                                    "Seurat-NB",
                                    "t-test",
                                    "KS test",
                                    "MAST",
                                    "MIMOSCA"))) |>
    mutate(dataset_rename = forcats::fct_recode(dataset,
                                                "Frangieh (Co Culture)" = "frangieh_co_culture_gene",
                                                "Frangieh (Control)" = "frangieh_control_gene",
                                                "Frangieh (IFN-\u03B3)" = "frangieh_ifn_gamma_gene",
                                                "Papalexi (Gene)" = "papalexi_eccite_screen_gene",
                                                "Papalexi (Protein)" = "papalexi_eccite_screen_protein",
                                                "Schraivogel" = "schraivogel_enhancer_screen",
                                                "Simulated" = "simulated_experiment_1_gene"
    ))
  }
