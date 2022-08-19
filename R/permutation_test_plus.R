#' Permutation test (plus)
#'
#' Carries out a more sophisticated permutation test by (i) regressing the expressions onto the technical factors, (ii) using a (score-based) z-score test statistic.
#'
#' @inherit abstract_interface
#' @param n_rep number of permutation replicates
#' @param return_permuted_test_stats return the permuted test statistics?
#' @param progress print progress messages?
#' @export
permutation_test_plus <- function(response_odm, grna_odm, response_grna_group_pairs, n_rep = 1000, return_permuted_test_stats = TRUE, progress = TRUE) {
  # get the cell covariates; define the model formula
  cell_covariate_df <- response_odm |> ondisc::get_cell_covariates()
  my_formula_str <- response_odm@misc$nb_regression_formula
  my_formula <- stats::as.formula(paste0("expression ", my_formula_str))

  # define the two sample test
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    df <- create_design_matrix(target_cells, control_cells, target_cell_indices, control_cell_indices, cell_covariate_df)
    # first, use aux function to estimate size
    gene_precomp_size <- max(estimate_size(df, my_formula), 0.01)
    # fit GLM of expressions on technical factors; extract the fitted values
    fit_nb <- VGAM::vglm(formula = my_formula, family = VGAM::negbinomial.size(gene_precomp_size), data = df)
    exp_gene_offsets <- as.numeric(fit_nb@fitted.values)

    # obtain the ground truth z-score
    expressions <- df$expression
    grna_indicators <- df$pert_indicator
    z_star <- sceptre:::compute_nb_test_stat_fast_score(y = expressions[grna_indicators == 1],
                                                        exp_o = exp_gene_offsets[grna_indicators == 1],
                                                        gene_precomp_size = gene_precomp_size)

    # compute permuted test statistics
    z_null <- replicate(n = n_rep, expr = {
      grna_indicators <- sample(df$pert_indicator)
      sceptre:::compute_nb_test_stat_fast_score(y = expressions[grna_indicators == 1],
                                                exp_o = exp_gene_offsets[grna_indicators == 1],
                                                gene_precomp_size = gene_precomp_size)
    })

    out <- run_skew_t_prepare_output(z_null, z_star, return_permuted_test_stats)
    return(out)
  }

  # run the permutation test
  cbind_res <- if (return_permuted_test_stats) TRUE else FALSE
  res <- abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress, cbind_res)
  return(res)
}


run_skew_t_prepare_output <- function(z_null, z_star, return_permuted_test_stats) {
  # compute the p-value via call to fit skew t (in sceptre)
  p_val_fit <- sceptre:::fit_skew_t(z_null, z_star, "both")
  if (return_permuted_test_stats) {
    resamples_df <- as.data.frame(matrix(z_null, nrow = 1, ncol = length(z_null)))
    colnames(resamples_df) <- paste0("z_null_", seq(1, length(z_null)))
    out <- cbind(data.frame(p_value = p_val_fit$out_p, z_value = z_star,
                            skew_t_fit_success = p_val_fit$skew_t_fit_success,
                            xi = p_val_fit$skew_t_mle[["xi"]], omega = p_val_fit$skew_t_mle[["omega"]],
                            alpha = p_val_fit$skew_t_mle[["alpha"]], nu = p_val_fit$skew_t_mle[["nu"]]), resamples_df)
  } else {
    out <- p_val_fit$out_p
  }
  return(out)
}
