#' NB regression
#'
#' Implements a negative binomial regression.
#'
#' @inherit abstract_interface
#' @export
nb_regression <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  # obtain cell covariate data frame
  cell_covariate_df <- response_odm |> ondisc::get_cell_covariates()
  my_formula_str <- response_odm@misc$nb_regression_formula
  my_formula <- stats::as.formula(paste0("expression ", my_formula_str, " + pert_indicator"))

  # define the NB regression test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    # construct the data matrix to pass to GLM
    df_left <- data.frame(expression = c(target_cells, control_cells),
                          pert_indicator = c(rep(1, length(target_cells)),
                                             rep(0, length(control_cells))))
    df_right <- rbind(cell_covariate_df[target_cell_indices,], cell_covariate_df[control_cell_indices,])
    df <- cbind(df_left, df_right)

    # first, use aux function to estimate size
    est_size <- max(estimate_size(df, my_formula), 0.01)

    # fit GLM
    fit_nb <- VGAM::vglm(formula = my_formula, family = VGAM::negbinomial.size(est_size), data = df)

    # extract p-value
    s <- VGAM::summary(fit_nb)
    p_val <- s@coef3["pert_indicator", "Pr(>|z|)"]
    return(p_val)
  }

  # run the NB regression on all the data
  res <- abstract_two_sample_test(response_odm, gRNA_odm, response_gRNA_group_pairs, two_sample_test)
  return(res)
}


estimate_size <- function(df, formula) {
  # second backup: method of moments on Poisson residuals
  backup_2 <- function(pois_fit) {
    MASS::theta.mm(y = df[["expression"]], mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
  }

  # first backup: MLE on poisson residuals
  backup <- function() {
    pois_fit <- stats::glm(formula = formula, data = df, family = stats::poisson())
    gene_precomp_size_out <- tryCatch({
      MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
    }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
    return(gene_precomp_size_out)
  }

  # try to fit a negative binomial GLM with unknown dispersion
  gene_precomp_size_out <- tryCatch({
    fit_nb <- MASS::glm.nb(formula = formula, data = df)
    gene_precomp_size_out <- fit_nb$theta
    return(gene_precomp_size_out)
  }, error = function(e) backup(), warning = function(w) backup())

  return(gene_precomp_size_out)
}
