#' NB regression
#'
#' Implements a negative binomial regression.
#'
#' @inherit abstract_interface
#' @param progress print progress messages?
#' @export
nb_regression <- function(response_odm, grna_odm, response_grna_group_pairs, progress = TRUE, with_covariates = TRUE) {
  if (is.character(progress)) progress <- as.logical(progress)
  # obtain cell covariate data frame
  cell_covariate_df <- response_odm |> ondisc::get_cell_covariates()

  # get the formula
  if (with_covariates) {
    my_formula_str <- response_odm@misc$nb_regression_formula
    my_formula <- stats::as.formula(paste0("expression ", my_formula_str))
  } else {
    my_formula <- stats::formula(expression ~ log(n_umis))
  }

  # define the NB regression test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices, response_id, grna_group) {
    # construct the data matrix to pass to GLM
    df <- create_design_matrix(target_cells, control_cells, target_cell_indices, control_cell_indices, cell_covariate_df)
    # first, use aux function to estimate size
    est_size <- max(estimate_size(df, my_formula), 0.01)

    # fit GLM
    p_val <- tryCatch({
      fit_nb <- stats::glm(formula = my_formula, family = MASS::negative.binomial(est_size),  data = df)
      z <- statmod::glm.scoretest(fit_nb, df$pert_indicator)
      2 * pnorm(-abs(z), lower.tail = TRUE)
    }, error = function(e) 1, warning = function(e) 1)

    return(p_val)
  }

  # run the NB regression on all the data
  res <- abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress)
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


create_design_matrix <- function(target_cells, control_cells, target_cell_indices, control_cell_indices, cell_covariate_df) {
  df_left <- data.frame(expression = c(target_cells, control_cells),
                        pert_indicator = c(rep(1, length(target_cells)),
                                           rep(0, length(control_cells))))
  df_right <- rbind(cell_covariate_df[target_cell_indices,], cell_covariate_df[control_cell_indices,])
  df <- cbind(df_left, df_right)
  return(df)
}


# helper functions: NB regression with and without covariates
#' @export
nb_regression_w_covariates <- function(response_odm, grna_odm, response_grna_group_pairs, progress = TRUE) {
  nb_regression(response_odm = response_odm,
                grna_odm = grna_odm,
                response_grna_group_pairs = response_grna_group_pairs,
                progress = progress,
                with_covariates = TRUE)
}


#' @export
nb_regression_no_covariates <- function(response_odm, grna_odm, response_grna_group_pairs, progress = TRUE) {
  nb_regression(response_odm = response_odm,
                grna_odm =grna_odm,
                response_grna_group_pairs = response_grna_group_pairs,
                progress = progress,
                with_covariates = FALSE)
}
