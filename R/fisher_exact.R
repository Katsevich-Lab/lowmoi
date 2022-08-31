#' Fisher's exact test
#'
#' Carries out Fisher's exact test.
#'
#' @inherit abstract_interface
#' @export
fisher_exact <- function(response_odm, grna_odm, response_grna_group_pairs, progress = TRUE) {
  # define the Fisher exact two sample test
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    # create the contingency table
    n_treatment_zero <- sum(target_cells == 0)
    n_treatment_exp <- sum(target_cells >= 1)
    n_control_zero <- sum(control_cells == 0)
    n_control_exp <- sum(control_cells >= 1)

    contingency_table <- matrix(data = c(n_treatment_zero, n_treatment_exp,
                                         n_control_zero, n_control_exp),
                                nrow = 2,
                                byrow = TRUE)
    colnames(contingency_table) <- c("no_expression", "expression")
    rownames(contingency_table) <- c("treatment", "control")
    fit <- fisher.test(contingency_table)
    fit$p.value
  }
  res <- abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress)
  return(res)
}

