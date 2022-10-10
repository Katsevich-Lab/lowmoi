#' Fisher's exact test
#'
#' Carries out Fisher's exact test.
#'
#' @param progress print progress statements?
#' @param threshold threshold to use to categorize cells into "highly expressed" and "lowly expressed" groups; one of "median" and "zero" (for now)
#' @inherit abstract_interface
#' @export
fisher_exact <- function(response_odm, grna_odm, response_grna_group_pairs, progress = TRUE, threshold = "median") {
  # define the Fisher exact two sample test
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    if (threshold == "median") {
      thresh <- median(c(target_cells, control_cells))
    } else if (threshold == "zero") {
      thresh <- 0
    }
    # create the contingency table
    n_treatment_zero <- sum(target_cells <= thresh)
    n_treatment_exp <- sum(target_cells > thresh)
    n_control_zero <- sum(control_cells <= thresh)
    n_control_exp <- sum(control_cells > thresh)

    contingency_table <- matrix(data = c(n_treatment_zero, n_treatment_exp,
                                         n_control_zero, n_control_exp),
                                nrow = 2,
                                byrow = TRUE)
    colnames(contingency_table) <- c("no_expression", "expression")
    rownames(contingency_table) <- c("treatment", "control")
    fit <- stats::fisher.test(contingency_table)
    fit$p.value
  }
  res <- abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress)
  return(res)
}

