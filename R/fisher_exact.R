#' Fisher's exact test
#'
#' Carries out Fisher's exact test.
#'
#' @param progress print progress statements?
#' @param control_group the control group to use: "ntc" for cells containing an NTC or "compliment" for the compliment of the treatment cells
#' @inherit abstract_interface
#' @export
fisher_exact <- function(response_odm, grna_odm, response_grna_group_pairs, control_group = "ntc", progress = TRUE) {
  # just set gene threshold to "median"
  threshold <- "median"
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
  res <- abstract_two_sample_test(response_odm = response_odm,
                                  grna_odm = grna_odm,
                                  response_grna_group_pairs = response_grna_group_pairs,
                                  two_sample_test = two_sample_test,
                                  progress = progress,
                                  control_group = control_group)
  return(res)
}


#' Fisher exact (thresholded)
#'
#' Carries out Fisher's exact test. `grna_odm` is assumed to be a gRNA expression matrix. `grna_thresholds`, meanwhile, is an integer vector of gRNA thresholds.
#' For each gRNA threshold (among those in `grna_thresholds`), assigns gRNAs to cells by thresholding the gRNA expression matrix at that level. Then, carries out
#' a fisher exact test. Cells that cannot be confidently assigned to a given gRNA are excluded.
#'
#' @inherit fisher_exact
#' @export
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_perturbseq/gene")
#' grna_odm <- load_dataset_modality("schraivogel/ground_truth_perturbseq/grna_expression")
#' response_grna_group_pairs <-
#'  expand.grid(grna_group = c("CCNE2", "CPQ"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 5))
#' res <- fisher_exact_thresholded(response_odm, grna_odm, response_grna_group_pairs)
#' }
fisher_exact_thresholded <- function(response_odm, grna_odm, response_grna_group_pairs, control_groups = c("ntc", "compliment"), progress = TRUE, grna_thresholds = c(1, 3, 5, 7, 8)) {
  grna_exp_mat <- as.matrix(load_whole_odm(grna_odm))
  res <- lapply(X = control_groups, FUN = function(control_group) {
    lapply(X = grna_thresholds, FUN = function(curr_grna_threshold) {
      # first, work with the grna odm, filtering out cells containing multiple grnas with expression level >= x
      n_grnas_assigned <- Matrix::colSums(grna_exp_mat >= curr_grna_threshold)
      ok_cells <- n_grnas_assigned == 1
      curr_grna_odm <- grna_odm[,ok_cells]
      curr_response_odm <- response_odm[,ok_cells]
      # run the fisher exact test
      fisher_exact(response_odm = curr_response_odm,
                   grna_odm = curr_grna_odm,
                   response_grna_group_pairs = response_grna_group_pairs,
                   control_group = control_group,
                   progress = progress) |>
        dplyr::mutate(grna_threshold = curr_grna_threshold,
                      control_group = control_group)
    }) |> data.table::rbindlist()
  }) |> data.table::rbindlist()
}
