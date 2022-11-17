
#' Abstract two sample test
#'
#' An abstract a two-sample test. Pass function `two_sample_test` to carry out a given two-sample test (e.g., a t-test or a permutation test).
#'
#' @inherit abstract_interface
#' @param two_sample_test a two-sample test; should take as arguments (i) vector of expressions of target cells, (ii) vector of expressions of control cells, (iii) the indices of cells receiving the targeting grna, and (iv) the indices of the cells receiving the NT grnas.
#' @param progress print progress messages?
#' @param control_group the control group to use: "ntc" for cells containing an NTC or "compliment" for the compliment of the treatment cells
#' @param cbind_res cbind the result to `response_grna_group_pairs` (TRUE) or simply return the `response_grna_group_pairs` data frame with the vector of p-values appended (FALSE)?
#' @export
#' @examples
#' \dontrun{
#' two_sample_test <- function(target_cells, control_cells, response_id,
#' target_cell_indices, control_cell_indices) {
#' t.test(target_cells, control_cells)$p.value
#' }
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' grna_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' response_grna_group_pairs <-
#'  expand.grid(grna_group = c("CCNE2", "HS2"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 50))
#' abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, TRUE)
#' }
abstract_two_sample_test <- function(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress, control_group = "ntc", cbind_res = FALSE) {
  set.seed(4)
  # get grna assignments and target assignments; obtain indices of NT cells
  grna_targets <- get_target_assignments_via_max_op(grna_odm)

  # loop through the pairs, calculating a p-value for each
  res <- apply(X = response_grna_group_pairs, MARGIN = 1, FUN = function(r) {
    grna_group <- as.character(r[["grna_group"]])
    target_cell_indices <- grna_targets == grna_group
    response_id <- as.character(r[["response_id"]])
    if (progress) print(paste0("Analyzing ", response_id, " and ", grna_group))
    if (control_group == "ntc") {
      control_cell_indices <- grna_targets == "non-targeting"
    } else if (control_group == "compliment") {
      control_cell_indices <- !target_cell_indices
    } else {
      stop("Control group not recognized.")
    }
    # get the target and control cells (NOTE: perhaps only load if necessary)
    target_cells <- response_odm[[response_id, target_cell_indices]] |> as.numeric()
    control_cells <- response_odm[[response_id, control_cell_indices]] |> as.numeric()
    two_sample_test(target_cells, control_cells, target_cell_indices,
                    control_cell_indices, response_id, grna_group)
  }, simplify = FALSE)
  if (cbind_res) {
    to_attach <- data.table::rbindlist(res)
    response_grna_group_pairs <- cbind(response_grna_group_pairs, to_attach)
  } else {
    response_grna_group_pairs$p_value <- unlist(res)
  }
  return(response_grna_group_pairs)
}
