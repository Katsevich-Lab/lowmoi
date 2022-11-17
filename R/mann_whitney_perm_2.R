#' Mann-whitney test (with permutations)
#'
#' Runs a Mann-Whitney test using permutations.
#'
#' @export
#' @inherit abstract_interface
#' @param B number of permutation replicates
#' @param progress print progress messages?
#' @param full_output return the full output (including the resampled statistics; TRUE) or a reduced output (FALSE)?
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("papalexi/eccite_screen/gene")
#' grna_odm <- load_dataset_modality("papalexi/eccite_screen/grna_assignment")
#' response_grna_group_pairs <-
#'  expand.grid(grna_group = c("CUL3", "CMTM6"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 5))
#' }
mann_whitney_perm <- function(response_odm, grna_odm, response_grna_group_pairs, B = 1000, progress = TRUE, full_output = FALSE) {
  # convert n_rep to integer type (if necessary)
  if (is.character(B)) B <- as.integer(B)
  if (is.character(progress)) progress <- as.logical(progress)
  if (is.character(full_output)) full_output <- as.logical(full_output)

  # obtain the library sizes
  lib_sizes <- get_library_sizes(response_odm)

  # obtain the random indexes
  unique_grna_groups <- as.character(unique(response_grna_group_pairs$grna_group))
  grna_targets <- get_target_assignments_via_max_op(grna_odm)
  grna_group_info <- sceptre2:::get_grna_group_info(grna_group_assignments = grna_targets,
                                                    input_grna_groups = unique_grna_groups)
  random_idxs <- lapply(X = unique_grna_groups, function(unique_grna_group) {
   cbind(matrix(data = seq(1, grna_group_info$n_cells_per_grna[[unique_grna_group]]),
                ncol = 1),
    sceptre2:::get_grna_permutation_idxs(n_cells_per_grna = grna_group_info$n_cells_per_grna,
                                         unique_grna = unique_grna_group,
                                         B = B))
  }) |> stats::setNames(unique_grna_groups)

  # define the permutation test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices, response_id, grna_group) {
    # normalize the data to create samples x (for target cells) and y (for control cells)
    x <- 1000 * target_cells/(lib_sizes[target_cell_indices])
    y <- 1000 * control_cells/(lib_sizes[control_cell_indices])
    n_x <- length(x)
    n_y <- length(y)
    combined <- c(x, y)

    # generate the permutation indices
    synthetic_treatment_indices <- random_idxs[[grna_group]]

    # compute z_null and z_star
    zs <- sapply(X = seq(1, B + 1), FUN = function(i) {
      if (i %% 1000 == 0) print(i)
      col <- synthetic_treatment_indices[,i]
      x_curr <- combined[col]
      y_curr <- combined[-col]
      curr_combined <- c(x_curr, y_curr)
      run_mw_test(x_curr, y_curr)
    })
    z_star <- zs[1]
    z_null <- zs[-1]

    # MW p-value
    lt_p <- stats::pnorm(z_star)
    p_value <- 2 * min(lt_p, 1 - lt_p)

    # empirical p-value
    set.seed(4)
    z_null_jitter <- z_null + runif(n = B, min = -1e-5, max = 1e-5)
    p_emp <- sceptre2:::compute_empirical_p_value(z_star, z_null_jitter, side = "both")

    # R's p-value
    p_r <- wilcox.test(x, y)$p.value

    # ks statistic for N(0,1) fit
    ks_stat <- stats::ks.test(z_null_jitter, stats::pnorm)$statistic[[1]]
    if (full_output) {
      out <- as.data.frame(matrix(data = c(z_star, z_null, p_value, p_emp, ks_stat), nrow = 1))
      colnames(out) <- c("z_star", paste0("z_", seq(1, B)), "p_value", "p_emp", "ks_stat")
    } else {
      out <- as.data.frame(matrix(data = c(z_star, p_value, p_emp, p_r, ks_stat), nrow = 1))
      colnames(out) <- c("z_star", "p_value", "p_emp", "p_r", "ks_stat")
    }
    return(out)
  }

  # run the permutation test
  res <- abstract_two_sample_test(response_odm, grna_odm,
                                  response_grna_group_pairs, two_sample_test,
                                  progress, "ntc", TRUE)
  return(res)
}
