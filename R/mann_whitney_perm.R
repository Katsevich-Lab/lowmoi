#' Mann-whitney test (with permutations)
#'
#' Runs a Mann-Whitney test using permutations.
#'
#' @inherit abstract_interface
#' @param B number of permutation replicates
#' @param progress print progress messages?
mann_whitney_perm <- function(response_odm, grna_odm, response_grna_group_pairs, B = 10, progress = TRUE) {
  # convert n_rep to integer type (if necessary)
  if (is.character(B)) B <- as.integer(B)
  if (is.character(progress)) progress <- as.logical(progress)

  # obtain the library sizes
  lib_sizes <- get_library_sizes(response_odm)

  # define the permutation test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    # compute the relative expressions
    expressions <- c(target_cells, control_cells)
    lib_sizes <- c(lib_sizes[target_cell_indices], lib_sizes[control_cell_indices])
    rm(target_cell_indices, control_cell_indices)
    rel_expressions <- log(100000 * expressions/lib_sizes + 1)
    rm(expressions, lib_sizes)

    # obtain the relevant sample sizes
    n_cells_curr_de <- length(rel_expressions)
    n_cells_treat <- length(target_cells)
    n_cells_control <- n_cells_curr_de - n_cells_treat

    # sample the synthetic treatment idxs
    synthetic_treatment_indices <- replicate(n = B, expr = sample.int(n = n_cells_curr_de, size = n_cells_treat))
    ground_truth_treatment_idxs <- seq(1, length(target_cells))

    # compute the permuted test statistics
    permutation_runs <- run_permutations_mw(rel_expressions, ground_truth_treatment_idxs, synthetic_treatment_indices)

    # compute the normalized statistics
    m_u <- n_cells_treat * n_cells_control / 2
    sigma_u <- sqrt(n_cells_treat * n_cells_control * (n_cells_treat +  n_cells_control + 1)/12)
    z <- (permutation_runs - m_u)/sigma_u
    p_emp <- sceptre2:::compute_empirical_p_value(z_star = z[1], z_null = z[-1], side = "both")

    out <- matrix(data = c(z, m_u, sigma_u, p_emp), nrow = 1)
    colnames(out) <- c("z_null", paste0("z_", seq(1, B)), "mu", "sigma", "p_emp")
    return(as.data.frame(out))
  }

  # run the permutation test
  res <- abstract_two_sample_test(response_odm, grna_odm, response_grna_group_pairs, two_sample_test, progress, TRUE)
  return(res)
}


run_permutations_mw <- function(rel_expressions, ground_truth_treatment_idxs, synthetic_treatment_indices) {
  u_star <- compute_mw_test_stat(rel_expressions, ground_truth_treatment_idxs)
  u_null <- apply(X = synthetic_treatment_indices, MARGIN = 2, FUN = function(curr_truth_treatment_idxs) {
    cat("*")
    compute_mw_test_stat(rel_expressions, curr_truth_treatment_idxs)
  })
  cat("\n")
  return(c(u_star, u_null))
}


compute_mw_test_stat <- function(rel_expressions, curr_truth_treatment_idxs) {
  y <- rel_expressions[curr_truth_treatment_idxs]
  x <- rel_expressions[-curr_truth_treatment_idxs]
  m <- outer(X = x, Y = y, FUN = function(x, y) {
    ifelse(x == y, 0.5, as.integer(x > y))
  })
  u <- sum(m)
}
