#' Mann-whitney test (with permutations)
#'
#' Runs a Mann-Whitney test using permutations.
#'
#' @export
#' @inherit abstract_interface
#' @param B number of permutation replicates
#' @param progress print progress messages?
mann_whitney_perm <- function(response_odm, grna_odm, response_grna_group_pairs, B = 100000, progress = FALSE) {
  # convert n_rep to integer type (if necessary)
  if (is.character(B)) B <- as.integer(B)
  if (is.character(progress)) progress <- as.logical(progress)

  # obtain the library sizes
  lib_sizes <- get_library_sizes(response_odm)

  # define the permutation test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices) {
    # normalize the data to create samples x (for target cells) and y (for control cells)
    x <- target_cells/(lib_sizes[target_cell_indices])
    y <- control_cells/(lib_sizes[control_cell_indices])
    combined <- c(x, y)

    # generate the permutation indices
    synthetic_treatment_indices <- replicate(n = B,
                                             expr = sample.int(n = length(combined), size = length(x), replace = FALSE))

    # compute the empirical null distribution
    z_null <- sapply(X = seq(1, B), FUN = function(i) {
      if (i %% 1000 == 0) print(i)
      col <- synthetic_treatment_indices[,i]
      x_curr <- combined[col]
      y_curr <- combined[-col]
      run_mw_test(x_curr, y_curr)
    })
    z_null_jitter <- z_null + runif(n = B, min = -1e-5, max = 1e-5)
    z_star <- run_mw_test(x, y)

    # MW p-value
    p_value <- 2 * min(pnorm(z_star), pnorm(z_star, lower.tail = FALSE))

    # empirical p-value
    p_emp <- sceptre2:::compute_empirical_p_value(z_star, z_null_jitter, side = "both")

    # ks statistic for N(0,1) fit
    ks_stat <- stats::ks.test(z_null_jitter, pnorm)$statistic[[1]]
    out <- as.data.frame(matrix(data = c(z_star, z_null, p_value, p_emp, ks_stat), nrow = 1))
    colnames(out) <- c("z_star", paste0("z_", seq(1, B)), "p_value", "p_emp", "ks_stat")

    return(out)
  }

  # run the permutation test
  res <- abstract_two_sample_test(response_odm, grna_odm,
                                  response_grna_group_pairs, two_sample_test,
                                  progress, "ntc", TRUE)
  return(res)
}


run_mw_test <- function(x, y) {
  r <- c(x, y)
  r <- rank(r)
  n.x <- as.double(length(x))
  n.y <- as.double(length(y))
  STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1)/2)
  NTIES <- table(r)
  z <- STATISTIC - n.x * n.y/2
  SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) -  sum(NTIES^3 - NTIES)/((n.x + n.y) * (n.x + n.y - 1))))
  CORRECTION <- sign(z) * 0.5
  z <- (z - CORRECTION)/SIGMA
  return(z[[1]])
}
