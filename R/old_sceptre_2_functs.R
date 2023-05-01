

get_grna_group_info <- function(grna_group_assignments, input_grna_groups) {
  unique_grnas <- c(input_grna_groups |> stats::setNames(input_grna_groups),
                    "non-targeting" = "non-targeting")
  # get the indices of each gRNA
  grna_specific_idxs <- lapply(unique_grnas, function(unique_grna) {
    which(grna_group_assignments == unique_grna)
  })
  # get the number of cells per gRNA; also record the number of NT cells
  n_cells_per_grna <- table(grna_group_assignments)[unique_grnas]
  return(list(grna_specific_idxs = grna_specific_idxs,
              n_cells_per_grna = n_cells_per_grna))
}


get_grna_permutation_idxs <- function(n_cells_per_grna, unique_grna, B) {
  n_nt_cells <- n_cells_per_grna[["non-targeting"]]
  n_cells_curr_grna_group <- n_cells_per_grna[[unique_grna]]
  n_cells_curr_de <- n_cells_curr_grna_group + n_nt_cells
  synthetic_treatment_idxs <- replicate(n = B, expr = sample.int(n = n_cells_curr_de,
                                                                 size = n_cells_curr_grna_group))
}


compute_empirical_p_value <- function(z_star, z_null, side) {
  out_p <- switch(side,
                  "left" = mean(c(-Inf, z_null) <= z_star),
                  "right" = mean(c(Inf, z_null) > z_star),
                  "both" = 2 * min(mean(c(-Inf, z_null) <= z_star),
                                   mean(c(Inf, z_null) > z_star)))
  return(out_p)
}


run_glm_perm_score_test_with_ingredients <- function(Z, working_resid, w, index_mat) {
  # compute Z^T w (to be used throughout)
  ZtW <- sapply(X = seq(1, length(w)), FUN = function(i) w[i] * Z[i,])

  # compute the precision matrix P = Z^t W Z
  P <- ZtW %*% Z

  # compute the spectral decomposition of P
  P_decomp <- eigen(P)

  # obtain U and Lambda^(-1/2)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1/sqrt(P_decomp$values)

  # compute the matrix B = Lambda^(1/2) U^t (Z^t W)
  B <- (Lambda_minus_half * t(U)) %*% ZtW

  # next, compute the vector W M (Y - mu_hat)
  a <- w * working_resid

  # compute the vector of z-scores
  z_scores <- low_level_score_test_vectorized(a = a, B = B, w = w, index_mat = index_mat)
  return(z_scores)
}

