#' Mann-whitney test (with permutations)
#'
#' Runs a Mann-Whitney test using permutations.
#'
#' @export
#' @inherit abstract_interface
#' @param B number of permutation replicates
#' @param progress print progress messages?
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("papalexi/eccite_screen/gene")
#' grna_odm <- load_dataset_modality("papalexi/eccite_screen/grna_assignment")
#' response_grna_group_pairs <-
#'  expand.grid(grna_group = c("CUL3", "CMTM6"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 5))
#' }
mann_whitney_perm <- function(response_odm, grna_odm, response_grna_group_pairs, B = 100000, progress = FALSE) {
  # convert n_rep to integer type (if necessary)
  if (is.character(B)) B <- as.integer(B)
  if (is.character(progress)) progress <- as.logical(progress)

  # obtain the library sizes
  lib_sizes <- get_library_sizes(response_odm)

  # define the permutation test function
  two_sample_test <- function(target_cells, control_cells, target_cell_indices, control_cell_indices, response_id, grna_group) {
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


#' Mann-whitney test (using MIMOSCA-type strategy to construct the null distribution)
#'
#' Runs the Mann-Whitney test, using a MIMOSCA-style strategy to construct the null distribution (i.e., pooling information across genes).
#'
#' @export
#' @inherit abstract_interface
#' @param B number of permutation replicates
#' @param control_group group of cells to use for control (either "ntc" or "compliment")
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("papalexi/eccite_screen/gene")
#' grna_odm <- load_dataset_modality("papalexi/eccite_screen/grna_assignment")
#' response_grna_group_pairs <-
#'  expand.grid(grna_group = c("CUL3", "CMTM6"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 5))
#' }
mann_whitney_perm_mimosca <- function(response_odm, grna_odm, response_grna_group_pairs, control_group = "ntc", B = 50) {
  exp_data <- load_whole_odm(response_odm)
  # obtain the library sizes
  lib_sizes <- get_library_sizes(response_odm)
  # obtain normalized exp matrix
  exp_data_norm <- Matrix::t(Matrix::t(exp_data)/lib_sizes)
  rm(exp_data)
  # assign grnas to cells; obtain the indexes of the negative control grnas
  grna_targets <- get_target_assignments_via_max_op(grna_odm)

  # cycle over unique gRNA groups
  unique_grna_grps <- as.character(unique(response_grna_group_pairs$grna_group))
  out <- lapply(X = unique_grna_grps, FUN = function(grna_grp) {
    target_cell_indices <- grna_targets == grna_grp
    if (control_group == "ntc") {
      control_cell_indices <- grna_targets == "non-targeting"
    } else if (control_group == "compliment") {
      control_cell_indices <- !target_cell_indices
    } else {
      stop("Control group not recognized.")
    }

    # create the "subset" mat
    target_cell_mat <- exp_data_norm[,target_cell_indices]
    control_cell_mat <- exp_data_norm[,control_cell_indices]
    n_target_cells <- ncol(target_cell_mat)
    n_control_cells <- ncol(control_cell_mat)
    n_cells <- n_target_cells + n_control_cells
    subset_mat <- cbind(target_cell_mat, control_cell_mat)
    rm(target_cell_mat, control_cell_mat)

    # generate the matrix of permutation indices
    synthetic_treatment_indices <- replicate(n = B,
                                             expr = sample.int(n = n_target_cells + n_control_cells,
                                                               size = n_target_cells,
                                                               replace = FALSE))

    # compute the null distribution in a double for loop: outer loop over permutations
    null_dist <- lapply(X = seq(1, B), FUN = function(i) {
      print(paste0("Working on permutation ", i, " of ", B))
      curr_treatment_idxs <- synthetic_treatment_indices[,i]
      # inner loop over genes
      gene_wise_test_stats <- apply(subset_mat, 1, function(r) {
        xs <- r[curr_treatment_idxs]
        ys <- r[-curr_treatment_idxs]
        run_mw_test(xs, ys)
      })
    }) |> unlist() |> setNames(NULL)

    # compute the gene-wise "ground truth" test statistics
    print("Working on ground truth statistics.")
    ground_truth_stats <- apply(X = subset_mat, MARGIN = 1, FUN = function(r) {
      xs <- r[seq(1, n_target_cells)]
      ys <- r[-seq(1, n_target_cells)]
      run_mw_test(xs, ys)
    })

    # finally, compute the gene-wise p-values
    p_vals <- sapply(X = ground_truth_stats, FUN = function(test_stat) {
      sceptre2:::compute_empirical_p_value(test_stat, null_dist, "both")
    })
    df <- data.frame(grna_group = grna_grp,
                     response_id = names(p_vals),
                     p_value = p_vals)
    row.names(df) <- NULL
    return(df)
  }) |> dplyr::bind_rows()

  # join the p-values with the input data frame
  ret <- dplyr::left_join(x = response_grna_group_pairs, y = out,
                          by = c("grna_group", "response_id"))
  return(ret)
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
