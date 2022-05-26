#' Dummy method 1
#'
#' A dummy method; samples uniformly distributed p-values
#'
#' @inherit abstract_interface
#' @export
dummy_method_1 <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  n_pairs <- nrow(response_gRNA_group_pairs)
  response_gRNA_group_pairs |>
    dplyr::mutate(p_value = stats::runif(n_pairs))
}


#' Dummy method 2
#'
#' A dummy method; samples beta distributed p-values
#'
#' @inherit abstract_interface
#' @export
dummy_method_2 <- function(response_odm, gRNA_odm, response_gRNA_group_pairs) {
  n_pairs <- nrow(response_gRNA_group_pairs)
  response_gRNA_group_pairs |>
    dplyr::mutate(p_value = stats::rbeta(n_pairs, 1, 0.5))
}
