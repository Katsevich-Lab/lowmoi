#' Dummy method 1
#'
#' A dummy method; samples uniformly distributed p-values
#'
#' @inherit abstract_interface
#' @param dist a character vector giving the distribution to generate random variates from
#' @export
dummy_method_1 <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, dist = "runif") {
  n_pairs <- nrow(response_gRNA_group_pairs)
  response_gRNA_group_pairs |>
    dplyr::mutate(p_value = do.call(what = dist, args = list(n = n_pairs)))
}


#' Dummy method 2
#'
#' A dummy method; samples beta distributed p-values
#'
#' @inherit abstract_interface
#' @param shape1 first shape of beta distribution; if character, cast to numeric
#' @param shape2 second shape of beta distribution; if character, cast to numeric
#' @export
dummy_method_2 <- function(response_odm, gRNA_odm, response_gRNA_group_pairs, shape1 = 5, shape2 = 10) {
  if (is.character(shape1)) shape1 <- as.numeric(shape1)
  if (is.character(shape2)) shape2 <- as.numeric(shape2)
  n_pairs <- nrow(response_gRNA_group_pairs)
  response_gRNA_group_pairs |>
    dplyr::mutate(p_value = stats::rbeta(n = n_pairs, shape1 = shape1, shape2 = shape2))
}
