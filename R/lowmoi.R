utils::globalVariables(c("NTC", "Pr(>Chisq)", "ci.hi", "ci.lo", "ci_high", "ci_low", "coef", "component", "contrast",
                         "gRNA_group", "gene", "response_id", "guide", "logFC", "p_val", "p_value", "primerid", "pvalue",
                         "target_type", "n_nonzero", "n_umis", "dataset", "method", "max_ram", "n_fragments", ".get_config_path", "n_rep"))

#' Abstract method interface
#'
#' @param response_odm an expression ODM of responses (typically genes)
#' @param gRNA_odm an ODM of gRNA expressions/indicators
#' @param response_gRNA_group_pairs a data frame with columns `response_id` and `gRNA_group` giving the response ID / gRNA group pairs to analyze
#' @name abstract_interface
#'
#' @return a data frame with columns `response_id`, `gRNA_group`, and `p_value`.
#' @examples
#' \dontrun{
#' gene_odm <- load_dataset_modality("frangieh/co_culture/gene")
#' gRNA_odm <- load_dataset_modality("frangieh/co_culture/grna")
#' response_gRNA_group_pairs <-
#' data.frame(response_id = sample(ondisc::get_feature_ids(gene_odm), 500, FALSE),
#' gRNA_group = sample(ondisc::get_feature_ids(gRNA_odm), 500, FALSE))
#' }
NULL
