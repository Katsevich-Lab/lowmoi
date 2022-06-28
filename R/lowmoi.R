utils::globalVariables(c("NTC", "Pr(>Chisq)", "ci.hi", "ci.lo", "ci_high", "ci_low", "coef", "component", "contrast",
                         "gRNA_group", "gene", "response_id", "guide", "logFC", "p_val", "p_value", "primerid", "pvalue",
                         "target_type", "n_nonzero", "n_umis", "dataset", "method", "max_ram", "n_fragments", ".get_config_path", "n_rep"))

#' Abstract method interface
#'
#' Abstract interface that all methods in the `lowmoi` package must satisfy.
#'
#' The `gRNA_group` column of `response_gRNA_group_pairs` contains the names of gRNA groups.
#' These gRNA groups are assumed to be present in the `target` column of the feature
#' covariate matrix of `gRNA_odm`. This column should contain entries "non-targeting"
#' indicating the non-targeting gRNAs. If `target_type` is present as a column of the
#' feature covariate matrix of `gRNA_odm`, `target_type` is ignored.
#'
#' @param response_odm an expression ODM of responses (typically genes)
#' @param gRNA_odm an ODM of either gRNA expressions (counts) or gRNA
#' assignments (logical)
#' @param response_gRNA_group_pairs a data frame with columns `response_id` and
#' `gRNA_group` giving the response ID / gRNA group pairs to analyze.
#' @name abstract_interface
#'
#' @return a data frame with columns `response_id`, `gRNA_group`, and `p_value`.
#' @examples
#' \dontrun{
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' gRNA_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' response_gRNA_group_pairs <-
#'  expand.grid(gRNA_group = c("CCNE2-TSS", "HS2-enh"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 50))
#' }
NULL
