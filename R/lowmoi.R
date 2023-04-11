utils::globalVariables(c("NTC", "Pr(>Chisq)", "ci.hi", "ci.lo", "ci_high", "ci_low", "coef", "component", "contrast",
                         "grna_group", "gene", "response_id", "guide", "logFC", "p_val", "p_value", "primerid", "pvalue",
                         "target_type", "n_nonzero", "n_umis", "dataset", "method", "max_ram", "n_fragments", ".get_config_path", "n_rep",
                         "expressions", "target", "batch", "gem_group", "cell_barcode", "UMI_count", "guide_target", "assigned_grna",
                         "ntc", "bio_rep", "num_ntcs", "est_size"))

#' Abstract method interface
#'
#' Abstract interface that all methods in the `lowmoi` package must satisfy.
#'
#' The `grna_group` column of `response_grna_group_pairs` contains the names of grna groups.
#' These grna groups are assumed to be present in the `target` column of the feature
#' covariate matrix of `grna_odm`. This column should contain entries "non-targeting"
#' indicating the non-targeting grnas. If `target_type` is present as a column of the
#' feature covariate matrix of `grna_odm`, `target_type` is ignored.
#'
#' @param response_odm an expression ODM of responses (typically genes)
#' @param grna_odm an ODM of either grna expressions (counts) or grna
#' assignments (logical)
#' @param response_grna_group_pairs a data frame with columns `response_id` and
#' `grna_group` giving the response ID / grna group pairs to analyze.
#' @name abstract_interface
#'
#' @return a data frame with columns `response_id`, `grna_group`, and `p_value`.
#' @examples
#' \dontrun{
#' # a schraivogel example
#' response_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/gene")
#' grna_odm <- load_dataset_modality("schraivogel/ground_truth_tapseq/grna_assignment")
#' response_grna_group_pairs <-
#'  expand.grid(grna_group = c("OXR1-TSS", "LRRCC1-TSS"),
#'              response_id = sample(ondisc::get_feature_ids(response_odm), 50))
#'
#' # a frangieh example
#' response_odm <- load_dataset_modality("frangieh/control/gene")
#' grna_odm <- load_dataset_modality("frangieh/control/grna_assignment")
#' response_grna_group_pairs <- data.frame(grna_group = "A2M",
#'                                         response_id = sample(ondisc::get_feature_ids(response_odm), 1))
#'
#' # a papalexi example
#' response_odm <- load_dataset_modality("papalexi/eccite_screen/gene")
#' grna_odm <- load_dataset_modality("papalexi/eccite_screen/grna_assignment")
#' response_grna_group_pairs <- data.frame(grna_group = "CUL3",
#'                                         response_id = sample(ondisc::get_feature_ids(response_odm), 4))
#' }
NULL
