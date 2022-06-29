test_that("MAST function works", {
  # read the reported results from the web
  reported_results <- readRDS(url("http://steinmetzlab.embl.de/TAPdata//Figure2_DE_TAP.RDS"))

  # read the Schraivogel data directory from .research_config, which is equivalent to
  # schraivogel_dir <- .get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")
  cmd <- paste0("source ~/.research_config; echo $", "LOCAL_SCHRAIVOGEL_2020_DATA_DIR")
  schraivogel_dir <- system(command = cmd, intern = TRUE)

  # NOTE: Want to use processed but not QC'd data for this check, since Schraivogel's
  # analysis did not include the restriction to one gRNA per cell.

  # read the gRNA ODM
  processed_gRNA_dir <- sprintf("%sprocessed/ground_truth_tapseq/grna_assignment", schraivogel_dir)
  grna_odm_fp <- sprintf("%s/raw_ungrouped.odm", processed_gRNA_dir)
  grna_metadata_fp <- sprintf("%s/raw_ungrouped_metadata.rds", processed_gRNA_dir)
  grna_odm <- ondisc::read_odm(odm_fp = grna_odm_fp, metadata_fp = grna_metadata_fp)

  # read the gene ODM
  processed_gene_dir <- sprintf("%sprocessed/ground_truth_tapseq/gene", schraivogel_dir)
  gene_odm_fp <- sprintf("%s/expression_matrix.odm", processed_gene_dir)
  gene_metadata_fp <- sprintf("%s/metadata.rds", processed_gene_dir)
  gene_odm <- ondisc::read_odm(odm_fp = gene_odm_fp, metadata_fp = gene_metadata_fp)

  # choose a gRNA to test, and pair it with all genes
  guide <- reported_results |>
    dplyr::pull(guide) |>
    unique() |>
    sample(size = 1)
  gene_gRNA_group_pairs <- reported_results |>
    dplyr::select(gene, guide) |>
    dplyr::rename(response_id = gene, gRNA_group = guide) |>
    dplyr::filter(gRNA_group == !!guide) |>
    tibble::as_tibble()

  # modify gRNA ODM to change target of selected gRNA
  new_targets <- grna_odm |>
    ondisc::get_feature_covariates() |>
    tibble::rownames_to_column(var = "guide") |>
    dplyr::mutate(target = ifelse(guide == !!guide, guide, target)) |>
    dplyr::pull(target)
  grna_odm_modified <- grna_odm |>
    ondisc::mutate_feature_covariates(target = new_targets)

  # compute the results based on the schraivogel_method function
  computed_results <- schraivogel_method(
    gene_odm,
    grna_odm_modified,
    gene_gRNA_group_pairs
  )

  # check whether computed and reported p-values match up to roughly machine precision
  expect_lt(
    dplyr::left_join(computed_results |>
      dplyr::rename(p_value_computed = p_value),
    reported_results |>
      dplyr::rename(
        response_id = gene,
        gRNA_group = guide,
        p_value_reported = p_val
      ) |>
      dplyr::select(response_id, gRNA_group, p_value_reported),
    by = c("response_id", "gRNA_group")
    ) |>
      dplyr::summarise(max(abs(p_value_computed - p_value_reported))) |>
      dplyr::pull(),
    1e-10
  )
})
