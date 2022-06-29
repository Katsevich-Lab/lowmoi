test_that("Weissman function works", {
  # reproduce results from https://github.com/thomasmaxwellnorman/perturbseq_demo/blob/master/perturbseq_demo.ipynb
  reticulate::py_run_file(system.file("python/reproduce_ipynb.py", package = "lowmoi"))

  # get the SCEPTRE2 directory from .research_config, which is equivalent to
  # sceptre2_dir <- .get_config_path("SCEPTRE2_DATA_DIR")
  cmd <- paste0("source ~/.research_config; echo $", "LOCAL_SCEPTRE2_DATA_DIR")
  sceptre2_dir <- system(command = cmd, intern = TRUE)

  # get the directory containing the Weissman data and original results
  weissman_check_dir <- paste0(sceptre2_dir, "results/weissman_check")

  # read the original results
  original_results <- readr::read_csv(paste0(weissman_check_dir, '/original_pvalues.csv'))
  original_results <- original_results |>
    dplyr::rename(gene_id = `...1`) |>
    dplyr::select(-DMSO_control) |>
    tidyr::pivot_longer(-gene_id,
                        names_to = "gRNA_group",
                        values_to = "pvalue") |>
    dplyr::rename(response_id = gene_id) |>
    dplyr::select(gRNA_group, response_id, pvalue)

  # read the Weissman data in CellPopulation format
  thaps_pop <- CellPopulation$from_hdf(paste0(weissman_check_dir, '/thaps_pop.hdf'))

  # convert the Weissman data to ODM format
  odms <- cell_pop_to_odm(thaps_pop)
  gene_odm <- odms$gene_odm
  grna_odm <- odms$grna_odm

  # select a subset of pairs to check
  response_gRNA_group_pairs <- original_results |>
    dplyr::select(gRNA_group, response_id) |>
    dplyr::slice_sample(n = 50)

  # compute the result based on our function
  computed_result <- weissman_method(gene_odm, grna_odm, response_gRNA_group_pairs)

  # check whether computed and original p-values match up to roughly machine precision
  expect_lt(
    computed_result |>
      dplyr::rename(computed_pvalue = p_value)  |>
      dplyr::left_join(original_results |>
                         dplyr::rename(original_pvalue = pvalue),
                       by = c("gRNA_group", "response_id")) |>
      dplyr::summarise(max(abs(computed_pvalue - original_pvalue))) |>
      dplyr::pull(),
    1e-10
  )
})
