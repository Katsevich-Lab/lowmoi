.onLoad <- function(libname, pkgname) {
  py_modules_path <- system.file("python", package = "lowmoi")

  # MIMOSCA
  mimosca_module <- reticulate::import_from_path(module = "mimosca",
                                                 path = py_modules_path,
                                                 delay_load = TRUE)
  assign('run_mimosca', mimosca_module$run_mimosca, envir = topenv())

  # Weissman (cell population class)
  cell_pop_module <- reticulate::import_from_path(module = "cell_population",
                                                  path = py_modules_path,
                                                  delay_load = TRUE)
  assign('CellPopulation', cell_pop_module$CellPopulation, envir = topenv())

  # Weissman (expression normalization)
  exp_norm_module <- reticulate::import_from_path(module = "expression_normalization",
                                                  path = py_modules_path,
                                                  delay_load = TRUE)
  assign('normalize_to_gemgroup_control', exp_norm_module$normalize_to_gemgroup_control, envir = topenv())

  # Weissman (differential expression)
  diff_exp_module <- reticulate::import_from_path(module = "differential_expression",
                                                  path = py_modules_path,
                                                  delay_load = TRUE)
  assign('ks_de', diff_exp_module$ks_de, envir = topenv())
}
