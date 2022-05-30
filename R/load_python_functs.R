.onLoad <- function(libname, pkgname) {
  py_modules_path <- system.file("python", package = "lowmoi")

  # MIMOSCA
  mimosca_module <- reticulate::import_from_path(module = "mimosca",
                                                 path = py_modules_path,
                                                 delay_load = TRUE)
  run_mimosca <<- mimosca_module$run_mimosca
  get_dense_array <<- mimosca_module$get_dense_array

  # Weissman (cell population class)
  cell_pop_module <- reticulate::import_from_path(module = "cell_population",
                                                  path = py_modules_path,
                                                  delay_load = TRUE)
  CellPopulation <<- cell_pop_module$CellPopulation

  # Weissman (expression normalization)
  exp_norm_module <- reticulate::import_from_path(module = "expression_normalization",
                                                  path = py_modules_path,
                                                  delay_load = TRUE)
  normalize_to_gemgroup_control <<- exp_norm_module$normalize_to_gemgroup_control

  # Weissman (differential expression)
  diff_exp_module <- reticulate::import_from_path(module = "differential_expression",
                                                  path = py_modules_path,
                                                  delay_load = TRUE)
  ks_de <<- diff_exp_module$ks_de
}
