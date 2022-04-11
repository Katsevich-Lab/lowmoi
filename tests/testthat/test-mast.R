test_that("MAST function works", {
  schraivogel_dir <- .get_config_path("LOCAL_SCHRAIVOGEL_2020_DATA_DIR")
  TAP.nods_filename <- sprintf("%sraw/ftp/TAP.nods.RDS", schraivogel_dir)
  TAP.per_filename <- sprintf("%sraw/ftp/TAP.per.RDS", schraivogel_dir)
  TAP.nods <- readRDS(TAP.nods_filename)
  TAP.per <- readRDS(TAP.per_filename) |> Matrix::t()

  scramble.cols <- grep("non-targeting", colnames(TAP.per))
  TAP.nods <- TAP.nods[!grepl("CROPseq", rownames(TAP.nods)), ]
  ngenes <- apply(TAP.nods > 0, 2, sum)

  out <- runSeuratTest(
    g = "MYC-B",
    DGE = TAP.nods,
    pert = TAP.per,
    covariate = ngenes,
    scrcols = scramble.cols,
    normfun = function(x) sum(x[x < quantile(x, probs = 0.9)]) + 1
  )

  reported_results <- readRDS(url("http://steinmetzlab.embl.de/TAPdata//Figure2_DE_TAP.RDS"))

  expect_equal(reported_results |>
                 dplyr::filter(guide == "MYC-B") |>
                 dplyr::arrange(gene) |>
                 tibble::as_tibble(),
               out |>
                 dplyr::arrange(gene) |>
                 tibble::as_tibble())
})
