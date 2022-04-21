#' MAST DE function
#'
#' Perform differential expression analysis between perturbed and control cells using MAST. This code
#' was drawn from https://github.com/argschwind/TAPseq_manuscript/blob/master/vignettes/functions/MAST.R
#' and only lightly edited.
#'
#' The first column in colData is expected to provide the perturbation status as a factor with two
#' levels, where the first level is non-perturbed cells. If this is not the case, the function will
#' attempt to convert the first column to a two level factor, with any risks that entails.
#' Any other columns in colData will be used as covariates during modelling, but only the
#' perturbation will be tested for signficance.
#'
#' @param pert_object SingleCellExperiment object containing gene expression data and cell groupings
#'   in colData. The perturbation to be tested is assumed to be the first column in colData!
de_MAST <- function(pert_object) {

  # normalize data using censored mean normalization with log1p as scale function
  counts <- SummarizedExperiment::assay(pert_object, "counts")
  logcounts <- normalize_cens_mean(counts, percentile = 0.9, norm_counts = FALSE, scale_fun = log1p)
  SummarizedExperiment::assay(pert_object, "logcounts") <- logcounts

  # add some row- and colData expected by MAST (not used in DE tests)
  SummarizedExperiment::rowData(pert_object) <- cbind(SummarizedExperiment::rowData(pert_object), primerid = rownames(pert_object))
  SummarizedExperiment::colData(pert_object) <- cbind(SummarizedExperiment::colData(pert_object), wellKey = colnames(pert_object))

  # coerce pert_object to SingleCellAssay, since MAST requires that as input
  sca <- MAST::SceToSingleCellAssay(pert_object)

  # prepare colData for model fitting
  colnames(SummarizedExperiment::colData(sca))[1] <- "pert"  # make sure 1st column is called pert for coeff testing
  pert <- as.factor(SummarizedExperiment::colData(sca)$pert)
  levels(pert) <- c(0, 1)
  SummarizedExperiment::colData(sca)$pert <- pert

  # fit hurdle model
  vars <- colnames(SummarizedExperiment::colData(sca))
  vars_to_test <- vars[vars != "wellKey"]
  model <- paste0("~", paste(vars_to_test, collapse = "+"))
  zlm_fit <- MAST::zlm(stats::as.formula(model), sca)

  # following command throws warning with MAST v1.8.1, so lrt test is performed separately
  summary_zlm_fit <- MAST::summary(zlm_fit, logFC = TRUE, doLRT = "pert1")

  pvalue_info <- summary_zlm_fit$datatable |>
    tibble::as_tibble() |>
    dplyr::filter(contrast == "pert1", component == "H") |>
    dplyr::select(primerid, `Pr(>Chisq)`)

  lfc_info <- summary_zlm_fit$datatable |>
    tibble::as_tibble() |>
    dplyr::filter(contrast == "pert1", component == "logFC") |>
    dplyr::select(primerid, ci.hi, ci.lo, coef)

  dplyr::left_join(pvalue_info,
                             lfc_info,
                             by = "primerid") |>
    dplyr::rename(gene = primerid, ci_high = ci.hi, ci_low = ci.lo, logFC = coef, pvalue = `Pr(>Chisq)`) |>
    dplyr::select(gene, logFC, ci_high, ci_low, pvalue)
}


#' Censored mean normalization
#'
#' Normalize DGE data based on censored mean, which excludes genes highly expressed genes when
#' calculating size factors. This code
#' was drawn from https://github.com/argschwind/TAPseq_manuscript/blob/master/vignettes/functions/MAST.R
#' and only lightly edited.
#'
#'
#' @param dge Matrix or data.frame containing digital gene expression data for all cells and genes
#'   to be tested. Rows are genes and columns are cells, with row names being gene ids and column
#'   names cell ids.
#' @param percentile Genes of each cell with expression from that percentile upwards will be
#'   excluded to calculate size factors.
#' @param norm_counts (logical) Return normalized counts, meaning normalized data is rounded to the
#'   nearest integer to resemble a raw count matrix required by certain DE methods such as DEsingle.
#'   Default = FALSE.
#' @param scale_fun Scaling function to be applied to censored mean normalized data. Ignored if
#'   norm_counts = TRUE.
normalize_cens_mean <- function(dge, percentile = 0.9, norm_counts = FALSE, scale_fun = log1p) {

  # function to calculate censored size factor for one cell
  calculate_sf <- function(x) {
    sum(x[x <= stats::quantile(x, probs = percentile)]) + 1
  }

  # calculate size factors for normalization
  size_factors <- apply(dge, MARGIN = 2, FUN = calculate_sf)

  # normalize data of each cell based on computed size factors
  dge_norm <- t(t(dge) / size_factors) * mean(size_factors)

  # round to nearest number to generate normalized counts or scale using specified function
  if (norm_counts == TRUE) {
    round(dge_norm)
  }else if (!is.null(scale_fun)) {
    scale_fun(dge_norm)
  }
}

# This code was drawn from https://github.com/argschwind/TAPseq_manuscript/blob/master/vignettes/functions/functions_runSeuratTest.R
# and only lightly edited.
runSeuratTest <- function(g = "DS42_eGATA1_D", DGE, pert, scrcols, scalefun= log1p, normfun = function(x) 1,normcounts=T, covariate) {

  is.scramble <- apply(pert[,scrcols],1,function(x) any(x > 0)) & apply(pert[,-scrcols],1,function(x) all(x == 0)) # <<- CHECK!
  is.genes <- !grepl("CROPseq", rownames(DGE))

  object <- Seurat::CreateSeuratObject(DGE[is.genes,]) #normalization causes problems, only log-transform
  sf <-  apply(as.matrix(DGE[is.genes,]) , 2, normfun)
  object <- Seurat::SetAssayData(object, new.data = scalefun(t(t(as.matrix(DGE[is.genes,])) /sf )))
  if (normcounts) object <- Seurat::SetAssayData(object, new.data = round(t(t(as.matrix(DGE[is.genes,])) /sf ) * mean(sf)), slot = "counts")  else object <- Seurat::SetAssayData(object, new.data = as.matrix(DGE[is.genes,]), slot = "counts")

  #is.scramble <- apply(pert[,1:4],1,function(x) any(x > 0)) & apply(pert[,5:ncol(pert)],1,function(x) all(x == 0)) # <<- CHECK! and remove
  Seurat::Idents(object) <- ifelse(pert[,g] == 1, g, ifelse(is.scramble, "scramble", "other"))
  object$covariate <- covariate
  a <- colnames(object)[Seurat::Idents(object) == g]
  b <- colnames(object)[Seurat::Idents(object) == "scramble"]

  use.a <- a
  use.b <- b
  testobject <- subset(object, cells = c(use.a, use.b))

  a.obj <- subset(object, cells = use.a)
  b.obj <- subset(object, cells = use.b)

  rowmeans.a <- apply(Seurat::GetAssayData(a.obj, slot="counts"),1,mean)
  rowmeans.b <- apply(Seurat::GetAssayData(b.obj, slot="counts"),1,mean)
  abs.diff <- rowmeans.a - rowmeans.b

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = Seurat::GetAssayData(testobject, slot = "counts")  ) , colData = data.frame(row.names = colnames(testobject), ident = Seurat::Idents(testobject), cov = testobject$covariate))
  MAST.result <- de_MAST(sce)
  m <- data.frame(row.names = MAST.result$gene, p_val = MAST.result$pvalue, p_val_adj = stats::p.adjust( MAST.result$pvalue, method = "fdr"), avg_logFC = MAST.result$logFC)

  #m <- FindMarkers(testobject, ident.1 = g, ident.2 = "scramble", test.use = "MAST", logfc.threshold =0 , min.pcr = 0)
  m$guide <- g
  m$gene <- rownames(m)
  m$ncells <- length(a)
  m$absdiff <- abs.diff[as.character(m$gene)]
  m
}
