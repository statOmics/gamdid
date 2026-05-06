#' @importFrom SummarizedExperiment assays assays<- assay rowData rowData<-

#' @include locationTest.R partialResiduals.R shapeTest.R aggregatePvals.R adjustOutput.R

compDistrSCE <- function(sce, covariate, nuisance = NULL, minObs = 50,
                         bins = 25, adjust = T, alpha = 0.05, ...){
  if (any(grepl(paste0("\\b",covariate, "\\b"),as.character(nuisance)))) stop("covariate of interest cannot be in nuisance formula")

  sce_original <- sce
  # ------------------------------------------------------------
  # Location analysis
  # ------------------------------------------------------------
  location <- locationTest(
    sce = sce,
    covariate = covariate,
    nuisance = nuisance,
    ...
  )

  # ------------------------------------------------------------
  # Shape analysis (on partial residuals)
  # ------------------------------------------------------------
  if (!is.null(nuisance)) {
    sce <- partialResiduals(
      sce = sce,
      covariate = covariate,
      msqrob_models = location$msqrob_models,
      msqrob_formula = location$msqrob_formula
      )
    }

  shape <- shapeTest(
    sce = sce,
    covariate = covariate,
    bins = bins,
    minObs = minObs
  )

  # ------------------------------------------------------------
  # Aggregate p-values
  # ------------------------------------------------------------
  pval <- aggregatePvals(
    pval_shape = shape$pvals,
    pval_location = location$pvals
  )

  # ------------------------------------------------------------
  # Output formatting
  # ------------------------------------------------------------
  out <- adjustOutput(
    pval = pval,
    pval_shape = shape$pvals,
    pval_location = location$pvals,
    adjust = adjust,
    alpha = alpha
  )
  # ------------------------------------------------------------
  # Assemble SCE
  # ------------------------------------------------------------
  sce_out <- sce_original
  if (!is.null(nuisance)) assays(sce_out)$partialResiduals <- assay(sce)
  rowData(sce_out) <- cbind(rowData(sce_out), out$pvals)
  metadata(sce_out)$covariate <- covariate
  metadata(sce_out)$comparisons <- setdiff(colnames(location$pvals), "omnibus")
  metadata(sce_out)$bins <- bins
  metadata(sce_out)$adj_alpha <- out$adj_alpha
  sce_out
}
