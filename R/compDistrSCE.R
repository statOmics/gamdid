compDistrSCE <- function(sce, covariate, nuisance = NULL, bins = 25,
                      adjust = T, alpha = 0.05, ...){
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
    sce_partres <- partialResiduals(
      sce = sce,
      covariate = covariate,
      msqrob_models = location$msqrob_models,
      msqrob_formula = location$msqrob_formula
      )
    }

  shape <- shapeTest(
    sce = sce,
    covariate = covariate,
    bins = bins
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
  if (!is.null(nuisance)) assays(sce_out)$assay_partialresiduals <- assay(sce)
  rowData(sce_out)$GAMmodels <- shape$models
  rowData(sce_out) <- cbind(rowData(sce_out), out$pvals)
  metadata(sce_out)$adj_alpha <- out$adj_alpha

  sce_out
}
