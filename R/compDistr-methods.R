#' Compare distributions
#'
#' This functions compares the distributions between all the levels a covariate
#' using GAMs. To gain power in detecting pure mean shifts, the result is
#' aggregated with a location analysis using msqrob2
#'
#'
#' @param object `SingleCellExperiment` or `QFeatures` instance
#' @param covariate Column in the colData of the object of which the levels
#' need to be compared
#' @param nuisance Formula for covariates you want to correct for
#' @param bins Number of bins in the histogram to which the GAMs are fitted
#' @param adjust if TRUE, return output adjusted for multiple testing
#' @param alpha level of significance
#' @export

setMethod(
  "compDistr", "SummarizedExperiment",
  function(object,
           covariate,
           nuisance = NULL,
           bins = 25,
           adjust = T,
           alpha = 0.05, ...){
    compDistrSCE(sce = object, covariate = covariate, nuisance = nuisance, bins = bins,
                 adjust = adjust, alpha = alpha, ...)
  }
)

#' @param i `character` or `integer` to specify the element of the `QFeatures` that
#'        contains the log expression intensities that will be modelled.
#'
#' @return A `SingleCellExperiment` or a `QFeatures` instance with the GAM objects and comparisons.
#' @export

setMethod(
  "compDistr", "QFeatures",
  function(object,
           i,
           covariate,
           nuisance = NULL,
           bins = 25,
           adjust = T,
           alpha = 0.05, ...){
    sce <- getWithColData(object, i)
    sce_out <- compDistrSCE(sce = sce, covariate = covariate, nuisance = nuisance, bins = bins,
                            adjust = adjust, alpha = alpha, ...)
    object[[i]] <- sce_out
    return(object)
  }
)
