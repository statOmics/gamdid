#' Compare distributions
#'
#' @description
#' This functions compares the distributions between all levels a covariate.
#' The sample space is descretized and Poisson GAMs are fitted for inference.
#' The resulting GAM p-value is aggregated with the with an msqrob2 p value.
#'
#' @include compDistrSCE.R
#'
#' @aliases compDistr compDistr,SummarizedExperiment-method compDistr,QFeatures-method
#'
#' @author Lucas Beerland
#'
#' @param object `SingleCellExperiment` or `QFeatures` instance
#'
#' @param covariate Column in the colData of the object of which the levels
#' need to be compared
#'
#' @param nuisance Formula for covariates you want to correct for
#'
#' @param minObs Minimal number of observation in each group needed to fit GAM
#'
#' @param bins Number of bins in the histogram to which the GAMs are fitted
#'
#' @param adjust if TRUE, return output adjusted for multiple testing
#'
#' @param alpha level of significance
#'
#' @param ... additional arguments passed to the msqrob2() function
#'
#' @export
#'
#' @rdname compDistr
#'
#' @importFrom MultiAssayExperiment getWithColData metadata metadata<-
#' @importFrom methods is

setMethod(
  "compDistr", "SummarizedExperiment",
  function(object,
           covariate,
           nuisance = NULL,
           minObs = 50,
           bins = 25,
           adjust = T,
           alpha = 0.05, ...){
    compDistrSCE(sce = object, covariate = covariate, nuisance = nuisance,
                 minObs = minObs, bins = bins, adjust = adjust, alpha = alpha, ...)
  }
)

#' @param i `character` or `integer` to specify the element of the `QFeatures` that
#'        contains the log expression intensities that will be modelled.
#'
#' @return A `SingleCellExperiment` or a `QFeatures` instance with the GAM objects and comparisons.
#' @export
#' @rdname compDistr

setMethod(
  "compDistr", "QFeatures",
  function(object,
           i,
           covariate,
           nuisance = NULL,
           minObs = 50,
           bins = 25,
           adjust = T,
           alpha = 0.05, ...){
    sce <- getWithColData(object, i)
    sce_out <- compDistrSCE(sce = sce, covariate = covariate, nuisance = nuisance,
                            minObs = minObs, bins = bins,
                            adjust = adjust, alpha = alpha, ...)
    object[[i]] <- sce_out
    return(object)
  }
)
