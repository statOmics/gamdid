#' @title Compare distributions across groups
#'
#' @description
#' Generic function for comparing distributions across levels of a covariate
#' in single-cell proteomics data. Dispatches to the appropriate method based
#' on the class of \code{object}. See
#' \code{\link{compDistr,SummarizedExperiment-method}} and
#' \code{\link{compDistr,QFeatures-method}} for full parameter documentation
#' and examples.
#'
#' @param object A \code{SummarizedExperiment} or \code{QFeatures} instance.
#' @param ... Additional arguments passed to the methods.
#'
#' @return See individual methods for return value details.
#'
#' @seealso \code{\link{visFit}} for visualizing the GAM fits after running
#'   \code{compDistr}.
#'
#' @export
setGeneric(
  "compDistr",
  function(object, ...) standardGeneric("compDistr")
)

#' @title Visualize GAM fits for a single feature
#'
#' @description
#' Generic function for visualizing the GAM fit and differential distribution
#' regions of a single feature after running \code{\link{compDistr}}.
#' Dispatches to the appropriate method based on the class of \code{object}.
#' See \code{\link{visFit,SummarizedExperiment-method}} and
#' \code{\link{visFit,QFeatures-method}} for full parameter documentation and
#' examples.
#'
#' @param object A \code{SummarizedExperiment} or \code{QFeatures} instance
#'   returned by \code{\link{compDistr}}.
#' @param ... Additional arguments passed to the methods.
#'
#' @return A \code{ggplot} object, or a named list of \code{ggplot} objects
#'   (one per pairwise comparison) when more than two groups are present.
#'   The plot is also printed as a side effect.
#'
#' @seealso \code{\link{compDistr}} for performing the distributional analysis
#'   that must precede visualization.
#'
#' @export
setGeneric(
    "visFit",
    function(object, ...) standardGeneric("visFit")
)

#' @title Retrieve top-ranked features after differential distribution analysis
#'
#' @description
#' Generic function for extracting and ranking the most significant features
#' from the output of \code{\link{compDistr}}. Dispatches to the appropriate
#' method based on the class of \code{object}. See
#' \code{\link{topFeatures,SummarizedExperiment-method}} and
#' \code{\link{topFeatures,QFeatures-method}} for full parameter documentation
#' and examples.
#'
#' @param object A \code{SummarizedExperiment} or \code{QFeatures} instance
#'   returned by \code{\link{compDistr}}.
#' @param ... Additional arguments passed to the methods.
#'
#' @return A \code{data.frame} with \code{n} rows sorted by the chosen
#'   p-value column, with features ordered from most to least significant.
#'
#' @seealso \code{\link{compDistr}} which must be run before
#'   \code{topFeatures}. \code{\link{visFit}} for visualizing individual
#'   features.
#'
#' @export
setGeneric(
    "topFeatures",
    function(object, ...) standardGeneric("topFeatures")
)
