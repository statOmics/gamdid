#' @include visFitSCE.R
NULL

#' Visualise the fit of one feature
#'
#' This functions visualizes the histogram and the fitted GAMs used for the
#' distributional analysis. If there is a significant difference in distribution,
#' the shaded regions pinpoint regions that drive the significant difference.
#'
#'
#' @param sce Single cell experiment object
#' @param feature The name of feature (or the number of the row of the feature)
#' that should be visualized
#' @param plot hist or diff, specify the plot you want
#' @param shade if TRUE differential regions are marked with shades
#' @return Interpretable visualisation
#' @export

setMethod(
  "visFit", "SummarizedExperiment",
  function(object,
           feature,
           plot = "hist",
           shade = TRUE,
           ...){
    visFitSCE(sce = object, feature = feature, plot = plot, shade = shade, ...)
  }
)

#' @param i `character` or `integer` to specify the element of the `QFeatures` that
#'        contains the log expression intensities that will be modelled.
#'
#' @return A `SingleCellExperiment` or a `QFeatures` instance with the GAM objects and comparisons.
#' @export

setMethod(
  "visFit", "QFeatures",
  function(object,
           i,
           feature,
           plot = "hist",
           shade = TRUE,
           ...){
    sce <- getWithColData(object, i)
    visFitSCE(sce = sce, feature = feature, plot = plot, shade = shade, ...)
  }
)
