#' Interpretable visualization of the GAM fit of one feature
#'
#' @description This function visualizes the histograms and the fitted GAMs used
#'              for inference. Shaded regions indicate differential regions in the distribution.
#'
#' @include visFitSCE.R
#'
#' @aliases visFit visFit,SummarizedExperiment-method visFit,QFeatures-method
#'
#' @author Lucas Beerland
#'
#' @param sce Single cell experiment object
#'
#' @param feature The name of feature (or the number of the row of the feature)
#' that should be visualized
#' @param plot hist or diff, specify the plot you want
#'
#' @param shade if TRUE differential regions are marked with shades
#'
#' @return ggplot object or list containing ggplot objects for pairwise comparisons
#'
#' @export
#'
#' @rdname visFit

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
#'
#' @export
#'
#' @rdname visFit

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
