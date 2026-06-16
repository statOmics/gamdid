#' Interpretable visualization of the GAM fit for a single feature
#'
#' @description
#' \code{visFit} produces a publication-ready plot for a single feature after
#' running \code{\link{compDistr}}. It overlays the empirical histograms and
#' the fitted Poisson GAM density curves for each group, and optionally shades
#' the regions of the distribution that differ significantly between groups.
#'
#' @details
#' Two plot types are available:
#'
#' \describe{
#'   \item{\code{"hist"}}{Displays the normalized histogram bars (empirical
#'     density) and the fitted GAM density curves side-by-side for the two
#'     groups. Grey shaded bands highlight bins where the pointwise confidence
#'     interval for the log-ratio of the two densities excludes zero, i.e.\
#'     regions where the groups differ significantly in density.}
#'
#'   \item{\code{"diff"}}{A two-panel plot. The top panel shows the same
#'     histogram and GAM density curves as \code{"hist"}. The bottom panel
#'     displays the log-ratio of the two fitted GAM densities (the
#'     \emph{difference smoother}) together with pointwise confidence bands.
#'     The dashed red horizontal line marks zero (no difference). Shaded
#'     regions indicate where the confidence band excludes zero.}
#' }
#'
#' Confidence bands are Bonferroni-corrected pointwise bands based on the
#' effective degrees of freedom of the GAM fit (Nychka 1988), ensuring that
#' the shaded regions have the correct simultaneous coverage. When
#' \code{shade = TRUE} and \code{alpha = NULL}, the significance threshold
#' is taken from \code{metadata(object)$adj_alpha}, the BH-adjusted threshold
#' stored by \code{\link{compDistr}}.
#'
#' When \code{nuisance} was specified in \code{\link{compDistr}}, the
#' partial-residual assay (nuisance-corrected intensities) is used for
#' visualization, so that group differences are shown on the nuisance-adjusted
#' scale.
#'
#' For more than two groups, one plot per pairwise comparison is produced and
#' returned as a named list.
#'
#' @include visFitSCE.R
#'
#' @aliases visFit visFit,SummarizedExperiment-method visFit,QFeatures-method
#'
#' @author Lucas Beerland
#'
#' @param object A \code{SummarizedExperiment} (including
#'   \code{SingleCellExperiment}) or \code{QFeatures} instance that has been
#'   processed by \code{\link{compDistr}}.
#'
#' @param feature A \code{character(1)} giving the row name of the feature to
#'   visualize, or an \code{integer(1)} giving its row index.
#'
#' @param plot A \code{character(1)} specifying the plot type: \code{"hist"}
#'   (default) for the histogram and GAM curves only, or \code{"diff"} for a
#'   two-panel plot that additionally shows the difference smoother and its
#'   confidence band.
#'
#' @param shade A \code{logical(1)}. If \code{TRUE} (default), regions where
#'   the two groups differ significantly are highlighted with grey shading.
#'
#' @param alpha A \code{numeric(1)} in \code{(0, 1)} specifying the
#'   significance level used to construct the confidence band for the
#'   difference smoother. If \code{NULL} (default), the adjusted threshold
#'   \code{metadata(object)$adj_alpha} set by \code{\link{compDistr}} is used.
#'
#' @return Invisibly returns a named list of \code{ggplot} objects, one per
#'   pairwise comparison. For a two-group comparison the list has a single
#'   element. The plots are also printed as a side effect.
#'
#' @seealso \code{\link{compDistr}} which must be run before \code{visFit}.
#'
#' @examples
#' # Load example data: a SummarizedExperiment with protein-level intensities
#' # for two cell types (Melanoma and Monocyte) across multiple MS runs.
#' data(pe)
#'
#' # Run distributional comparison
#' pe_dd <- compDistr(pe, covariate = "SampleType", nuisance = ~ Set)
#'
#' # Histogram overlay with differential region shading (default)
#' visFit(pe_dd, feature = 20)
#'
#' # Two-panel plot with difference smoother
#' visFit(pe_dd, feature = 20, plot = "diff")
#'
#' # Use a custom significance threshold for the confidence band
#' visFit(pe_dd, feature = 20, alpha = 0.01)
#'
#' @export
#'
#' @rdname visFit
#'
#' @importFrom MultiAssayExperiment getWithColData metadata metadata<-
#' @importFrom methods is

setMethod(
    "visFit", "SummarizedExperiment",
    function(object, feature, plot = "hist", shade = TRUE, alpha = NULL) {
        visFitSCE(sce = object, feature = feature,
            plot = plot, shade = shade, alpha = alpha)
    }
)

#' @param i A \code{character(1)} or \code{integer(1)} specifying which assay
#'   within the \code{QFeatures} object to visualize.
#'
#' @export
#'
#' @rdname visFit

setMethod(
    "visFit", "QFeatures",
    function(object, i, feature, plot = "hist", shade = TRUE, alpha = NULL) {
        sce <- getWithColData(object, i)
        visFitSCE(sce = sce, feature = feature,
            plot = plot, shade = shade, alpha = alpha)
    }
)
