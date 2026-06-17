#' Compare distributions across groups in single-cell proteomics data
#'
#' @description
#' \code{compDistr} detects differentially distributed features between two or
#' more experimental groups. The sample space of each feature is discretized
#' into a histogram, and a Poisson generalized additive model (GAM) is fitted
#' per group using \pkg{mgcv}. A Wald-type chi-squared test on the
#' group-specific smoothing splines yields a GAM p-value for any type of
#' distributional difference. To increase power for mean shifts, the GAM
#' p-value is aggregated with a p-value from the \pkg{msqrob2} framework via
#' the harmonic mean p-value.
#'
#' @details
#' The analysis proceeds in three steps:
#'
#' \enumerate{
#'   \item \strong{GAM analysis.} When a \code{nuisance} formula is provided,
#'     partial residuals are computed by removing the nuisance effect from the
#'     assay values before discretizing. A Poisson GAM with group-specific
#'     smoothing splines is then fitted to the histogram counts. A Wald-type
#'     chi-squared test on the spline coefficients yields a p-value for any
#'     distributional difference (both in location and in shape).
#'
#'   \item \strong{msqrob2 analysis.} A linear model is fitted using
#'   \code{\link[msqrob2]{msqrob}}, yielding a p-value for a mean difference
#'   between groups (one per pairwise contrast, plus an omnibus p-value when
#'   there are more than two groups).
#'
#'   \item \strong{P-value aggregation.} The GAM and msqrob2 p-values are
#'     combined per contrast using the harmonic mean p-value
#'     (\code{\link[harmonicmeanp]{hmp.stat}}), which is valid under arbitrary
#'     dependence and provides power gains when both tests detect signal.
#' }
#'
#' Multiple testing correction is applied across features after aggregation.
#'
#' @include compDistrSCE.R
#'
#' @aliases compDistr compDistr,SummarizedExperiment-method compDistr,QFeatures-method
#'
#' @author Lucas Beerland
#'
#' @param object A \code{SummarizedExperiment} (including
#'   \code{SingleCellExperiment}) or \code{QFeatures} instance. The default
#'   assay should contain log-transformed protein or peptide intensities.
#'
#' @param i A \code{character(1)} or \code{integer(1)} specifying which assay
#'    within the \code{QFeatures} object to analyse.
#'
#' @param covariate A \code{character(1)} naming a column in the
#'   \code{colData} of \code{object}. The levels of this column define the
#'   groups to be compared. Must be a \code{character} or \code{factor}
#'   variable.
#'
#' @param nuisance An optional one-sided formula (e.g. \code{~ (1|run)})
#'   specifying nuisance covariates to account for in the analyses.
#'   Fixed effects and random effects (specified with
#'   \code{lme4} syntax) are supported. The \code{covariate} of interest must
#'   not appear in this formula.
#'
#' @param minObs A \code{numeric(1)} giving the minimum number of non-missing
#'   observations required \emph{per group} to attempt GAM fitting for a
#'   feature. Features with fewer observations in any group are assigned
#'   \code{NA} p-values. Default is \code{50}.
#'
#' @param bins A \code{numeric(1)} giving the number of histogram bins used to
#'   discretize the sample space for the GAM shape test. Finer grids capture
#'   more detailed shape differences but require more observations per bin.
#'   Default is \code{25}.
#'
#' @param adjust A \code{character(1)} specifying the multiple testing
#'   correction strategy:
#'   \describe{
#'     \item{\code{"stagewise"}}{(default) For two groups, applies
#'       Benjamini–Hochberg FDR correction. For more than two groups, applies
#'       a two-stage procedure via \pkg{stageR}: in the first stage, the
#'       omnibus p-value is adjusted across features with Benjamini–Hochberg
#'       FDR control. Only features that pass this screening stage proceed to
#'       the second stage, where the pairwise comparisons are adjusted for
#'       multiple testing across the comparisons via the Holm method, within
#'       each significant feature.}
#'     \item{\code{"contrastwise"}}{Applies Benjamini–Hochberg FDR correction
#'       independently to each contrast (omnibus and all pairwise comparisons).}
#'     \item{\code{"none"}}{Returns raw (unadjusted) p-values in three columns
#'       per contrast: \code{_GAM} (GAM test), \code{_msqrob2} (msqrob2 test),
#'       and \code{_aggregated} (harmonic mean of the two), without any
#'       multiple testing correction.}
#'   }
#'
#' @param alpha A \code{numeric(1)} specifying the overall significance level.
#'   Used by stagewise testing and to determine the adjusted significance
#'   threshold stored in \code{metadata(result)$adj_alpha}. Default is
#'   \code{0.05}.
#'
#' @param ... Additional arguments passed to \code{\link[msqrob2]{msqrob}},
#'   e.g. \code{robust = TRUE} to use robust regression for the location
#'   analysis.
#'
#' @return For the \code{QFeatures} method: the input \code{QFeatures} object
#'   with the assay specified by \code{i} replaced by an updated
#'   \code{SummarizedExperiment} containing the same additions as described
#'   below for the \code{SummarizedExperiment} method.
#'
#' @seealso
#' \code{\link{visFit}} for visualizing the GAM fits.
#' \code{\link[msqrob2]{msqrob}} for the underlying location model.
#' \code{\link[stageR]{stageWiseAdjustment}} for stagewise testing.
#'
#' @examples
#' # Load example data: a SummarizedExperiment with protein-level intensities
#' # for two cell types (Melanoma and Monocyte) across multiple MS runs.
#' data(pe)
#'
#' # Two-group comparison, correcting for MS run (batch) effect
#' pe_dd <- compDistr(pe, covariate = "SampleType", nuisance = ~ Set)
#'
#' # View the adjusted p-values for the pairwise comparison
#' SummarizedExperiment::rowData(pe_dd)$SampleTypeMelanoma_VS_SampleTypeMonocyte
#'
#' # Retrieve the BH-adjusted significance threshold
#' MultiAssayExperiment::metadata(pe_dd)$adj_alpha
#'
#' @export
#'
#' @rdname compDistr
#'
#' @importFrom MultiAssayExperiment getWithColData metadata metadata<-
#' @importFrom methods is

setMethod(
  "compDistr", "QFeatures",
  function(object, i, covariate, nuisance = NULL, minObs = 50,
           bins = 25, adjust = "stagewise", alpha = 0.05, ...) {
    sce <- getWithColData(object, i)
    sce_out <- compDistrSCE(sce = sce, covariate = covariate,
                            nuisance = nuisance, minObs = minObs, bins = bins,
                            adjust = adjust, alpha = alpha, ...)
    object[[i]] <- sce_out
    object
  }
)


#' @return For the \code{SummarizedExperiment} method: the input
#'   \code{SummarizedExperiment} augmented with:
#'   \describe{
#'     \item{\code{rowData}}{One or more columns of (adjusted) p-values. For
#'       two groups, a single column named
#'       \code{<level1>_VS_<level2>}. For more than two groups, a column
#'       \code{omnibus} plus one column per pairwise comparison. When
#'       \code{adjust = "none"}, three columns per contrast are added:
#'       \code{_GAM}, \code{_msqrob2}, and \code{_aggregated}.}
#'     \item{\code{assays}}{If \code{nuisance} is provided, a
#'       \code{partialResiduals} assay is added containing nuisance-corrected
#'       intensities used for visualization.}
#'     \item{\code{metadata}}{Entries \code{covariate}, \code{comparisons},
#'       \code{bins}, and \code{adj_alpha} are added, which are used by
#'       \code{\link{visFit}}.}
#'   }
#'
#' @export
#' @rdname compDistr

setMethod(
  "compDistr", "SummarizedExperiment",
  function(object, covariate, nuisance = NULL, minObs = 50,
           bins = 25, adjust = "stagewise", alpha = 0.05, ...) {
    compDistrSCE(sce = object, covariate = covariate,
                 nuisance = nuisance, minObs = minObs, bins = bins,
                 adjust = adjust, alpha = alpha, ...)
  }
)
