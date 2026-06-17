#' Retrieve top-ranked features after differential distribution analysis
#'
#' @description
#' \code{topFeatures} extracts the \code{n} most significant features from
#' a \code{SummarizedExperiment} or \code{QFeatures} object returned by
#' \code{\link{compDistr}}, returning them as a \code{data.frame} sorted by
#' an (adjusted) p-value column from \code{rowData}.
#'
#' @details
#' The column used for sorting is determined as follows:
#'
#' \enumerate{
#'   \item If \code{sort.by} is supplied, that column is used directly.
#'   \item If \code{sort.by = NULL} (default) and the output was produced with
#'     \code{adjust = "stagewise"} or \code{adjust = "contrastwise"}:
#'     \itemize{
#'       \item For a \strong{two-group} comparison, the single pairwise
#'         comparison column (e.g. \code{"SampleTypeMelanoma_VS_SampleTypeMonocyte"})
#'         is used automatically.
#'       \item For a \strong{multi-group} comparison, the \code{"omnibus"}
#'         column is used automatically.
#'     }
#'   \item If \code{sort.by = NULL} and the output was produced with
#'     \code{adjust = "none"}, an error is thrown because there is no
#'     unambiguous default column. The user must specify \code{sort.by}
#'     explicitly, choosing from the \code{_GAM}, \code{_msqrob2}, or
#'     \code{_aggregated} columns (e.g. \code{sort.by = "omnibus_aggregated"}).
#' }
#'
#' Features for which the chosen p-value column is \code{NA} (typically
#' because fewer than \code{minObs} cells were observed in at least one group)
#' are placed at the end of the table.
#'
#' @include topFeaturesSCE.R
#'
#' @aliases topFeatures topFeatures,SummarizedExperiment-method topFeatures,QFeatures-method
#'
#' @author Lucas Beerland
#'
#' @param object A \code{SummarizedExperiment} (including
#'   \code{SingleCellExperiment}) or \code{QFeatures} instance returned by
#'   \code{\link{compDistr}}.
#'
#' @param sort.by A \code{character(1)} naming the \code{rowData} column to
#'   sort by. If \code{NULL} (default), the column is chosen automatically
#'   from the analysis metadata — see Details.
#'
#' @param n A \code{integer(1)} specifying the maximum number of features to
#'   return. Default is \code{10}.
#'
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{data.frame} with at most \code{n} rows, one per feature,
#'   sorted in ascending order of the chosen p-value column (most significant
#'   first). All columns from \code{rowData} are included.
#'
#' @seealso \code{\link{compDistr}} for producing the input object.
#'   \code{\link{visFit}} for visualizing individual top features.
#'
#' @examples
#' # Load example data and run distributional comparison
#' data(pe)
#' pe_dd <- compDistr(pe, covariate = "SampleType", nuisance = ~ Set)
#'
#' # Default: automatically sorts by the single pairwise comparison column
#' topFeatures(pe_dd)
#'
#' # Return the top 5 features
#' topFeatures(pe_dd, n = 5)
#'
#' # Specify the sort column explicitly
#' topFeatures(pe_dd, sort.by = "SampleTypeMelanoma_VS_SampleTypeMonocyte")
#'
#' @export
#'
#' @rdname topFeatures
#'
#' @importFrom MultiAssayExperiment getWithColData metadata
#' @importFrom methods is

setMethod(
    "topFeatures", "SummarizedExperiment",
    function(object, sort.by = NULL, n = 10) {
        topFeaturesSCE(sce = object, sort.by = sort.by, n = n)
    }
)

#' @param i A \code{character(1)} or \code{integer(1)} specifying which assay
#'   within the \code{QFeatures} object to retrieve results from.
#'
#' @export
#' @rdname topFeatures

setMethod(
    "topFeatures", "QFeatures",
    function(object, i, sort.by = NULL, n = 10) {
        sce <- getWithColData(object, i)
        topFeaturesSCE(sce = sce, sort.by = sort.by, n = n)
    }
)
