#' Example single-cell proteomics dataset
#'
#' @description
#' A \code{SingleCellExperiment} object containing a subset of protein-level
#' data from Leduc et al. (2022), generated with prioritized Single-Cell
#' Proteomics by Mass Spectrometry (pSCoPE). The dataset is included to
#' illustrate the \code{\link{compDistr}} and \code{\link{visFit}} workflow
#' in the \pkg{gamdid} package.
#'
#' @format A \code{SingleCellExperiment} object with rows corresponding to
#'   proteins and columns to cells/channels. Key slots:
#'   \describe{
#'     \item{\code{assay}}{Log-transformed, normalized protein intensity
#'       matrix (proteins \eqn{\times} cells).}
#'     \item{\code{rowData}}{Feature-level metadata, including protein
#'       identifiers.}
#'     \item{\code{colData}}{Sample-level metadata with the following columns:
#'       \describe{
#'         \item{\code{SampleType}}{Cell type label: \code{"Melanoma"} or
#'           \code{"Monocyte"}.}
#'         \item{\code{Set}}{MS run identifier, used as a batch covariate in
#'           examples.}
#'       }
#'     }
#'   }
#'
#' @details
#' The full Leduc 2022 pSCoPE dataset is available from the
#' \pkg{scpdata} Bioconductor package via \code{scpdata::leduc2022_pSCoPE()}.
#' The \code{pe} object is a small subset (200 proteins, a selection of cells)
#' retained at the \code{proteins_processed} level for fast, reproducible
#' examples and unit tests.
#'
#' @usage data(pe)
#'
#' @source
#' Leduc A, et al. (2022). Prioritized single-cell proteomics reveals
#' molecular and functional polarity in macrophages.
#' \emph{Molecular Systems Biology}.
#' \doi{10.15252/msb.202110798}
#'
#' @references
#' Leduc A, Huffman RG, Cantlon JD, Khan S, Slavov N (2022).
#' Prioritized single-cell proteomics reveals molecular and functional
#' polarity in macrophages. \emph{Mol Syst Biol} \strong{18}, e10798.
#'
#' @seealso \code{\link{compDistr}}, \code{\link{visFit}}
#'
#' @keywords datasets
"pe"
