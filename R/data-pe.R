#' Example single-cell proteomics dataset
#'
#' @description
#' A \code{SummarizedExperiment} object containing a subset of protein-level
#' data from Leduc et al. (2022), generated with prioritized Single-Cell
#' Proteomics by Mass Spectrometry (pSCoPE). The dataset is included to
#' illustrate the \code{\link{compDistr}} and \code{\link{visFit}} workflow
#' in the \pkg{gamdid} package.
#'
#' @format A \code{SummarizedExperiment} object with 100 proteins (rows) and
#'   1531 single cells (columns). Key slots:
#'   \describe{
#'     \item{\code{assay}}{Assay \code{"intensity"}: log2-transformed
#'       protein intensities, centered by the median of each cell
#'       (proteins \eqn{\times} cells), with missing values.}
#'     \item{\code{rowData}}{Feature-level metadata, including protein
#'       identifiers.}
#'     \item{\code{colData}}{Sample-level metadata, including the following
#'       columns:
#'       \describe{
#'         \item{\code{SampleType}}{Cell type label: \code{"Melanoma"} or
#'           \code{"Monocyte"}.}
#'         \item{\code{Set}}{MS run identifier (130 runs), used as a batch
#'           covariate in examples.}
#'         \item{\code{Field}}{\code{Field} annotation of the original data
#'           (1 to 4), used to illustrate the comparison of more than two
#'           groups.}
#'       }
#'     }
#'   }
#'
#' @details
#' The data were obtained from the \pkg{scpdata} Bioconductor package via
#' \code{scpdata::leduc2022_pSCoPE()}. Starting from the PSM data, low-quality
#' PSMs and cells were removed, the PSMs were summarized to proteins within
#' each MS run by median sweeping, the runs were joined and each cell was
#' centered by its median. \code{pe} contains the first 100 proteins for all
#' 1531 cells that passed quality control, for fast, reproducible examples.
#' The script that creates \code{pe} is available in the \code{scripts}
#' directory of the installed package:
#' \code{system.file("scripts", "make-pe.R", package = "gamdid")}.
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
