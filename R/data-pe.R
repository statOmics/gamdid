#' Example proteomics dataset
#'
#' A proteomics dataset subsampled from Leduc et al. (2022) generated using
#' prioritized Single-Cell Proteomics by Mass Spectrometry (pSCoPE).
#' The data contains quantitative measurements of proteins across
#' single cells and carrier channels acquired with a tandem mass tag (TMT)
#' labeling strategy.
#'
#' The dataset is provided as a \code{SingleCellExperiment} object and is commonly used
#' to illustrate preprocessing, normalization, and differential analysis
#' workflows in single-cell proteomics.
#'
#' @format A \code{SingleCellExperiment} object with:
#' \describe{
#'   \item{assays}{Quantitative protein intensity data}
#'   \item{rowData}{Feature-level metadata such as protein identifiers}
#'   \item{colData}{Sample-level metadata including cell type, carrier status, and experimental conditions}
#' }
#'
#' @details
#' The dataset includes both single-cell and carrier channels.
#'
#' @source
#' Leduc A, et al. (2022). \emph{Prioritized single-cell proteomics reveals
#' molecular and functional polarity in macrophages}. Molecular Systems Biology.
#'
#' @references
#' Leduc A, Huffman RG, Cantlon JD, et al. (2022).
#' Prioritized single-cell proteomics reveals molecular and functional
#' polarity in macrophages. Mol Syst Biol.
"pe"
