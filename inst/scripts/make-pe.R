## Creation of the `pe` example dataset of gamdid
##
## `pe` contains the protein intensities of the first 100 proteins of the
## pSCoPE data of Leduc et al. (2022) for 1531 single cells (melanoma cells
## and monocytes), as distributed by the scpdata package
## (scpdata::leduc2022_pSCoPE(), ExperimentHub record EH8301).
##
## The PSM data are minimally processed: quality control of PSMs and cells,
## then summarization to proteins within each MS run by median sweeping
## (log2-transform, subtract the median of each PSM, take the median per
## sample), joining all runs and centering each cell by its median.
##
## Packages: scpdata, scp, QFeatures, MsCoreUtils, matrixStats,
## SummarizedExperiment (Bioconductor 3.22/3.23).

library(scpdata)
library(scp)
library(QFeatures)
library(SummarizedExperiment)
library(matrixStats)

## Data -----------------------------------------------------------------------

leduc <- leduc2022_pSCoPE()

## Remove the assays processed by the authors: only the PSM data (one assay
## per MS run) are kept
leduc <- removeAssay(leduc, c("peptides", "peptides_log",
                              "proteins_norm2", "proteins_processed"))
psmAssays <- names(leduc)

## Zeros are missing values
leduc <- zeroIsNA(leduc, i = names(leduc))

## PSM quality control ----------------------------------------------------------

leduc <- computeSCR(leduc, names(leduc),
                    colvar = "SampleType",
                    samplePattern = "Mel|Macro",
                    carrierPattern = "Carrier",
                    sampleFUN = "mean",
                    rowDataName = "MeanSCR")
leduc <- filterFeatures(leduc, ~ Reverse != "+" &
                            Potential.contaminant != "+" &
                            !grepl("REV|CON", Leading.razor.protein) &
                            !is.na(PIF) & PIF > 0.6 &
                            dart_qval < 0.01 &
                            !is.na(MeanSCR) & MeanSCR < 0.05)

## Cell quality control ---------------------------------------------------------

## Keep single cells with a median CV below 0.6
leduc <- medianCVperCell(leduc,
                         i = names(leduc),
                         groupBy = "Leading.razor.protein.symbol",
                         nobs = 3,
                         na.rm = TRUE,
                         colDataName = "MedianCV",
                         norm = "sum")
leduc <- subsetByColData(leduc, !is.na(leduc$MedianCV) &
                             leduc$MedianCV < 0.6 &
                             grepl("Mono|Mel", leduc$SampleType))

## Remove cells with more than 95% missing peptides
peptideAssays <- paste0("peptides_", psmAssays)
leduc <- aggregateFeaturesOverAssays(leduc,
                                     i = psmAssays,
                                     fcol = "Sequence",
                                     name = peptideAssays,
                                     fun = colMedians,
                                     na.rm = TRUE)
leduc <- joinAssays(leduc, i = peptideAssays, name = "peptides")
leduc <- logTransform(leduc, i = "peptides", name = "peptides_log")
pmissc <- colMeans(is.na(assay(leduc[["peptides_log"]])))
leduc <- subsetByColData(leduc, pmissc < 0.95)

## Protein summarization by median sweeping -------------------------------------

sumFunctionMedianSweep <- function(mx) {
    mx <- log2(mx)
    mx <- sweep(mx, MARGIN = 1, FUN = "-",
                STATS = rowMedians(mx, na.rm = TRUE))
    colMedians(mx)
}
proteinAssays <- paste0("proteinSweepNonNorm_", psmAssays)
leduc <- aggregateFeaturesOverAssays(leduc,
                                     i = psmAssays,
                                     fcol = "Leading.razor.protein.id",
                                     name = proteinAssays,
                                     fun = sumFunctionMedianSweep)
leduc <- joinAssays(leduc, i = proteinAssays, name = "proteinSweepNonNorm")

## Center each cell by its median
leduc <- sweep(leduc,
               i = "proteinSweepNonNorm",
               name = "proteinSweep",
               MARGIN = 2,
               FUN = "-",
               STATS = colMedians(assay(leduc[["proteinSweepNonNorm"]]),
                                  na.rm = TRUE))
leduc <- filterNA(leduc, i = "proteinSweep", pNA = 0.95)

## Example dataset --------------------------------------------------------------

## First 100 proteins and all cells, stored as a SummarizedExperiment
## (the constructor is used because as(x, "SummarizedExperiment") drops the
## row names of a SingleCellExperiment)
sce <- getWithColData(leduc, "proteinSweep")[1:100, ]
pe <- SummarizedExperiment(assays = list(intensity = assay(sce)),
                           rowData = rowData(sce),
                           colData = colData(sce))

save(pe, file = "data/pe.rda", compress = "xz")
