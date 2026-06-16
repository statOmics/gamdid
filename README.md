# gamdid

<!-- badges: start -->
<!-- badges: end -->

`gamdid` is an R/Bioconductor package for differential distribution (DD)
analysis in single-cell proteomics data. It detects features whose entire
intensity distribution differs between experimental conditions — not just the
mean — using Poisson generalized additive models (GAMs) combined with
msqrob2 location tests via the harmonic mean p-value.

## Installation

Install the release version from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("gamdid")
```

Install the development version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("lbeerland/gamdid", dependencies = TRUE)
```

Setting `dependencies = TRUE` ensures all packages listed in `Suggests:`,
including those needed to reproduce the vignette, are installed.

## Usage

```r
library(gamdid)

## Load example data
data(pe)

## Compare distributions between two groups, correcting for batch
pe_dd <- compDistr(pe, covariate = "SampleType", nuisance = ~ Set)

## Visualize the most significant feature
rd <- as.data.frame(SummarizedExperiment::rowData(pe_dd))
top <- rownames(rd)[which.min(rd$SampleTypeMelanoma_VS_SampleTypeMonocyte)]
visFit(pe_dd, feature = top)
```

## Getting help

Please report bugs and feature requests at
<https://github.com/lbeerland/gamdid/issues>.
