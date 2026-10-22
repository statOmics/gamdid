# gamdid

<!-- badges: start -->
<!-- badges: end -->

`gamdid` (generalized additive models for differential distributions) is an
R/Bioconductor package for differential distribution (DD) analysis in
single-cell proteomics data. It detects features whose entire intensity
distribution differs between experimental conditions — not just the mean — and
shows *where* in the distribution the groups differ.

Conventional differential abundance (DA) methods only test for a shift in the
mean. Single-cell data, however, also reveal differences in variance,
modality, skewness, or in the relative size of cell subpopulations. `gamdid`
is designed to detect all of these, while keeping the power of established DA
methods for pure mean shifts.

## How it works

For each feature (protein or peptide), `gamdid`:

1. **Estimates a smooth distribution per group.** The intensities are binned
   into a histogram with the same bins for every group, and a Poisson
   generalized additive model (GAM, fitted with
   [mgcv](https://cran.r-project.org/package=mgcv)) turns the bin counts into a
   smooth density estimate per group.
2. **Tests whether the distributions differ.** A Wald-type test on the
   group-specific smoothing splines gives a GAM p-value that is sensitive to
   any distributional difference.
3. **Combines it with a mean-shift test.** The GAM p-value is aggregated with
   a [msqrob2](https://bioconductor.org/packages/msqrob2) p-value using the
   harmonic mean p-value, which retains power for pure mean shifts.
4. **Corrects for multiple testing** across features (Benjamini–Hochberg) and,
   for more than two groups, across pairwise comparisons via stagewise testing
   ([stageR](https://bioconductor.org/packages/stageR)).
5. **Localizes the differences.** A difference smoother with simultaneous
   confidence bands highlights the regions of the distribution where the
   groups differ.

## Installation

`gamdid` will be submitted to Bioconductor. Until it is available there,
install the development version from GitHub:

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("statOmics/gamdid", dependencies = TRUE)
```

Setting `dependencies = TRUE` ensures all packages listed in `Suggests:`,
including those needed to build the vignette, are installed. `gamdid` requires
R >= 4.6.0.

Once the package is available on Bioconductor, it can be installed with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("gamdid")
```

## Input data

`gamdid` works on a `SummarizedExperiment` (including `SingleCellExperiment`)
or on an assay of a `QFeatures` object, where:

- rows are features (proteins or peptides) and columns are single cells;
- the assay contains **log-transformed, normalized intensities**; missing
  values (`NA`) are allowed;
- the `colData` contains the grouping variable to compare (a `character` or
  `factor` column), and optionally nuisance variables such as MS run or batch.

A feature is only tested when it has at least `minObs` (default 50)
non-missing observations in every group, since a distribution cannot be
estimated reliably from fewer cells.

## Main functions

| Function | Purpose |
|----------|---------|
| `compDistr()` | Runs the differential distribution analysis and adds (adjusted) p-values to the `rowData`. |
| `topFeatures()` | Returns the most significant features as a `data.frame`, sorted by (adjusted) p-value. |
| `visFit()` | Plots the histograms and fitted GAM densities of one feature, with the differential regions shaded in grey. |

Each function has a help page with details and examples, e.g. `?compDistr`.

## Quick start

```r
library(gamdid)

## Example data: protein-level intensities of melanoma cells and monocytes
## acquired across several MS runs
data(pe)

## Compare the distributions of the two cell types, accounting for MS run
pe_dd <- compDistr(pe, covariate = "SampleType", nuisance = ~ Set)

## The 10 most significant features
top <- topFeatures(pe_dd, n = 10)
top

## Visualize the most significant feature, with the difference smoother
## (visFit returns one plot for two groups, or a named list with one plot
## per comparison for more groups)
visFit(pe_dd, feature = rownames(top)[1], plot = "diff")
```

`compDistr()` returns the input object with the results added:

- **`rowData`**: one column of adjusted p-values per comparison (for two
  groups, a column named `<level1>_VS_<level2>`; for more than two groups,
  an `omnibus` column plus one column per pairwise comparison);
- **`assays`**: when `nuisance` is given, a `partialResiduals` assay with the
  nuisance-corrected intensities;
- **`metadata`**: the settings and the adjusted significance threshold used
  by `visFit()`.

Further options include:

- **More than two groups:** an omnibus test across all groups, followed by
  post-hoc pairwise comparisons (`adjust = "stagewise"`, the default, or
  `adjust = "contrastwise"`).
- **Nuisance variables:** fixed effects and random intercepts in `lme4`
  syntax, e.g. `nuisance = ~ (1 | run)`.
- **Separate p-values:** `adjust = "none"` returns the unadjusted GAM,
  msqrob2 and aggregated p-values for each comparison.

## Documentation

The vignette provides a complete walkthrough, covering two-group and
multi-group analyses, the adjustment strategies, ranking features and
visualization. After installing the package with `dependencies = TRUE`, open
it with:

```r
browseVignettes("gamdid")
```

## Citation

If you use `gamdid` in your research, please cite:

> Beerland L, Vandenbulcke S, Vanderaa C, Martens L, Clement L (2026).
> gamdid: generalized additive models for differential distributions in single
> cell experiments. *bioRxiv*.
> [doi:10.64898/2026.06.18.733106](https://doi.org/10.64898/2026.06.18.733106)

The citation is also available from R with `citation("gamdid")`, and in the
[`CITATION.cff`](CITATION.cff) file of this repository.

## License

`gamdid` is released under the
[Artistic License 2.0](LICENSE.md).

## Getting help

Please report bugs and feature requests at
<https://github.com/statOmics/gamdid/issues>.
