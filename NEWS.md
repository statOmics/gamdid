# gamdid 0.99.0

## NEW FEATURES

* Initial submission to Bioconductor.
* `compDistr()` performs differential distribution analysis on
  `SummarizedExperiment` and `QFeatures` objects using Poisson GAMs combined
  with msqrob2 location tests via the harmonic mean p-value.
* `visFit()` produces interpretable histogram overlays and difference smoother
  plots, highlighting differential regions of the distribution with grey
  shading.
* Support for nuisance correction via fixed effects and random effects (lme4
  syntax) in both the location and shape analyses.
* Omnibus testing across more than two groups with post-hoc pairwise
  comparisons corrected by stagewise testing (`stageR`) or Benjamini-Hochberg
  FDR.
