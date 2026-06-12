# `gamdid`

<!-- badges: start -->
<!-- badges: end -->

`gamddid` is a package for differential distribution analysis in single cell experiments

## Installation

`gamdid` can be easily installed from GitHub, you can use:

```r
remotes::install_github("statomics/gamdid", dependencies = TRUE)
```

Setting `dependencies = TRUE` should ensure that all packages, including the 
ones in the `Suggests:` field of the `DESCRIPTION`, are installed - 
this can be essential if you want to reproduce the code in the vignette, for example.
