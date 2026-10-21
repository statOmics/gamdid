#' @importFrom SummarizedExperiment assays assays<- assay rowData rowData<-

#' @include locationTest.R partialResiduals.R shapeTest.R
#' @include aggregatePvals.R adjustOutput.R

compDistrSCE <- function(sce, covariate, nuisance = NULL, minObs = 50,
    bins = 25, adjust = "stagewise", alpha = 0.05, ...) {
    checkCompDistrArgs(sce, covariate, nuisance, minObs, bins, adjust)

    sce_original <- sce

    ## Location analysis
    location <- locationTest(
        sce = sce,
        covariate = covariate,
        nuisance = nuisance,
        ...
    )

    ## Shape analysis (on partial residuals when nuisance is specified)
    if (!is.null(nuisance)) {
        sce <- partialResiduals(
            sce = sce,
            covariate = covariate,
            msqrob_models = location$msqrob_models,
            msqrob_formula = location$msqrob_formula
        )
    }

    shape <- shapeTest(
        sce = sce,
        covariate = covariate,
        bins = bins,
        minObs = minObs
    )

    ## Aggregate shape and location p-values
    pval <- aggregatePvals(
        pval_shape = shape$pvals,
        pval_location = location$pvals
    )

    ## Apply multiple testing correction
    out <- adjustOutput(
        pval = pval,
        pval_shape = shape$pvals,
        pval_location = location$pvals,
        adjust = adjust,
        alpha = alpha
    )

    assembleOutput(
        sce_original, sce, nuisance, out, location, covariate, bins, alpha
    )
}

## Output object: the input with p values in rowData, the partial residuals
## (when nuisance is specified) and the settings used by visFit in metadata
assembleOutput <- function(sce_original, sce, nuisance, out, location,
    covariate, bins, alpha) {
    sce_out <- sce_original
    if (!is.null(nuisance)) assays(sce_out)$partialResiduals <- assay(sce)
    rowData(sce_out) <- cbind(rowData(sce_out), out$pvals)
    metadata(sce_out)$covariate <- covariate
    metadata(sce_out)$comparisons <- setdiff(
        colnames(location$pvals), "omnibus"
    )
    metadata(sce_out)$bins <- bins
    metadata(sce_out)$adj_alpha <- out$adj_alpha
    metadata(sce_out)$alpha <- alpha
    metadata(sce_out)$n_tests <- out$n_tests
    sce_out
}

## Argument checks of compDistr
checkCompDistrArgs <- function(sce, covariate, nuisance, minObs, bins, adjust) {
    if (!is.character(adjust) || length(adjust) != 1 ||
        !adjust %in% c("stagewise", "contrastwise", "none")) {
        stop(
            "'adjust' must be one of \"stagewise\", \"contrastwise\" ",
            "or \"none\""
        )
    }
    if (!is.character(covariate) || length(covariate) != 1 ||
        !covariate %in% colnames(colData(sce))) {
        stop("'covariate' must be the name of one column of colData")
    }
    ## Each group is modelled with a smoothing spline of 10 basis functions,
    ## which requires at least 10 bins
    if (!is.numeric(bins) || length(bins) != 1 || is.na(bins) ||
        bins < 10 || bins %% 1 != 0) {
        stop("'bins' must be a single whole number of at least 10")
    }
    if (!is.numeric(minObs) || length(minObs) != 1 || is.na(minObs) ||
        minObs < 0) {
        stop("'minObs' must be a single non-negative number")
    }
    if (covariate %in% all.vars(nuisance)) {
        stop("covariate of interest cannot be in nuisance formula")
    }
    ## Only random intercepts such as (1 | Set) are supported
    intercept_only <- vapply(
        findbars(nuisance),
        function(b) identical(b[[2]], 1), logical(1)
    )
    if (!all(intercept_only) ||
        any(grepl("||", deparse(nuisance), fixed = TRUE))) {
        stop(
            "only random intercepts such as (1 | Set) are supported in ",
            "'nuisance', not random slopes or '||'"
        )
    }
    if (!is.character(sce[[covariate]]) && !is.factor(sce[[covariate]])) {
        stop("covariate of interest must be a character or factor")
    }
}
