#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats poisson setNames
#'
#' @include shapeContrast.R waldTest.R

shapeTest <- function(sce, covariate, minObs = 50, bins = 25, tol = 1e-8) {
    colDat <- colData(sce)
    X <- assay(sce)

    cov_values <- colDat[[covariate]]
    cov_levels <- unique(cov_values)
    n_levels <- length(cov_levels)

    model_list <- fitShapeModels(X, colDat, covariate, minObs, bins)
    contrast <- shapeContrastFromFits(model_list, cov_levels)
    pvals <- shapeWaldTests(
        model_list, contrast$ctr, contrast$nbasis, n_levels, tol
    )
    colnames(pvals) <- renameShapeContrasts(colnames(pvals), covariate)

    list(models = model_list, pvals = pvals)
}

## Poisson GAM per feature, fitted to its histogram counts (NA if a feature
## has too few observations)
fitShapeModels <- function(X, colDat, covariate, minObs, bins) {
    ## Build histogram data per feature
    hist_dat <- apply(X,
        MARGIN = 1, FUN = createCounts,
        colDat = colDat, covariate = covariate,
        minObs = minObs, bins = bins
    )

    ## Fit a Poisson GAM per feature
    lapply(hist_dat, function(dat) {
        if (all(is.na(dat)) || any(is.infinite(dat$carrier))) {
            return(NA)
        }
        mgcv::gam(
            counts ~ covariate + s(midsC, by = covariate) + offset(carrier),
            data = dat, family = poisson(), method = "REML"
        )
    })
}

## Contrast matrix, derived from the first successfully fitted model
shapeContrastFromFits <- function(model_list, cov_levels) {
    n_levels <- length(cov_levels)
    is_gam <- vapply(model_list, function(x) inherits(x, "gam"), logical(1))
    idx <- which(is_gam)[1]
    ref_fit <- model_list[[idx]]

    if (is.null(ref_fit)) {
        stop("No valid GAM fits available to perform inference.")
    }

    p <- length(ref_fit$coefficients)
    nbasis <- (p - n_levels) / n_levels
    coef_names <- names(ref_fit$coefficients[seq_len(n_levels)])

    ## Replace "(Intercept)" with the reference level name
    coef_names_full <- paste0("covariate", cov_levels)
    if (!all(coef_names_full %in% coef_names)) {
        coef_names[coef_names == "(Intercept)"] <-
            coef_names_full[which(!(coef_names_full %in% coef_names))]
    }
    list(ctr = shapeContrast(n_levels, nbasis, coef_names), nbasis = nbasis)
}

## Wald tests per feature: one comparison, or omnibus + pairwise comparisons
shapeWaldTests <- function(model_list, ctr, nbasis, n_levels, tol) {
    if (ncol(ctr$comps) < 2) {
        ## Single pairwise comparison
        pvals <- vapply(model_list, function(fit) {
            if (all(is.na(fit))) {
                return(NA_real_)
            }
            edf <- summary(fit)$edf
            df_omni <- sum(sort(edf, decreasing = TRUE)[seq_len(n_levels - 1)])
            res <- waldTest(fit, ctr$L, df = df_omni, tol = tol)
            names(res) <- colnames(ctr$comps)
            res
        }, FUN.VALUE = numeric(1))
        pvals <- matrix(pvals,
            nrow = 1,
            dimnames = list(colnames(ctr$comps), NULL)
        )
    } else {
        ## Omnibus test + pairwise comparisons
        pvals <- vapply(model_list, function(fit) {
            if (all(is.na(fit))) {
                return(setNames(
                    rep(NA_real_, ncol(ctr$comps) + 1),
                    c("omnibus", colnames(ctr$comps))
                ))
            }
            edf <- summary(fit)$edf
            df_omni <- sum(sort(edf, decreasing = TRUE)[seq_len(n_levels - 1)])
            res <- waldTest(fit, ctr$L, df = df_omni, tol = tol)
            names(res)[names(res) == "pval"] <- "omnibus"
            c(res, pairwiseShapeTests(
                fit, ctr$L, ctr$comps, nbasis,
                tol = tol
            ))
        }, FUN.VALUE = numeric(ncol(ctr$comps) + 1))
    }

    if (is.matrix(pvals)) {
        pvals <- t(pvals)
    }
    pvals
}

## Rename the GAM factor "covariate" to the covariate of interest, only at
## the start of each level name of a contrast
renameShapeContrasts <- function(contrasts, covariate) {
    vapply(
        strsplit(contrasts, "_VS_", fixed = TRUE),
        function(parts) {
            parts <- ifelse(startsWith(parts, "covariate"),
                paste0(covariate, substring(parts, nchar("covariate") + 1)),
                parts
            )
            paste(parts, collapse = "_VS_")
        }, character(1)
    )
}
