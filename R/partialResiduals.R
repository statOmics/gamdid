#' @importFrom SummarizedExperiment colData assay assay<-
#' @importFrom Matrix cbind2
#' @include mixedDesign.R
#' @importFrom stats model.matrix
#' @importFrom msqrob2 getCoef
#' @importFrom reformulas findbars

partialResiduals <- function(sce, covariate, msqrob_models, msqrob_formula) {
    ## Design matrix columns of the covariate of interest (exact names)
    covariate_cols <- paste0(covariate, unique(colData(sce)[[covariate]]))

    if (!is.null(findbars(msqrob_formula))) {
        ## Mixed model: X and Z design matrices
        design <- mixedDesign(msqrob_formula, colData(sce))
        X <- design$X
        Z <- design$Z

        ## Zero out covariate columns; keep intercept to preserve overall mean
        X[, colnames(X) %in% covariate_cols] <- 0
        X[, 1] <- 1

        XZ <- cbind2(X, Z)

        full_coef <- c(colnames(X), design$colnames_Z)

        coefs <- nuisanceCoefs(msqrob_models, full_coef)

        assay(sce) <- assay(sce) - t(as.matrix(XZ %*% coefs))
    } else {
        ## Fixed-effects model
        X <- model.matrix(msqrob_formula, colData(sce))
        X[, colnames(X) %in% covariate_cols] <- 0
        X[, 1] <- 1

        full_coef <- colnames(X)

        coefs <- nuisanceCoefs(msqrob_models, full_coef)

        assay(sce) <- assay(sce) - t(X %*% coefs)
    }

    sce
}

## Coefficients of all models in the order of full_coef (one column per model)
nuisanceCoefs <- function(msqrob_models, full_coef) {
    vapply(msqrob_models, function(m) {
        cfs <- getCoef(m)
        ## Pad missing coefficients with zero for matrix multiplication
        if (all(is.na(cfs))) {
            return(rep(NA_real_, length(full_coef)))
        }
        cfs[setdiff(full_coef, names(cfs))] <- 0
        ## Non-estimable coefficients (e.g. a nuisance level without
        ## observations for this feature) are NA: treat them as zero
        cfs[is.na(cfs)] <- 0
        cfs[full_coef]
    }, FUN.VALUE = numeric(length(full_coef)))
}
