#' @importFrom SummarizedExperiment colData assay assay<-
#' @importFrom Matrix cbind2
#' @importFrom lme4 ranef getME lmer
#' @importFrom stats model.matrix
#' @importFrom msqrob2 getCoef
#' @importFrom reformulas findbars

partialResiduals <- function(sce, covariate, msqrob_models, msqrob_formula) {
    if (!is.null(findbars(msqrob_formula))) {
        ## Mixed model: use lme4 to recover X and Z design matrices
        mock <- data.frame(colData(sce), y = 0)
        mock_fit <- lmer(
            paste0("y~", as.character(msqrob_formula)[2]), data = mock
        )
        X <- getME(mock_fit, "X")
        Z <- getME(mock_fit, "Z")

        ## Zero out covariate columns; keep intercept to preserve overall mean
        X[, grepl(covariate, colnames(X))] <- 0
        X[, 1] <- 1

        XZ <- cbind2(X, Z)

        ## Construct Z column names matching msqrob output
        colnames_Z <- paste0("(Intercept)",
            unlist(lapply(names(ranef(mock_fit)), function(nm) {
                df <- ranef(mock_fit)[[nm]]
                paste0(nm, rownames(df))
            }))
        )
        full_coef <- c(colnames(X), colnames_Z)

        coefs <- vapply(msqrob_models, function(m) {
            cfs <- getCoef(m)
            ## Pad missing coefficients with zero for matrix multiplication
            if (all(is.na(cfs))) return(rep(NA_real_, length(full_coef)))
            cfs[setdiff(full_coef, names(cfs))] <- 0
            cfs[full_coef]
        }, FUN.VALUE = numeric(length(full_coef)))

        assay(sce) <- assay(sce) - t(as.matrix(XZ %*% coefs))

    } else {
        ## Fixed-effects model
        X <- model.matrix(msqrob_formula, colData(sce))
        X[, grepl(covariate, colnames(X))] <- 0
        X[, 1] <- 1

        full_coef <- colnames(X)

        coefs <- vapply(msqrob_models, function(m) {
            cfs <- getCoef(m)
            ## Pad missing coefficients with zero for matrix multiplication
            if (all(is.na(cfs))) return(rep(NA_real_, length(full_coef)))
            cfs[setdiff(full_coef, names(cfs))] <- 0
            cfs[full_coef]
        }, FUN.VALUE = numeric(length(full_coef)))

        assay(sce) <- assay(sce) - t(X %*% coefs)
    }

    sce
}
