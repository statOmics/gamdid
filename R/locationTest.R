#' @importFrom stats as.formula model.matrix
#' @importFrom msqrob2 msqrob hypothesisTest makeContrast StatModel
#' @importFrom SummarizedExperiment colData rowData assay
#' @importFrom reformulas findbars
#' @importFrom utils combn
#' @include mixedDesign.R

locationTest <- function(sce, covariate, nuisance, ...) {
    ## msqrob2 errors on features without any observation (mixed models),
    ## so only fit the observed features and return NA for the others
    keep <- rowSums(!is.na(assay(sce))) > 0
    if (!any(keep)) {
        stop("all features contain only missing values")
    }
    if (!all(keep)) {
        out <- locationTest(sce[keep, ], covariate, nuisance, ...)
        return(expandLocation(out, keep, rownames(sce)))
    }

    msqrob_formula <- msqrobFormula(covariate, nuisance)

    ## Fit location models
    pe <- msqrob(sce, formula = msqrob_formula, ...)

    coef_names <- coefNames(msqrob_formula, sce)

    ## Identify coefficients belonging to the covariate of interest: exact
    ## names, so that e.g. a nuisance "SetGroup" is not taken for "Group"
    levels <- unique(colData(sce)[[covariate]])
    coef <- coef_names[coef_names %in% paste0(covariate, levels)]
    colnames_L <- levelNames(coef_names, covariate, levels)

    ## Set up contrasts
    if (length(levels) < 3) {
        pvals <- twoGroupPvals(pe, coef, coef_names, colnames_L)
    } else {
        pvals <- multiGroupPvals(
            pe, coef_names, colnames_L, length(levels), nrow(sce)
        )
    }

    ## Return list; include model objects only when nuisance correction is
    ## needed
    if (is.null(nuisance)) {
        list(pvals = pvals)
    } else {
        list(
            pvals = pvals,
            msqrob_models = rowData(pe)$msqrobModels,
            msqrob_formula = msqrob_formula
        )
    }
}

## msqrob formula: covariate of interest plus the nuisance terms
msqrobFormula <- function(covariate, nuisance) {
    if (is.null(nuisance)) {
        msqrob_formula <- as.formula(paste0("~", covariate))
    } else {
        f <- as.character(nuisance)
        msqrob_formula <- as.formula(paste0(f[1], covariate, "+", f[2]))
    }
    msqrob_formula
}

## Coefficient names of the msqrob model, in the correct order
coefNames <- function(msqrob_formula, sce) {
    if (!is.null(findbars(msqrob_formula))) {
        design <- mixedDesign(msqrob_formula, colData(sce))
        coef_names <- c(colnames(design$X), design$colnames_Z)
    } else {
        coef_names <- colnames(model.matrix(msqrob_formula, colData(sce)))
    }
    coef_names
}

## Coefficient names with "(Intercept)" replaced by the reference level name
levelNames <- function(coef_names, covariate, levels) {
    colnames_L <- coef_names
    coef_names_full <- paste0(covariate, levels)
    if (!all(coef_names_full %in% coef_names)) {
        colnames_L[colnames_L == "(Intercept)"] <-
            coef_names_full[which(!(coef_names_full %in% coef_names))]
    }
    colnames_L
}

## Two groups: single pairwise comparison
twoGroupPvals <- function(pe, coef, coef_names, colnames_L) {
    L <- if (length(coef) < 2) {
        makeContrast(paste0(coef, "=0"), parameterNames = coef_names)
    } else {
        makeContrast(
            paste0(coef[1], "-", coef[2], "=0"),
            parameterNames = coef_names
        )
    }

    pe <- hypothesisTest(pe,
        contrast = L, adjust.method = "none", overwrite = TRUE
    )
    pvals <- matrix(rowData(pe)[[colnames(L)]]$pval, ncol = 1)
    colnames(pvals) <- paste0(colnames_L[1], "_VS_", colnames_L[2])
    pvals
}

## More than two groups: omnibus test + pairwise comparisons
multiGroupPvals <- function(pe, coef_names, colnames_L, n_levels,
                            n_features) {
    comps <- combn(n_levels, 2)
    L <- matrix(0, length(coef_names), ncol(comps),
        dimnames = list(coef_names, NULL)
    )
    colnames(L) <- apply(comps, 2, function(x) {
        paste0(colnames_L[x[1]], "_VS_", colnames_L[x[2]])
    })

    for (i in seq_len(ncol(comps))) {
        L[comps[1, i], i] <- -1
        L[comps[2, i], i] <- 1
    }
    if ("(Intercept)" %in% coef_names) L["(Intercept)", ] <- 0

    p_global <- topFeaturesOmnibus(
        models = rowData(pe)$msqrobModels,
        contrasts = L,
        adjust.method = "none",
        sort = FALSE
    )$pval

    p_pc <- vapply(seq_len(ncol(L)), function(i) {
        pe_pc <- hypothesisTest(pe,
            contrast = L[, i, drop = FALSE],
            adjust.method = "none", overwrite = TRUE
        )
        rowData(pe_pc)[[colnames(L)[i]]]$pval
    }, FUN.VALUE = numeric(n_features))
    ## keep a matrix when there is a single feature
    p_pc <- matrix(p_pc,
        nrow = n_features, dimnames = list(NULL, colnames(L))
    )

    cbind(omnibus = p_global, p_pc)
}

## Expand locationTest output for a subset of features (keep) back to all
## features: NA p-values and fitError models for the features not fitted
expandLocation <- function(out, keep, feature_names) {
    pvals <- matrix(NA_real_,
        nrow = length(keep), ncol = ncol(out$pvals),
        dimnames = list(NULL, colnames(out$pvals))
    )
    pvals[keep, ] <- out$pvals
    out$pvals <- pvals

    if (!is.null(out$msqrob_models)) {
        models <- rep(list(StatModel(
            type = "fitError",
            params = list(
                coefficients = NA, vcovUnscaled = NA, sigma = NA,
                df.residual = NA, w = NA
            ),
            varPosterior = NA_real_, dfPosterior = NA_real_
        )), length(keep))
        names(models) <- feature_names
        models[keep] <- out$msqrob_models
        out$msqrob_models <- models
    }
    out
}
