#' @importFrom msqrob2 getContrast getFitMethod getCoef getVcovUnscaled
#' @importFrom msqrob2 getVarPosterior getDfPosterior
#' @importFrom stats na.exclude pf

topFeaturesOmnibus <- function(models, contrasts, adjust.method = "BH",
    sort = TRUE, alpha = 1) {
    qrL <- qr(contrasts)
    Lh <- qr.Q(qrL)[, seq_len(qrL$rank)]
    logFCs <- apply(contrasts, 2, function(models, L) {
        vapply(models, getContrast, L = L, FUN.VALUE = numeric(1))
    }, models = models)
    ## keep one row per model when there is a single model
    logFCs <- matrix(logFCs,
        nrow = length(models),
        dimnames = list(names(models), colnames(contrasts))
    )

    Fstat <- omnibusFstat(models, rownames(contrasts), Lh, qrL$rank)

    df <- vapply(models, getDfPosterior, FUN.VALUE = numeric(1))
    pval <- pf(Fstat, qrL$rank, df, lower.tail = FALSE)
    adjPval <- p.adjust(pval, method = adjust.method)
    out <- data.frame(logFCs, df1 = qrL$rank, df2 = df, Fstat, pval, adjPval)

    if (sort) {
        if (alpha < 1) {
            ids <- adjPval < alpha
            out <- na.exclude(out[ids, ])
        }
        return(out[order(out$pval), ])
    }
    out
}

## Omnibus F statistic of each model for the contrasts Lh of rank r
omnibusFstat <- function(models, full_coef, Lh, r) {
    full_varcov <- matrix(0,
        nrow = length(full_coef), ncol = length(full_coef),
        dimnames = list(full_coef, full_coef)
    )

    vapply(models, function(model, Lh, r) {
        if (getFitMethod(model) == "fitError") {
            return(NA_real_)
        }
        cfs <- getCoef(model)
        present_cfs <- intersect(names(cfs), full_coef)
        varcov <- full_varcov
        ## Fill in variance-covariance for present coefficients only
        varcov[present_cfs, present_cfs] <- as.matrix(getVcovUnscaled(model))
        cfs[setdiff(full_coef, names(cfs))] <- 0
        cfs <- cfs[full_coef]
        est <- cfs %*% Lh
        varEst <- t(Lh) %*% (varcov * getVarPosterior(model)) %*% Lh
        if (inherits(try(solve(varEst), silent = TRUE), "try-error")) {
            return(NA_real_)
        }
        as.numeric((est %*% solve(varEst) %*% t(est)) / r)
    }, Lh = Lh, r = r, FUN.VALUE = numeric(1))
}
