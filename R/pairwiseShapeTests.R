pairwiseShapeTests <- function(fit, L, comps, nbasis, tol = tol) {
    edf <- summary(fit)$edf
    n_comps <- ncol(comps)
    pvals <- numeric(n_comps)

    for (q in seq_len(n_comps)) {
        idx <- seq_len(nbasis) + (q - 1) * nbasis
        Lq <- L[, idx, drop = FALSE]
        lev <- comps[, q]
        df_q <- max(edf[lev])
        res <- waldTest(fit, Lq, df = df_q, tol = tol)
        pvals[q] <- res["pval"]
        names(pvals)[q] <- colnames(comps)[q]
    }
    pvals
}
