#' @importFrom stats pchisq

waldTest <- function(fit, L, df, tol = 1e-8) {
    coef <- fit$coefficients
    Vp <- fit$Vp

    est <- t(L) %*% coef
    Sigma <- t(L) %*% Vp %*% L

    eig <- eigen(Sigma, symmetric = TRUE)
    r <- sum(eig$values > tol)

    ## No testable contrast left: same output format as below, NA p-value
    if (r == 0) {
        return(c(pval = NA_real_))
    }

    halfCovInv <- eig$vectors[, seq_len(r), drop = FALSE] %*%
        diag(1 / sqrt(eig$values[seq_len(r)]), nrow = r)
    halfStat <- t(est) %*% halfCovInv
    stat <- crossprod(t(halfStat))
    pval <- pchisq(stat, df = df, lower.tail = FALSE)

    c("pval" = pval)
}
