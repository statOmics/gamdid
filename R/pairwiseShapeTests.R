pairwiseShapeTests <- function(fit, L, comps, nbasis, tol = tol) {
  edf <- summary(fit)$edf
  coef_names <- names(fit$coefficients)
  
  n_comps <- ncol(comps)
  pvals <- numeric(n_comps)
  
  # loop over all possible pairwise comparisons
  for (q in seq_len(n_comps)) {
    idx <- ((q - 1) * nbasis + 1):(q * nbasis)
    Lq <- L[, idx, drop = FALSE]
    lev <- comps[, q]
    
    df_q <- max(edf[lev])
    
    res <- waldTest(fit, Lq, df = df_q, tol = tol)
    pvals[q] <- res["pval"]
    
    names(pvals)[q] <- colnames(comps)[q]
  }
  pvals
}