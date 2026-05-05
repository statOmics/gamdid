partialResiduals <- function(sce, covariate, msqrob_models, msqrob_formula) {
  # Fixed or mixed model?
  if (!is.null(lme4::findbars(msqrob_formula))) {
    # mock example to create model matrix
    mock <- data.frame(colData(sce), y = 0)
    mock_fit <- lme4::lmer(paste0("y~",as.character(msqrob_formula)[2]), data = mock)

    X <- getME(mock_fit, "X")
    Z <- getME(mock_fit, "Z")
    ## remove the covariate effect
    X[, grepl(covariate, colnames(X))] <- 0
    X[, 1] <- 1

    XZ <- cbind2(X, Z)

    # Match coefficient names to msqrob output
    colnames_Z <- paste0("(Intercept)",
                         unlist(lapply(names(ranef(mock_fit)),
                                       function(nm){
                                         df <- ranef(mock_fit)[[nm]]
                                         paste0(nm, rownames(df))
                                         }
                                       )
                                )
                         )
    full_coef <- c(colnames(X), colnames_Z)

    # retrieve coefficients from msqrob
    coefs <- sapply(msqrob_models, function(m) {
      cfs <- getCoef(m)
      ## fill out missing coefficients to allow matrix multiplication
      if (all(is.na(cfs))) return(rep(NA, length(full_coef)))
      cfs[setdiff(full_coef, names(cfs))] <- 0
      cfs[full_coef]
    })

    test <- sapply(msqrob_models, function(m) {
      cfs <- getCoef(m)})

    # create partial residuals
    assay(sce) <- assay(sce) - t(as.matrix(XZ %*% coefs))

  } else {
    # create model matrix
    X <- model.matrix(msqrob_formula, colData(sce))
    ## remove covariate effect
    X[, grepl(covariate, colnames(X))] <- 0
    X[, 1] <- 1

    full_coef <- colnames(X)

    # retrieve coefficients from msqrob
    coefs <- sapply(msqrob_models, function(m) {
      cfs <- getCoef(m)
      ## fill out missing coefficients to allow matrix multiplication
      if (all(is.na(cfs))) return(rep(NA, length(full_coef)))
      cfs[setdiff(full_coef, names(cfs))] <- 0
      cfs[full_coef]
    })

    # create partial residuals
    assay(sce) <- assay(sce) - t(X %*% coefs)
  }

  sce
}
