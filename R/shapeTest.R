#' @include shapeContrast.R waldTest.R
NULL

shapeTest <- function(sce, covariate, minObs = 50, bins = 25, tol = 1e-8){
  colDat <- colData(sce)
  X <- assay(sce)

  cov_values <- colDat[[covariate]]
  cov_levels <- unique(cov_values)
  n_levels <- length(cov_levels)

  # ------------------------------------------------------------------
  # Histogram construction per feature
  # ------------------------------------------------------------------
  hist_dat <- apply(
    X,
    MARGIN = 1,
    FUN = createCounts,
    colDat = colDat,
    covariate = covariate,
    minObs = minObs,
    bins = bins
  )

  # ------------------------------------------------------------------
  # GAM fitting
  # ------------------------------------------------------------------
  model_list <- lapply(hist_dat, function(dat) {
    if (all(is.na(dat)) || any(is.infinite(dat$carrier))) {
      return(NA)
    }

    mgcv::gam(
      counts ~ covariate + s(midsC, by = covariate) + offset(carrier),
      data = dat,
      family = poisson(),
      method = "REML"
    )
  })

  # --------------------------------------------------------------
  # Build contrast (from first valid model)
  # --------------------------------------------------------------
  idx <- which(vapply(model_list, function(x) inherits(x, "gam"), logical(1)))[1]
  ref_fit <- model_list[[idx]]

  if (is.null(ref_fit)) {
    stop("No valid GAM fits available to perform inference.")
  }

  p <- length(ref_fit$coefficients)
  nbasis <- (p - n_levels) / n_levels
  coef_names <- names(ref_fit$coefficients[1:n_levels])
  # replace "(Intercept)" with reference class name
  coef_names_full <- paste0("covariate", cov_levels)
  if (!all(coef_names_full %in% coef_names)){
    coef_names[coef_names == "(Intercept)"] <-
      coef_names_full[which(!(coef_names_full %in% coef_names))]
  }
  ctr <- shapeContrast(n_levels, nbasis, coef_names)



  # ------------------------------------------------------------------
  # Inference
  # ------------------------------------------------------------------
  if (ncol(ctr$comps) < 2){
    pvals <- vapply(model_list, function(fit) {
      if (all(is.na(fit))) return(NA)

      # Wald test
      edf <- summary(fit)$edf
      df_omni <- sum(sort(edf, decreasing = TRUE)[seq_len(n_levels - 1)])
      res <- waldTest(fit, ctr$L, df = df_omni, tol = tol)
      names(res) <- colnames(ctr$comps)
      res
    }, FUN.VALUE = numeric(1))
    pvals <- matrix(pvals, nrow = 1,dimnames = list(colnames(ctr$comps), NULL))
  } else {
    pvals <- vapply(model_list, function(fit){
      if (all(is.na(fit))) return(setNames(rep(NA_real_, ncol(ctr$comps)+1), c("omnibus", colnames(ctr$comps)) ))

      # Omnibus Wald test
      edf <- summary(fit)$edf
      df_omni <- sum(sort(edf, decreasing = TRUE)[seq_len(n_levels - 1)])
      res <- waldTest(fit, ctr$L, df = df_omni, tol = tol)

      # Pairwise comparisons
      names(res)[which(names(res) == "pval")] <- "omnibus"
      res <- c(res,
               pairwiseShapeTests(fit, ctr$L, ctr$comps, nbasis, tol = tol)
      )
      res
    }, FUN.VALUE = numeric(ncol(ctr$comps)+1))
  }

  # ------------------------------------------------------------------
  # Output harmonisation
  # ------------------------------------------------------------------
  if (is.matrix(pvals)) {
    pvals <- t(pvals)
  } #else {
  #   all_names <- unique(unlist(lapply(pvals, names)))
  #   pvals <- t(vapply(pvals, function(x) {
  #     out <- setNames(rep(NA_real_, length(all_names)), all_names)
  #     out[names(x)] <- x
  #     out
  #   }, numeric(length(all_names))))
  # }
  colnames(pvals) <- gsub("covariate",covariate,colnames(pvals))

  return(list(models = model_list, pvals = pvals))
}
