#' @importFrom stats as.formula model.matrix
#' @importFrom msqrob2 msqrob hypothesisTest makeContrast
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom lme4 ranef
#' @importFrom reformulas findbars
#' @importFrom utils combn


locationTest <- function(sce, covariate, nuisance, ...) {
  # Create the msqrob formula
  if (is.null(nuisance)) {
    msqrob_formula <- as.formula(paste0("~", covariate))
  } else {
    f <- as.character(nuisance)
    msqrob_formula <- as.formula(paste0(f[1], covariate, "+", f[2]))
  }

  # fit the models
  pe <- msqrob(sce, formula = msqrob_formula, ...)
  # retrieve all coefficient names in correct order
  if (!is.null(findbars(msqrob_formula))){
    mock <- data.frame(colData(sce), y = 0)
    mock_fit <- lme4::lmer(paste0("y~",as.character(msqrob_formula)[2]), data = mock)
    X <- lme4::getME(mock_fit, "X")
    Z <- lme4::getME(mock_fit, "Z")
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
    coef_names <- c(colnames(X), colnames_Z)
  } else {
    coef_names <- colnames(model.matrix(msqrob_formula, colData(sce)))
  }

  ## retrieve all the covariate coefficients
  coef <- coef_names[grepl(covariate, coef_names)]

  # Set up contrasts
  levels <- unique(colData(sce)[[covariate]])
  if (length(levels) < 3) {
    # 2 groups: 1 comparison
    L <- if (length(coef) < 2) {
      makeContrast(paste0(coef, "=0"), parameterNames = coef_names)
    } else {
      makeContrast(paste0(coef[1], "-", coef[2], "=0"),
                   parameterNames = coef_names)
    }

    colnames_L <- coef_names
    # replace "(Intercept)" with reference class name
    coef_names_full <- paste0(covariate, levels)
    if (!all(coef_names_full %in% coef_names)){
      colnames_L[colnames_L == "(Intercept)"] <-
        coef_names_full[which(!(coef_names_full %in% coef_names))]
    }

    pe <- hypothesisTest(pe, contrast = L, adjust.method = "none", overwrite = TRUE)
    pvals <- matrix(rowData(pe)[[colnames(L)]]$pval, ncol = 1)
    colnames(pvals) <- paste0(colnames_L[1], "_VS_", colnames_L[2])

  } else {
    # > 2 groups: global test and pairwise comparisons
    comps <- combn(length(levels), 2)
    L <- matrix(0, length(coef_names), ncol(comps),
                dimnames = list(coef_names, NULL))

    colnames_L <- coef_names
    # replace "(Intercept)" with reference class name
    coef_names_full <- paste0(covariate, levels)
    if (!all(coef_names_full %in% coef_names)){
      colnames_L[colnames_L == "(Intercept)"] <-
        coef_names_full[which(!(coef_names_full %in% coef_names))]
    }
    colnames(L) <- apply(comps, 2, function(x)
      paste0(colnames_L[x[1]], "_VS_", colnames_L[x[2]])
      )

    for (i in seq_len(ncol(comps))) {
      L[comps[1, i], i] <- -1
      L[comps[2, i], i] <-  1
    }

    if ("(Intercept)" %in% coef_names) L["(Intercept)", ] <- 0

    p_global <- topFeaturesOmnibus(
      models = rowData(pe)$msqrobModels,
      contrasts = L,
      adjust.method = "none",
      sort = FALSE
    )$pval

    p_pc <- sapply(seq_len(ncol(L)), function(i) {
      pe_pc <- hypothesisTest(pe, contrast = L[, i, drop = FALSE],
                              adjust.method = "none", overwrite = TRUE)
      rowData(pe_pc)[[colnames(L)[i]]]$pval
    })
    colnames(p_pc) <- colnames(L)

    pvals <- cbind(omnibus = p_global, p_pc)
  }

  # create the output
  if (is.null(nuisance)){
    out <- list(pvals = pvals)
  } else {
    out <- list(pvals = pvals, msqrob_models = rowData(pe)$msqrobModels, msqrob_formula = msqrob_formula)
  }
}
