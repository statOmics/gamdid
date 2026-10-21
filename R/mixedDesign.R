#' @importFrom lme4 lFormula
#' @importFrom stats as.formula

## Fixed (X) and random (Z) effects design matrices of a mixed model formula
## with random intercepts, and the names msqrob2 gives to the random
## intercepts. lFormula() builds the matrices without fitting a model.
mixedDesign <- function(formula, colDat) {
    mock <- data.frame(colDat, y = 0)
    lf <- lFormula(as.formula(paste0("y~", as.character(formula)[2])),
        data = mock
    )
    fl <- lf$reTrms$flist
    colnames_Z <- paste0(
        "(Intercept)",
        unlist(lapply(names(fl), function(nm) paste0(nm, levels(fl[[nm]]))))
    )
    list(X = lf$X, Z = Matrix::t(lf$reTrms$Zt), colnames_Z = colnames_Z)
}
