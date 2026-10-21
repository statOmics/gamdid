#' @importFrom harmonicmeanp hmp.stat

aggregatePvals <- function(pval_shape, pval_location) {
    ## Both tests must report the same contrasts, otherwise contrasts would
    ## silently be dropped
    stopifnot(
        "shape and location tests report different contrasts" =
            setequal(colnames(pval_shape), colnames(pval_location))
    )
    common <- colnames(pval_shape)
    n_features <- nrow(pval_shape)

    pval <- vapply(common, function(nm) {
        apply(cbind(pval_shape[, nm], pval_location[, nm]), 1, hmp.stat)
    }, FUN.VALUE = numeric(n_features))

    ## keep a matrix when there is a single feature
    matrix(pval,
        nrow = n_features,
        dimnames = list(rownames(pval_shape), common)
    )
}
