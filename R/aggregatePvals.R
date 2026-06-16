#' @importFrom harmonicmeanp hmp.stat

aggregatePvals <- function(pval_shape, pval_location) {
    common <- intersect(colnames(pval_shape), colnames(pval_location))
    n_features <- nrow(pval_shape)

    pval <- vapply(common, function(nm) {
        apply(cbind(pval_shape[, nm], pval_location[, nm]), 1, hmp.stat)
    }, FUN.VALUE = numeric(n_features))

    rownames(pval) <- rownames(pval_shape)
    pval
}
