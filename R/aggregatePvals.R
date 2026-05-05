#' @importFrom harmonicmeanp hmp.stat

aggregatePvals <- function(pval_shape, pval_location) {
  common <- intersect(colnames(pval_shape), colnames(pval_location))

  pval <- sapply(common, function(nm) {
    apply(cbind(pval_shape[, nm], pval_location[, nm]), 1, hmp.stat)
  })

  rownames(pval) <- rownames(pval_shape)
  return(pval)
}
