#' @importFrom stats p.adjust
#' @importFrom stageR stageR stageWiseAdjustment getAdjustedPValues


adjustOutput <- function(pval, pval_shape, pval_location,
                          adjust = T, alpha) {
  if (adjust) {
    if (ncol(pval) < 2) {
      p_adj <- apply(pval, 2, p.adjust, method = "fdr")
    } else {
      sr <- stageR(pScreen = pval[, "omnibus"],
                   pConfirmation = pval[, -1])
      sr_adj <- stageWiseAdjustment(sr, method = "holm", alpha = alpha, allowNA = T)
      p_adj_raw <- getAdjustedPValues(sr_adj, onlySignificantGenes = FALSE, order = F)
      # keep the NA p-values that are removed
      p_adj <- matrix(NA_real_, nrow = nrow(pval), ncol = ncol(p_adj_raw),
                               dimnames = list(rownames(pval), colnames(p_adj_raw)))
      p_adj[rownames(p_adj_raw),] <- p_adj_raw
      colnames(p_adj)[1] <- "omnibus"
    }
  } else {
    colnames(pval_shape) <- paste0(colnames(pval_shape), "_shape")
    colnames(pval_location) <- paste0(colnames(pval_location), "_location")
    colnames(pval) <- paste0(colnames(pval), "_aggregated")
    p_adj <- cbind(pval_shape, pval_location, pval)
  }

  if (ncol(pval) < 2){
    adj_alpha <- alpha
  } else {
    # calculate the adjusted alpha level
    if (adjust){
      pval_vec <- pval[,"omnibus"]
    } else{
      pval_vec <- pval[,"omnibus_aggregated"]
    }
    m <- length(pval_vec[!is.na(pval_vec)])
    p_sorted <- sort(pval_vec)
    thresholds <- (1:m) / m * 0.05
    k <- max(which(p_sorted <= thresholds), 0)
    adj_alpha <- if (k > 0) thresholds[k] else NA
  }

  list(pvals = p_adj, adj_alpha = adj_alpha)
}
