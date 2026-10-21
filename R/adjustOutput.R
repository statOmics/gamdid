#' @importFrom stats p.adjust
#' @importFrom stageR stageR stageWiseAdjustment getAdjustedPValues

adjustOutput <- function(pval, pval_shape, pval_location,
    adjust = "stagewise", alpha) {
    if (adjust == "contrastwise") {
        p_adj <- adjustColumns(pval)
    } else if (adjust == "stagewise") {
        if (ncol(pval) < 2) {
            p_adj <- adjustColumns(pval)
        } else {
            p_adj <- stagewiseAdjust(pval, alpha)
        }
    } else if (adjust == "none") {
        colnames(pval_shape) <- paste0(colnames(pval_shape), "_GAM")
        colnames(pval_location) <- paste0(colnames(pval_location), "_msqrob2")
        colnames(pval) <- paste0(colnames(pval), "_aggregated")
        p_adj <- cbind(pval_shape, pval_location, pval)
    }

    threshold <- bhThreshold(pval, adjust, alpha)
    list(
        pvals = p_adj, adj_alpha = threshold$adj_alpha,
        n_tests = threshold$n_tests
    )
}

## Stagewise testing with stageR: omnibus screening, then Holm correction of
## the pairwise comparisons within each significant feature
stagewiseAdjust <- function(pval, alpha) {
    sr <- stageR(
        pScreen = setNames(pval[, "omnibus"], rownames(pval)),
        pConfirmation = pval[, -1, drop = FALSE]
    )
    sr_adj <- stageWiseAdjustment(sr,
        method = "holm", alpha = alpha, allowNA = TRUE
    )
    p_adj_raw <- getAdjustedPValues(sr_adj,
        onlySignificantGenes = FALSE, order = FALSE
    )
    ## preserve NA p-values that stageR drops for non-significant
    ## features
    p_adj <- matrix(NA_real_,
        nrow = nrow(pval), ncol = ncol(p_adj_raw),
        dimnames = list(rownames(pval), colnames(p_adj_raw))
    )
    p_adj[rownames(p_adj_raw), ] <- p_adj_raw
    colnames(p_adj)[1] <- "omnibus"
    p_adj
}

## BH threshold on the raw p values of the comparison (two groups) or the
## omnibus test: the largest k / m * alpha with p_(k) <= k / m * alpha
bhThreshold <- function(pval, adjust, alpha) {
    if (ncol(pval) < 2) {
        pval_vec <- pval
    } else {
        if (adjust == "none") {
            pval_vec <- pval[, "omnibus_aggregated"]
        } else {
            pval_vec <- pval[, "omnibus"]
        }
    }
    m <- length(pval_vec[!is.na(pval_vec)])
    p_sorted <- sort(pval_vec)
    thresholds <- seq_len(m) / m * alpha
    k <- max(which(p_sorted <= thresholds), 0)
    adj_alpha <- if (k > 0) thresholds[k] else NA
    list(adj_alpha = adj_alpha, n_tests = m)
}

## Benjamini-Hochberg correction of each column; assigning into p keeps the
## matrix when there is a single feature
adjustColumns <- function(p) {
    p[] <- apply(p, 2, p.adjust, method = "fdr")
    p
}
