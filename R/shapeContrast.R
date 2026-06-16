#' @importFrom utils combn

shapeContrast <- function(n_levels, nbasis, coef_names) {
    comps <- combn(n_levels, 2) - 1
    colnames(comps) <- apply(comps, 2, function(x) {
        paste0(coef_names[x[1] + 1], "_VS_", coef_names[x[2] + 1])
    })
    n_comps <- ncol(comps)

    p <- n_levels + n_levels * nbasis
    L <- matrix(0, nrow = p, ncol = n_comps * nbasis)

    for (i in seq_len(n_comps)) {
        a <- comps[, i]
        for (j in seq_len(nbasis)) {
            L[j + n_levels + nbasis * a[1], j + (i - 1) * nbasis] <- -1
            L[j + n_levels + nbasis * a[2], j + (i - 1) * nbasis] <-  1
        }
    }

    list(L = L, comps = comps)
}
