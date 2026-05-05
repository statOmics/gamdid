#' @importFrom graphics hist
#' @importFrom utils globalVariables

utils::globalVariables(c("counts", "n"))

createCounts <- function(x, colDat, covariate, minObs = 50, bins = 25){

  if (all(is.na(x))) {
    return(NA)
  }

  cov_values <- colDat[[covariate]]
  cov_levels <- unique(cov_values)

  # Check if at least minObs in each group
  obsPerLvl <- tapply(!is.na(x), colDat[[covariate]], sum)
  if (! all(obsPerLvl > minObs)) return(NA)

  # Define histogram grid
  min_x <- floor(min(x, na.rm = TRUE) * 10) / 10
  max_x <- ceiling(max(x, na.rm = TRUE) * 10) / 10
  binwidth <- (max_x - min_x) / bins
  breaks <- seq(min_x, max_x, by = binwidth)

  # Define histogram
  histograms <- lapply(cov_levels, function(lvl) {
    hist(
      x[cov_values == lvl],
      breaks = breaks,
      plot = FALSE
    )
  })

  # generate output
  out <- data.frame(
    midsC = unlist(lapply(histograms, `[[`, "mids")),
    counts = unlist(lapply(histograms, `[[`, "counts")),
    covariate = rep(cov_levels, each = length(breaks) - 1)
  )

  out |>
    dplyr::group_by(covariate) |>
    dplyr::mutate(
      n = sum(counts),
      covariate = as.factor(covariate),
      carrier = log(n * binwidth)
    ) |>
    dplyr::ungroup()
}
