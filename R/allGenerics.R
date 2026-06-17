#' @title Compare distributions across groups
#'
#' @export
setGeneric(
  "compDistr",
  function(object, ...) standardGeneric("compDistr")
)

#' @title Visualize GAM fits for a single feature
#'
#' @export
setGeneric(
    "visFit",
    function(object, ...) standardGeneric("visFit")
)

#' @title Retrieve top-ranked features after differential distribution analysis
#'
#' @export
setGeneric(
    "topFeatures",
    function(object, ...) standardGeneric("topFeatures")
)
