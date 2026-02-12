setMethod(
  "visFit", "SummarizedExperiment",
  function(object,
           feature,
           plot = "hist",
           shade = TRUE,
           ...){
    visFitSCE(sce = object, feature = feature, plot = plot, shade = shade, ...)
  }
)


setMethod(
  "visFit", "QFeatures",
  function(object,
           i,
           feature,
           plot = "hist",
           shade = TRUE,
           ...){
    sce <- getWithColData(object, i)
    visFitSCE(sce = sce, feature = feature, plot = plot, shade = shade, ...)
  }
)
