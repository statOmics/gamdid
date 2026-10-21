test_that("compDistr checks its arguments", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    for (adjust in list("BH", c("stagewise", "none"), NA)) {
        expect_error(
            compDistr(sce, covariate = "SampleType", adjust = adjust),
            "'adjust' must be one of"
        )
    }
    expect_error(
        compDistr(sce, covariate = "Sampletype"),
        "'covariate' must be the name of one column of colData"
    )
    expect_error(
        compDistr(sce, covariate = c("SampleType", "Set")),
        "'covariate' must be the name of one column of colData"
    )
    expect_error(
        compDistr(sce, covariate = 1),
        "'covariate' must be the name of one column of colData"
    )
    for (bins in list(0, 9, 12.5, "25", c(10, 20), NA)) {
        expect_error(
            compDistr(sce, covariate = "SampleType", bins = bins),
            "'bins' must be a single whole number of at least 10"
        )
    }
    for (minObs in list(-5, "50", c(10, 20), NA)) {
        expect_error(
            compDistr(sce, covariate = "SampleType", minObs = minObs),
            "'minObs' must be a single non-negative number"
        )
    }
})

test_that("compDistr accepts valid boundary values", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    expect_no_error(compDistr(sce,
        covariate = "SampleType", bins = 10,
        minObs = 0, adjust = "none"
    ))
})

test_that("visFit checks the plot type", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    out <- compDistr(sce, covariate = "SampleType")
    expect_error(
        visFit(out, feature = "E9PAV3", plot = "histogram"),
        "'plot' must be \"hist\" or \"diff\""
    )
})
