## Half-width of the confidence band of the difference smoother
bandHalfWidth <- function(p) {
    band <- Filter(function(l) "lwr" %in% names(l$data), p$layers)[[1]]$data
    band$diff - band$lwr
}

test_that("all comparisons use the same alpha", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    out <- compDistr(object = sce, covariate = "Field", minObs = 20)
    ## The half-width is a critical value (depending on alpha only) times a
    ## standard error (independent of alpha), so the ratio of the half-widths
    ## at two alpha levels must be the same for every comparison
    p1 <- visFit(out, feature = "E9PAV3", plot = "diff", alpha = 0.01)
    p2 <- visFit(out, feature = "E9PAV3", plot = "diff", alpha = 0.1)
    ratios <- unlist(lapply(names(p1), function(comp) {
        bandHalfWidth(p1[[comp]]) / bandHalfWidth(p2[[comp]])
    }))
    expect_length(p1, 6)
    expect_equal(ratios, rep(ratios[1], length(ratios)))
})

test_that("alpha defaults to the adjusted threshold of compDistr", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    out <- compDistr(object = sce, covariate = "Field", minObs = 20)
    p_default <- visFit(out, feature = "E9PAV3", plot = "diff")
    p_explicit <- visFit(out,
        feature = "E9PAV3", plot = "diff",
        alpha = metadata(out)$adj_alpha
    )
    for (comp in names(p_default)) {
        expect_equal(
            bandHalfWidth(p_default[[comp]]),
            bandHalfWidth(p_explicit[[comp]])
        )
    }
})

test_that("compDistr stores alpha and the number of tests", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    out <- compDistr(
        object = sce, covariate = "SampleType", adjust = "none",
        alpha = 0.1
    )
    p <- rowData(out)$SampleTypeMelanoma_VS_SampleTypeMonocyte_aggregated
    expect_identical(metadata(out)$alpha, 0.1)
    expect_identical(metadata(out)$n_tests, sum(!is.na(p)))
})

test_that("without significant features the most stringent BH level is used", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    ## Nothing is significant at this FDR level
    out <- compDistr(object = sce, covariate = "SampleType", alpha = 1e-4)
    expect_true(is.na(metadata(out)$adj_alpha))
    expect_message(
        p_default <- visFit(out, feature = "E9PAV3", plot = "diff"),
        "most stringent BH level"
    )
    p_explicit <- visFit(out,
        feature = "E9PAV3", plot = "diff",
        alpha = 1e-4 / metadata(out)$n_tests
    )
    hw <- bandHalfWidth(p_default)
    expect_false(anyNA(hw))
    expect_equal(hw, bandHalfWidth(p_explicit))
})

test_that("error when the BH level cannot be derived", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    out <- compDistr(object = sce, covariate = "SampleType", alpha = 1e-4)
    ## Object created before alpha and n_tests were stored
    metadata(out)$alpha <- NULL
    metadata(out)$n_tests <- NULL
    expect_error(visFit(out, feature = "E9PAV3"), "please specify 'alpha'")
    ## An explicit alpha still works
    expect_no_error(visFit(out, feature = "E9PAV3", alpha = 0.05))
})
