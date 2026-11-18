test_that("compute a matrix", {
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    expect_true(is(out, "matrix"))
})
test_that("compute the harmonic mean of the shape and location p values", {
    pval_shape <- matrix(c(0.01, 0.5, NA),
        ncol = 1,
        dimnames = list(c("a", "b", "c"), "A_VS_B")
    )
    pval_location <- matrix(c(0.04, 0.5, 0.01),
        ncol = 1,
        dimnames = list(NULL, "A_VS_B")
    )
    out <- aggregatePvals(pval_shape = pval_shape, pval_location = pval_location)
    expect_equal(out, matrix(c(2 / (1 / 0.01 + 1 / 0.04), 0.5, NA),
        ncol = 1,
        dimnames = list(c("a", "b", "c"), "A_VS_B")
    ))
})

test_that("aggregate the p values of the shape and location tests", {
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~Set)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    expect_identical(dimnames(out), list(NULL, "SampleTypeMelanoma_VS_SampleTypeMonocyte"))
    expect_equal(out[, 1], apply(cbind(out_shape$pvals[, 1], out_location$pvals[, 1]), 1, hmp.stat),
        ignore_attr = TRUE
    )
})

test_that("compute the correct p values", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    expect_equal(out, structure(c(
        NA, NA, NA, NA, NA, NA, 1.7814572489290558e-08, 0.49898268700186771,
        0.027326600058983585, 0.26612527104841677, NA, NA, NA, NA, NA,
        NA, 6.8961171945156473e-09, 0.16440797663485862, 0.80173161956156258,
        0.84897218875173808, NA, NA, NA, NA, NA, NA, 7.6178909328397806e-05,
        0.41889939495096046, 0.0070196168389546317, 0.17886948118175355,
        NA, NA, NA, NA, NA, NA, 0.0023515755576724059, 0.29751497661733772,
        0.21593024174642575, 0.73193789213276772, NA, NA, NA, NA, NA,
        NA, 0.080284047356523322, 0.57210566035671317, 0.023544400096216594,
        0.052389537317230409, NA, NA, NA, NA, NA, NA, 0.016583016313481178,
        0.52092792245892039, 0.22794713477546721, 0.67602472708658867,
        NA, NA, NA, NA, NA, NA, 0.41002327173880804, 0.59670175623173194,
        0.019350120793173407, 0.33319969030909974
    ), dim = c(10L, 7L), dimnames = list(
        c(
            "A0AVT1", "A2A368", "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2",
            "E9PAV3", "E9PRG8", "O00154", "O00159"
        ), c(
            "omnibus", "Field1_VS_Field2",
            "Field1_VS_Field3", "Field1_VS_Field4", "Field2_VS_Field3",
            "Field2_VS_Field4", "Field3_VS_Field4"
        )
    )), tolerance = 1e-6)
})
