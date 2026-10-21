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
        NA, NA, NA, NA, NA, NA, 8.3575303343265779e-09, 0.44860621529663919,
        0.018185656766666329, 0.26354290748796699, NA, NA, NA, NA, NA,
        NA, 7.2688920685479539e-10, 0.13626464563042967, 0.79451456257764164,
        0.84809240948143216, NA, NA, NA, NA, NA, NA, 4.2353590614914295e-05,
        0.39096330978667698, 0.0046454554329421307, 0.17502843412285474,
        NA, NA, NA, NA, NA, NA, 0.0015212631745834093, 0.25939518256584998,
        0.18612667034540964, 0.73094821669126109, NA, NA, NA, NA, NA,
        NA, 0.063837726581496307, 0.54837491507015745, 0.013225102278790028,
        0.052085340458161059, NA, NA, NA, NA, NA, NA, 0.011567197794244883,
        0.51228392531882438, 0.20875470952186209, 0.67459620961723954,
        NA, NA, NA, NA, NA, NA, 0.3999045693317706, 0.59160411842514971,
        0.019195402618452102, 0.32222652996441853
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
