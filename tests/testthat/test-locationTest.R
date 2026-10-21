test_that("returns a list", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
    expect_true(is(out, "list"))
})

test_that("list with pvals", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
    expect_identical(names(out), "pvals")
})

test_that("pvals object is a matrix", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
    expect_true(is(out$pvals, "matrix"))
})

test_that("list with pvals matrix, msqrob_models en msqrob_formula", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~ (1 | Set))
    expect_identical(names(out), c("pvals", "msqrob_models", "msqrob_formula"))
})

test_that("list with pvals matrix, msqrob_models en msqrob_formula (2)", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~ -1 + Set)
    expect_identical(names(out), c("pvals", "msqrob_models", "msqrob_formula"))
})

test_that("pvals object is a matrix (for nuisance correction)", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~ -1 + (1 | Set))
    expect_true(is(out$pvals, "matrix"))
})

test_that("pvals object is a matrix (for nuisance correction 2)", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~Set)
    expect_true(is(out$pvals, "matrix"))
})

test_that("msqrobModels contain StatModel objects", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~Set)
    expect_true(is(out$msqrob_models[[1]], "StatModel"))
})

test_that("msqrob_formula contain formula objects", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~Set)
    expect_true(is(out$msqrob_formula, "formula"))
})

test_that("compute correct p values", {
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~ (1 | Set))
    expect_equal(out$pvals, structure(c(
        0.237037266770739, 0.981984915759283, 0.996259580193254,
        0.000181307193181744, 0.344405576858417, 0.000708840124653451,
        3.35828989214574e-08, 0.226269016259427, 0.235860756436294, 4.9922259329987e-05
    ), dim = c(10L, 1L), dimnames = list(NULL, "SampleTypeMelanoma_VS_SampleTypeMonocyte")), tolerance = 1e-6)
})

test_that("correct dimensions", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = ~ (1 | Set))
    expect_equal(dim(out$pvals), c(10, ncol(combn(nlevels(leduc_testdata$Field), 2)) + 1))
})


test_that("features without observations get NA p values", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    assay(sce)[1, ] <- NA
    out <- locationTest(sce = sce, covariate = "SampleType", nuisance = ~ (1 | Set))
    ref <- locationTest(sce = sce[-1, ], covariate = "SampleType", nuisance = ~ (1 | Set))
    expect_true(all(is.na(out$pvals[1, ])))
    expect_equal(out$pvals[-1, , drop = FALSE], ref$pvals)
    expect_identical(names(out$msqrob_models), rownames(sce))
    expect_identical(msqrob2::getFitMethod(out$msqrob_models[[1]]), "fitError")
})

test_that("features without observations get NA p values (multiple groups)", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    assay(sce)[1, ] <- NA
    out <- locationTest(sce = sce, covariate = "Field", nuisance = ~ (1 | Set))
    ref <- locationTest(sce = sce[-1, ], covariate = "Field", nuisance = ~ (1 | Set))
    expect_true(all(is.na(out$pvals[1, ])))
    expect_equal(out$pvals[-1, , drop = FALSE], ref$pvals)
})

test_that("error when all features contain only missing values", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    assay(sce)[] <- NA
    expect_error(
        locationTest(sce = sce, covariate = "SampleType", nuisance = NULL),
        "only missing values"
    )
})
