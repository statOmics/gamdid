test_that("returns a list", {
    out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    expect_true(is(out, "list"))
})

test_that("list with GAM models and pvals", {
    out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    expect_identical(names(out), c("models", "pvals"))
})

test_that("pvals object is a matrix", {
    out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    expect_true(is(out$pvals, "matrix"))
})

test_that("models contain StatModel objects", {
    out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    expect_true(is(out$models[[1]], "gam"))
})

test_that("compute correct p values", {
    out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    expect_equal(out$pvals, structure(
        c(
            0.0183805817026951, 0.000386393734274719, 0.5383450730135,
            0.00546027282729357, 0.595862779628328, 0.0150442734089543, 4.80544360393079e-05,
            0.0199681251740768, 0.335987280739658, 0.00189898572291178
        ),
        dim = c(10L, 1L), dimnames = list(NULL, "SampleTypeMelanoma_VS_SampleTypeMonocyte")
    ), tolerance = 1e-6)
})

test_that("correct dimensions", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
    expect_equal(dim(out$pvals), c(10, ncol(combn(nlevels(leduc_testdata$Field), 2)) + 1))
})


test_that("pairwise tests use the edf of the two compared groups", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    out <- shapeTest(sce = sce, covariate = "Field")
    fit <- out$models[["E9PAV3"]]
    edf <- summary(fit)$edf
    labels <- vapply(fit$smooth, `[[`, character(1), "label")
    lvls <- levels(sce$Field)
    nbasis <- (length(coef(fit)) - length(lvls)) / length(lvls)
    ctr <- shapeContrast(length(lvls), nbasis, paste0("covariate", lvls))
    for (q in seq_len(ncol(ctr$comps))) {
        pair <- lvls[ctr$comps[, q] + 1]
        df_q <- max(edf[labels %in% paste0("s(midsC):covariate", pair)])
        idx <- seq_len(nbasis) + (q - 1) * nbasis
        expected <- waldTest(fit, ctr$L[, idx, drop = FALSE], df = df_q)
        expect_equal(
            out$pvals["E9PAV3", paste0("Field", pair[1], "_VS_Field", pair[2])],
            unname(expected)
        )
    }
})

test_that("non-testable contrasts give NA p values instead of an error", {
    ## A huge tolerance treats every eigenvalue as zero (rank 0)
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    out <- shapeTest(sce = sce, covariate = "SampleType", tol = 1e10)
    expect_identical(dim(out$pvals), c(10L, 1L))
    expect_true(all(is.na(out$pvals)))

    sce$Field <- as.factor(sce$Field)
    out <- shapeTest(sce = sce, covariate = "Field", tol = 1e10)
    expect_identical(dim(out$pvals), c(10L, 7L))
    expect_true(all(is.na(out$pvals)))
})

test_that("waldTest returns a single named p value", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    fit <- shapeTest(sce = sce, covariate = "SampleType")$models[["E9PAV3"]]
    L <- diag(length(coef(fit)))[, 3:4, drop = FALSE]
    expect_named(waldTest(fit, L, df = 2), "pval")
    expect_named(waldTest(fit, L * 0, df = 2), "pval")
    expect_true(is.na(waldTest(fit, L * 0, df = 2)))
})

test_that("waldTest handles contrasts of rank 1", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    fit <- shapeTest(sce = sce, covariate = "SampleType")$models[["E9PAV3"]]
    j <- 3
    b <- coef(fit)[j]
    v <- fit$Vp[j, j]
    expected <- pchisq(b^2 / v, df = 1, lower.tail = FALSE)
    ## Scale the contrast so that 1 / sqrt(variance) is 5, 1.5 and 0.5
    for (target in c(5, 1.5, 0.5)) {
        L <- matrix(0, nrow = length(coef(fit)), ncol = 1)
        L[j, 1] <- 1 / (target * sqrt(v))
        expect_equal(unname(waldTest(fit, L, df = 1)), unname(expected))
    }
})
