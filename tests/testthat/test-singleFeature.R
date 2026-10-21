test_that("compDistr works for a single feature", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    one <- sce["E9PAV3", ]
    for (adjust in c("stagewise", "contrastwise", "none")) {
        out <- compDistr(one, covariate = "SampleType", adjust = adjust)
        expect_identical(nrow(rowData(out)), 1L)
        out <- compDistr(one, covariate = "Field", adjust = adjust)
        expect_identical(nrow(rowData(out)), 1L)
        expect_true("Field1_VS_Field2" %in% sub(
            "_(GAM|msqrob2|aggregated)$", "",
            colnames(rowData(out))
        ))
    }
})

test_that("single feature GAM p values equal those of the full analysis", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    full <- compDistr(sce, covariate = "Field", adjust = "none")
    one <- compDistr(sce["E9PAV3", ], covariate = "Field", adjust = "none")
    gam <- grep("_GAM$", colnames(rowData(one)), value = TRUE)
    expect_equal(
        unlist(as.list(rowData(one)[gam])),
        unlist(as.list(rowData(full)["E9PAV3", gam]))
    )
})

test_that("adjustOutput keeps a matrix for a single feature", {
    pval <- matrix(c(0.01, 0.02, 0.03),
        nrow = 1,
        dimnames = list("a", c("omnibus", "A_VS_B", "A_VS_C"))
    )
    for (adjust in c("stagewise", "contrastwise", "none")) {
        out <- adjustOutput(
            pval = pval, pval_shape = pval, pval_location = pval,
            adjust = adjust, alpha = 0.05
        )
        expect_true(is.matrix(out$pvals))
        expect_identical(nrow(out$pvals), 1L)
    }
    ## BH with a single test leaves the p values unchanged
    out <- adjustOutput(
        pval = pval, pval_shape = pval, pval_location = pval,
        adjust = "contrastwise", alpha = 0.05
    )
    expect_equal(out$pvals, pval)
})
