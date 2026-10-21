test_that("topFeatures sorts by the comparison for two groups", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    out <- compDistr(sce, covariate = "SampleType")
    col <- "SampleTypeMelanoma_VS_SampleTypeMonocyte"
    top <- topFeatures(out)
    expect_s3_class(top, "data.frame")
    expect_identical(nrow(top), 10L)
    expect_false(is.unsorted(top[[col]], na.rm = TRUE))
    ## NA p values come last
    expect_false(is.unsorted(is.na(top[[col]])))
    expect_identical(rownames(topFeatures(out, n = 3)), rownames(top)[1:3])
    expect_identical(topFeatures(out, sort.by = col), top)
})

test_that("topFeatures sorts by the omnibus test for more than two groups", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    out <- compDistr(sce, covariate = "Field", minObs = 20)
    top <- topFeatures(out)
    expect_identical(top, topFeatures(out, sort.by = "omnibus"))
    expect_false(is.unsorted(top$omnibus, na.rm = TRUE))
    ## another column can be chosen
    pw <- topFeatures(out, sort.by = "Field1_VS_Field2")
    expect_false(is.unsorted(pw$Field1_VS_Field2, na.rm = TRUE))
})

test_that("topFeatures errors", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    out <- compDistr(sce, covariate = "SampleType")
    expect_error(topFeatures(out, sort.by = "nonexisting"), "not found in rowData")
    none <- compDistr(sce, covariate = "SampleType", adjust = "none")
    expect_error(topFeatures(none), "specify 'sort.by'")
    expect_no_error(topFeatures(
        none,
        sort.by = "SampleTypeMelanoma_VS_SampleTypeMonocyte_aggregated"
    ))
})

test_that("topFeatures works on QFeatures objects", {
    ## getWithColData() warns about colData columns that the test data stores
    ## twice; this is a property of the test data, not of gamdid
    suppressWarnings({
        out <- compDistr(leduc_testdata,
            i = "proteins_norm", covariate = "SampleType"
        )
        expect_identical(
            topFeatures(out, i = "proteins_norm"),
            topFeatures(getWithColData(out, "proteins_norm"))
        )
    })
})
