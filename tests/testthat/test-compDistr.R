test_that("returns a QFeatures object", {
    out <- compDistr(object = leduc_testdata, i = "peptides_log", covariate = "SampleType", nuisance = ~ -1 + Set)
    expect_true(is(out, "QFeatures"))
})

## TODO: re-enable once msqrob2 >= 1.21.1 (estimability handling of
## non-estimable nuisance levels) is available on Bioconductor devel.
## Fails with msqrob2 1.21.0, which returns empty fits for these models.
# test_that("fixed effect nuisance without intercept", {
#     out <- compDistr(object = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~ -1 + Set, adjust = "none")
#     expect_false(anyNA(rowData(out)$SampleTypeMelanoma_VS_SampleTypeMonocyte_GAM))
# })

test_that("compute the correct comparisons", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out <- compDistr(object = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = ~ -1 + (1 | Set))
    expect_true(all(c(
        "omnibus", "Field1_VS_Field2", "Field1_VS_Field3",
        "Field1_VS_Field4", "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"
    ) %in% colnames(rowData(out))))
})
test_that("compute the correct pvals", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out <- compDistr(object = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = ~ -1 + (1 | Set), minObs = 20)
    expect_equal(
        rowData(out)[, c(
            "Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
            "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"
        )],
        new("DFrame",
            rownames = c(
                "A0AVT1", "A2A368", "A6H8Y1", "A6NCS6",
                "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154", "O00159"
            ),
            nrows = 10L, elementType = "ANY", elementMetadata = NULL,
            metadata = list(), listData = list(Field1_VS_Field2 = c(
                NA,
                0.51623135963743494, NA, 1, 0.17151291999168761, NA, 2.2363545678252193e-09,
                NA, NA, NA
            ), Field1_VS_Field3 = c(
                NA, 1, NA, 0.061301669528861374,
                0.30317146061423539, NA, 0.0001979436755479931, NA, NA, NA
            ), Field1_VS_Field4 = c(
                NA, 0.51623135963743494, NA, 1, 1,
                NA, 0.042301390102330187, NA, NA, NA
            ), Field2_VS_Field3 = c(
                NA,
                0.41751933826266591, NA, 0.011417411242869107, 1, NA, 1,
                NA, NA, NA
            ), Field2_VS_Field4 = c(
                NA, 0.24930433452072756,
                NA, 1, 0.27320189759487123, NA, 0.17838620784364811, NA,
                NA, NA
            ), Field3_VS_Field4 = c(
                NA, 0.51623135963743494, NA,
                0.17004780784039408, 0.27184830110984626, NA, 1, NA, NA,
                NA
            ))
        ),
        tolerance = 1e-6
    )
})

test_that("don't allow covariate to be in nuisance formula", {
    expect_error(compDistr(object = leduc_testdata, i = "peptides_log", covariate = "SampleType", nuisance = ~ -1 + SampleType:Set))
})

## TODO: re-enable once msqrob2 >= 1.21.1 (estimability handling of
## non-estimable nuisance levels) is available on Bioconductor devel.
## Fails with msqrob2 1.21.0, which returns empty fits for these models.
# test_that("fixed effect nuisance does not remove the GAM test", {
#     sce <- getWithColData(leduc_testdata, "proteins_norm")
#     out <- compDistr(
#         object = sce, covariate = "SampleType", nuisance = ~Set,
#         adjust = "none", minObs = 20
#     )
#     ref <- compDistr(
#         object = sce, covariate = "SampleType", adjust = "none",
#         minObs = 20
#     )
#     col <- "SampleTypeMelanoma_VS_SampleTypeMonocyte_GAM"
#     expect_true(any(!is.na(rowData(out)[[col]])))
#     expect_identical(is.na(rowData(out)[[col]]), is.na(rowData(ref)[[col]]))
# })
