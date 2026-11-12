## Results must not depend on how the covariate and nuisance variables are
## named, e.g. a nuisance name that contains the covariate name

renameRun <- function(nuisance, k, adjust = "none", covariate = "Group") {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    field <- as.character(sce$Field)
    sce[[covariate]] <- if (k == 2) sce$SampleType else factor(ifelse(field == "4", "3", field))
    sce$Block <- sce$Set
    sce$SetGroup <- sce$Set
    sce$CellXtype <- sce$Set
    out <- compDistr(sce,
        covariate = covariate, nuisance = nuisance,
        adjust = adjust, minObs = 20
    )
    rd <- rowData(out)
    as.matrix(as.data.frame(rd[grep("_GAM$|_msqrob2$", colnames(rd))]))
}

test_that("a nuisance name containing the covariate name gives the same result", {
    expect_equal(renameRun(~SetGroup, 2), renameRun(~Block, 2))
    expect_equal(renameRun(~ (1 | SetGroup), 2), renameRun(~ (1 | Block), 2))
    expect_equal(renameRun(~SetGroup, 3), renameRun(~Block, 3))
    expect_equal(renameRun(~ (1 | SetGroup), 3), renameRun(~ (1 | Block), 3))
})

test_that("covariate names are not treated as regular expressions", {
    ## "." in "Cell.type" matched any character, e.g. the X in "CellXtype"
    expect_equal(
        renameRun(~CellXtype, 2, covariate = "Cell.type"),
        renameRun(~Block, 2, covariate = "Cell.type")
    )
})

test_that("level names containing the covariate name are handled", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Cell <- ifelse(sce$SampleType == "Melanoma", "CellA", "CellB")
    out <- compDistr(sce, covariate = "Cell")
    expect_true("CellCellA_VS_CellCellB" %in% colnames(rowData(out)))
    p <- visFit(out, feature = "E9PAV3")
    groups <- unique(as.character(p$layers[[1]]$data$Cell))
    expect_setequal(groups, c("CellA", "CellB"))
})

test_that("only the covariate itself is forbidden in the nuisance formula", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Sample <- sce$SampleType
    sce$Sample.batch <- sce$Set
    expect_no_error(compDistr(sce, covariate = "Sample", nuisance = ~Sample.batch))
    expect_error(
        compDistr(sce, covariate = "Sample", nuisance = ~ Set + Sample:Set),
        "covariate of interest cannot be in nuisance formula"
    )
})

test_that("aggregatePvals refuses tests with different contrasts", {
    shape <- matrix(0.1, nrow = 2, ncol = 1, dimnames = list(NULL, "A_VS_B"))
    location <- matrix(0.1, nrow = 2, ncol = 1, dimnames = list(NULL, "B_VS_A"))
    expect_error(aggregatePvals(shape, location), "different contrasts")
})
