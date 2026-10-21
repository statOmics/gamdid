## The adjust options of compDistr, checked against the unadjusted p values
adjustRuns <- function(covariate) {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    sce$Field <- as.factor(sce$Field)
    lapply(
        c(stagewise = "stagewise", contrastwise = "contrastwise", none = "none"),
        function(a) {
            rowData(compDistr(sce,
                covariate = covariate, adjust = a,
                minObs = 20
            ))
        }
    )
}

test_that("adjust = 'none' returns three p values per contrast", {
    rd <- adjustRuns("Field")$none
    contrasts <- c("omnibus", combn(paste0("Field", 1:4), 2, paste, collapse = "_VS_"))
    expected <- paste0(rep(contrasts, 3), "_", rep(c("GAM", "msqrob2", "aggregated"), each = 7))
    expect_setequal(intersect(expected, colnames(rd)), expected)
})

test_that("adjust = 'contrastwise' is BH on each aggregated column", {
    for (covariate in c("SampleType", "Field")) {
        r <- adjustRuns(covariate)
        agg <- grep("_aggregated$", colnames(r$none), value = TRUE)
        for (a in agg) {
            expect_equal(
                r$contrastwise[[sub("_aggregated$", "", a)]],
                p.adjust(r$none[[a]], method = "BH")
            )
        }
    }
})

test_that("adjust = 'stagewise' equals contrastwise for two groups", {
    r <- adjustRuns("SampleType")
    col <- "SampleTypeMelanoma_VS_SampleTypeMonocyte"
    expect_identical(r$stagewise[[col]], r$contrastwise[[col]])
})

test_that("adjust = 'stagewise' equals stageR for more than two groups", {
    r <- adjustRuns("Field")
    p <- as.matrix(as.data.frame(r$none[grep("_aggregated$", colnames(r$none))]))
    colnames(p) <- sub("_aggregated$", "", colnames(p))
    ok <- !is.na(p[, "omnibus"])
    sr <- stageR::stageR(
        pScreen = p[ok, "omnibus"],
        pConfirmation = p[ok, -1, drop = FALSE]
    )
    ref <- stageR::getAdjustedPValues(
        stageR::stageWiseAdjustment(sr,
            method = "holm", alpha = 0.05,
            allowNA = TRUE
        ),
        onlySignificantGenes = FALSE, order = FALSE
    )
    got <- as.matrix(as.data.frame(r$stagewise[colnames(p)]))
    expect_equal(unname(got[ok, ]), unname(ref))
    expect_true(all(is.na(got[!ok, ])))
})
