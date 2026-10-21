## mixedDesign() must give the same design matrices and random effect names
## as fitting the model with lme4::lmer() on a dummy response
dummyLmerDesign <- function(formula, colDat) {
    fit <- lme4::lmer(paste0("y~", as.character(formula)[2]),
        data = data.frame(colDat, y = 0)
    )
    re <- lme4::ranef(fit)
    list(
        X = lme4::getME(fit, "X"), Z = lme4::getME(fit, "Z"),
        colnames_Z = paste0("(Intercept)", unlist(lapply(names(re), function(nm) {
            paste0(nm, rownames(re[[nm]]))
        })))
    )
}

test_that("mixedDesign equals the design of a fitted lmer model", {
    cd <- colData(getWithColData(leduc_testdata, "proteins_norm"))
    cd$Field <- as.factor(cd$Field)
    formulas <- list(
        ~ SampleType + (1 | Set), ~ SampleType + -1 + (1 | Set),
        ~ SampleType + lcbatch + (1 | Set), ~ Field + (1 | Set),
        ~ SampleType + (1 | Set) + (1 | lcbatch),
        ~ SampleType + (1 | lcbatch / Set)
    )
    for (f in formulas) {
        expect_identical(mixedDesign(f, cd), dummyLmerDesign(f, cd))
    }
})

test_that("only random intercepts are allowed in the nuisance formula", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    for (nuisance in list(
        ~ (Diameter | Set), ~ (1 + Diameter | Set),
        ~ (0 + Diameter | Set), ~ (1 || Set),
        ~ (Diameter || Set)
    )) {
        expect_error(
            compDistr(sce, covariate = "SampleType", nuisance = nuisance),
            "only random intercepts"
        )
    }
    expect_no_error(compDistr(sce,
        covariate = "SampleType",
        nuisance = ~ (1 | Set) + (1 | lcbatch)
    ))
    expect_no_error(compDistr(sce,
        covariate = "SampleType",
        nuisance = ~ lcbatch + (1 | Set)
    ))
})
