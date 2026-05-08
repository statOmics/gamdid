test_that("returns a list", {
  out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  expect_true(is(out, "list"))
}
)

test_that("list with GAM models and pvals", {
  out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  expect_identical(names(out), c("models", "pvals"))
}
)

test_that("pvals object is a matrix", {
  out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  expect_true(is(out$pvals, "matrix"))
}
)

test_that("models contain StatModel objects", {
  out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  expect_true(is(out$models[[1]], "gam"))
}
)

test_that("compute correct p values", {
  out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  expect_equal(out$pvals, structure(c(0.0183805817026951, 0.000386393734274719, 0.5383450730135,
                                      0.00546027282729357, 0.595862779628328, 0.0150442734089543, 4.80544360393079e-05,
                                      0.0199681251740768, 0.335987280739658, 0.00189898572291178),
                                    dim = c(10L, 1L), dimnames = list(NULL, "SampleTypeMelanoma_VS_SampleTypeMonocyte")),tolerance = 1e-6)
}
)

test_that("correct dimensions", {
  leduc_testdata$Field <- as.factor(leduc_testdata$Field)
  out <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
  expect_equal(dim(out$pvals), c(10, ncol(combn(nlevels(leduc_testdata$Field),2))+1))
}
)

