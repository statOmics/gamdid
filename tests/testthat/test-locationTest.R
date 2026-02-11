test_that("returns a list", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  expect_true(is(out, "list"))
  }
)

test_that("list with pvals", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  expect_identical(names(out), "pvals")
  }
)

test_that("pvals object is a matrix", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  expect_true(is(out$pvals, "matrix"))
  }
)

test_that("list with pvals matrix, msqrob_models en msqrob_formula", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ (1|Set))
  expect_identical(names(out), c("pvals", "msqrob_models", "msqrob_formula"))
  }
)

test_that("list with pvals matrix, msqrob_models en msqrob_formula (2)", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ -1 + Set)
  expect_identical(names(out), c("pvals", "msqrob_models", "msqrob_formula"))
  }
)

test_that("pvals object is a matrix (for nuisance correction)", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ -1 + (1|Set))
  expect_true(is(out$pvals, "matrix"))
  }
)

test_that("pvals object is a matrix (for nuisance correction 2)", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ Set)
  expect_true(is(out$pvals, "matrix"))
  }
)

test_that("msqrobModels contain StatModel objects", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ Set)
  expect_true(is(out$msqrob_models[[1]], "StatModel"))
  }
)

test_that("msqrob_formula contain formula objects", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ Set)
  expect_true(is(out$msqrob_formula, "formula"))
  }
)

test_that("compute correct p values", {
  leduc <- get_data()
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ (1|Set))
  expect_equal(out$pvals, structure(c(0.237037266770739, 0.981984915759283, 0.996259580193254,
                                        0.000181307193181744, 0.344405576858417, 0.000708840124653451,
                                        3.35828989214574e-08, 0.226269016259427, 0.235860756436294, 4.9922259329987e-05
  ), dim = c(10L, 1L), dimnames = list(NULL, "SampleTypeMelanoma_VS_SampleTypeMonocyte")))
  }
)

test_that("correct dimensions", {
  leduc <- get_data()
  leduc$Field <- as.factor(leduc$Field)
  out <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = ~ (1|Set))
  expect_equal(dim(out$pvals), c(10, ncol(combn(nlevels(leduc$Field),2))+1))
  }
)

