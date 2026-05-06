test_that("returns a SingleCellExperiment", {
  inp <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"),
                      covariate = "SampleType", nuisance = ~-1+lcbatch)
  out <- partialResiduals(sce = getWithColData(leduc_testdata, "proteins_norm"),
                          covariate = "SampleType",
                          msqrob_models = inp$msqrob_models,
                          msqrob_formula = inp$msqrob_formula)
  expect_true(is(out, "SingleCellExperiment"))
}
)

test_that("compute the correct partial residuals", {
  inp <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"),
                      covariate = "SampleType", nuisance = ~-1+lcbatch)
  out <- partialResiduals(sce = getWithColData(leduc_testdata, "proteins_norm"),
                          covariate = "SampleType",
                          msqrob_models = inp$msqrob_models,
                          msqrob_formula = inp$msqrob_formula)
  expect_equal(assay(out)[1:10,1:10], structure(c(NA, NA, NA, NA, NA, NA, -1.02744960212803, NA, -0.577829257901651,
                                                  NA, NA, NA, NA, NA, NA, NA, -0.192588429358822, NA, 0.496489500072798,
                                                  NA, NA, NA, NA, NA, NA, NA, 0.212044067589839, NA, 0.230113948960856,
                                                  0.162283310969211, NA, NA, NA, NA, NA, NA, -0.530926538143929,
                                                  NA, -0.280171670944339, NA, NA, NA, NA, NA, NA, NA, -1.16414339987,
                                                  NA, -0.0145558889736888, NA, NA, NA, NA, NA, NA, NA, -0.445782235122565,
                                                  NA, -0.925957675074437, NA, NA, NA, NA, NA, NA, NA, -0.26631850777565,
                                                  NA, -0.936704604043534, NA, NA, NA, NA, NA, NA, NA, -1.20229498423774,
                                                  NA, -1.10033860586228, -1.32688737687657, NA, NA, NA, NA, NA,
                                                  NA, -0.234381798468667, NA, -0.58737107593997, -0.51924004054175,
                                                  NA, NA, NA, NA, NA, NA, -1.11889998757941, NA, -1.43347756004521,
                                                  NA), dim = c(10L, 10L), dimnames = list(c("A0AVT1", "A2A368",
                                                                                            "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154",
                                                                                            "O00159"), c("eAL00219RI5", "eAL00219RI6", "eAL00219RI7", "eAL00219RI8",
                                                                                                         "eAL00219RI9", "eAL00219RI10", "eAL00219RI11", "eAL00219RI13",
                                                                                                         "eAL00219RI14", "eAL00219RI15"))))
}
)

test_that("compute the correct partial residuals (2)", {
  inp <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"),
                      covariate = "SampleType", nuisance = ~(1|Set))
  out <- partialResiduals(sce = getWithColData(leduc_testdata, "proteins_norm"),
                          covariate = "SampleType",
                          msqrob_models = inp$msqrob_models,
                          msqrob_formula = inp$msqrob_formula)
  expect_equal(assay(out)[1:10,1:10], structure(c(NA, NA, NA, NA, NA, NA, -0.32428438976326, NA, -0.389761602365137,
                                                  NA, NA, NA, NA, NA, NA, NA, 0.510576783005946, NA, 0.684557155609313,
                                                  NA, NA, NA, NA, NA, NA, NA, 0.915209279954607, NA, 0.418181604497371,
                                                  0.21827795810478, NA, NA, NA, NA, NA, NA, 0.172238674220838,
                                                  NA, -0.0921040154078245, NA, NA, NA, NA, NA, NA, NA, -0.460978187505234,
                                                  NA, 0.173511766562826, NA, NA, NA, NA, NA, NA, NA, 0.257382977242202,
                                                  NA, -0.737890019537923, NA, NA, NA, NA, NA, NA, NA, 0.436846704589118,
                                                  NA, -0.748636948507019, NA, NA, NA, NA, NA, NA, NA, -0.499129771872977,
                                                  NA, -0.912270950325765, -1.270892729741, NA, NA, NA, NA, NA,
                                                  NA, 0.468783413896101, NA, -0.399303420403455, -0.463245393406181,
                                                  NA, NA, NA, NA, NA, NA, -0.415734775214642, NA, -1.2454099045087,
                                                  NA), dim = c(10L, 10L), dimnames = list(c("A0AVT1", "A2A368",
                                                                                            "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154",
                                                                                            "O00159"), c("eAL00219RI5", "eAL00219RI6", "eAL00219RI7", "eAL00219RI8",
                                                                                                         "eAL00219RI9", "eAL00219RI10", "eAL00219RI11", "eAL00219RI13",
                                                                                                         "eAL00219RI14", "eAL00219RI15"))))
}
)

test_that("assay dimensions are preserved", {
  inp <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"),
                      covariate = "SampleType", nuisance = ~(1|Set))
  out <- partialResiduals(sce = getWithColData(leduc_testdata, "proteins_norm")[1:10,],
                          covariate = "SampleType",
                          msqrob_models = inp$msqrob_models,
                          msqrob_formula = inp$msqrob_formula)
  expect_identical(dim(assay(out)), dim(assay(getWithColData(leduc_testdata, "proteins_norm"))))
}
)

test_that("Handles NA values", {
  leduc_error <- getWithColData(leduc_testdata, "proteins_norm")
  assay(leduc_error)[1,] <- NA
  inp_error <- locationTest(sce = leduc_error,
                      covariate = "SampleType", nuisance = ~(1|Set))
  out <- partialResiduals(sce = leduc_error,
                          covariate = "SampleType",
                          msqrob_models = inp_error$msqrob_models,
                          msqrob_formula = inp_error$msqrob_formula)
  expect_true(all(is.na(assay(out)[1,])))
}
)





