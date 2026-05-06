test_that("return a singleCellExperiment", {
  out <- compDistr(object = leduc_testdata, i = "peptides_log", covariate = "SampleType", nuisance = ~ -1 + Set)
  expect_true(is(out, "QFeatures"))
}
)

test_that("return error", {
  expect_error(compDistr(object = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = ~ -1 + Set))
}
)

test_that("compute the correct comparisons", {
  leduc_testdata$Field <- as.factor(leduc_testdata$Field)
  out <- compDistr(object = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = ~ -1 + (1|Set))
  expect_true(all(c("omnibus", "Field1_VS_Field2", "Field1_VS_Field3",
                    "Field1_VS_Field4", "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"
  ) %in% colnames(rowData(out))))
}
)
test_that("compute the correct pvals", {
  leduc_testdata$Field <- as.factor(leduc_testdata$Field)
  out <- compDistr(object = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = ~ -1 + (1|Set), minObs = 20)
  expect_equal(rowData(out)[,c("Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
                               "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4")],
               new("DFrame", rownames = c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6",
                                          "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154", "O00159"),
                   nrows = 10L, elementType = "ANY", elementMetadata = NULL,
                   metadata = list(), listData = list(Field1_VS_Field2 = c(NA,
                                                                           NA, NA, 1, 0.203274571842, NA, 2.65049430261989e-09, NA,
                                                                           NA, NA), Field1_VS_Field3 = c(NA, NA, NA, 0.0683449020924028,
                                                                                                         0.292534737658064, NA, 9.55707522003746e-05, NA, NA, NA),
                                                      Field1_VS_Field4 = c(NA, NA, NA, 1, 1, NA, 0.0257698993492636,
                                                                           NA, NA, NA), Field2_VS_Field3 = c(NA, NA, NA, 0.0121476840456777,
                                                                                                             1, NA, 1, NA, NA, NA), Field2_VS_Field4 = c(NA, NA, NA,
                                                                                                                                                         1, 0.320509846016792, NA, 0.355553099693652, NA, NA,
                                                                                                                                                         NA), Field3_VS_Field4 = c(NA, NA, NA, 0.201538142625652,
                                                                                                                                                                                   0.292534737658064, NA, 1, NA, NA, NA))))
}
)

test_that("don't allow covariate to be in nuisance formula", {
  expect_error(compDistr(object = leduc_testdata, i = "peptides_log", covariate = "SampleType", nuisance = ~ -1 + SampleType:Set))
}
)
