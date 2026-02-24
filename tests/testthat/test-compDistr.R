test_that("return a singleCellExperiment", {
  leduc <- get_data()
  out <- compDistr(object = leduc, i = "peptides_log", covariate = "SampleType", nuisance = ~ -1 + Set)
  expect_true(is(out, "QFeatures"))
}
)

test_that("return a singleCellExperiment", {
  leduc <- get_data()
  out <- compDistr(object = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~ -1 + Set)
  expect_true(is(out, "SingleCellExperiment"))
}
)

test_that("compute the correct comparisons", {
  leduc <- get_data()
  leduc$Field <- as.factor(leduc$Field)
  out <- compDistr(object = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = ~ -1 + (1|Set))
  expect_true(all(c("omnibus", "Field1_VS_Field2", "Field1_VS_Field3",
                    "Field1_VS_Field4", "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"
  ) %in% colnames(rowData(out))))
}
)
test_that("compute the correct pvals", {
  leduc <- get_data()
  leduc$Field <- as.factor(leduc$Field)
  out <- compDistr(object = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = ~ -1 + (1|Set), minObs = 20)
  expect_equal(rowData(out)[,c("Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
                               "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4")],
               new("DFrame", rownames = c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6",
                                          "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154", "O00159"),
                   nrows = 10L, elementType = "ANY", elementMetadata = NULL,
                   metadata = list(), listData = list(Field1_VS_Field2 = c(NA,
                                                                           NA, NA, 1, 0.203274571958196, NA, 2.65049431870263e-09, NA,
                                                                           NA, NA), Field1_VS_Field3 = c(NA, NA, NA, 0.0683449020923757,
                                                                                                         0.292534736624142, NA, 9.55707495296236e-05, NA, NA, NA),
                                                      Field1_VS_Field4 = c(NA, NA, NA, 1, 1, NA, 0.0257698986956148,
                                                                           NA, NA, NA), Field2_VS_Field3 = c(NA, NA, NA, 0.0121476840456722,
                                                                                                             1, NA, 1, NA, NA, NA), Field2_VS_Field4 = c(NA, NA, NA,
                                                                                                                                                         1, 0.320509846707837, NA, 0.355553104763042, NA, NA,
                                                                                                                                                         NA), Field3_VS_Field4 = c(NA, NA, NA, 0.201538142625581,
                                                                                                                                                                                   0.292534736624142, NA, 1, NA, NA, NA))))
}
)

test_that("don't allow covariate to be in nuisance formula", {
  leduc <- get_data()
  expect_error(compDistr(object = leduc, i = "peptides_log", covariate = "SampleType", nuisance = ~ -1 + SampleType:Set))
}
)
