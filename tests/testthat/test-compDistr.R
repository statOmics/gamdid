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
  out <- compDistr(object = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = ~ -1 + (1|Set))
  expect_equal(rowData(out)[,c("Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
                               "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4")],
               new("DFrame", rownames = c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6",
                                          "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154", "O00159"),
                   nrows = 10L, elementType = "ANY", elementMetadata = NULL,
                   metadata = list(), listData = list(Field1_VS_Field2 = c(NA,
                                                                           0.457083203173877, 0.047789282551562, 1, 0.152455928968647,
                                                                           NA, 1.98787073902697e-09, NA, NA, NA), Field1_VS_Field3 = c(NA,
                                                                                                                                       0.497584737355003, 0.163793066536348, 0.0512586765692817,
                                                                                                                                       0.219401052468107, NA, 7.16780621472177e-05, NA, NA, NA),
                                                      Field1_VS_Field4 = c(NA, 0.457083203173877, 0.163793066536348,
                                                                           1, 1, NA, 0.0193274240217111, NA, NA, NA), Field2_VS_Field3 = c(NA,
                                                                                                                                           0.248368490743068, 1, 0.00911076303425413, 1, NA, 1,
                                                                                                                                           NA, NA, NA), Field2_VS_Field4 = c(NA, 0.221170928708685,
                                                                                                                                                                             0.708927612796853, 1, 0.240382385030878, NA, 0.266664828572281,
                                                                                                                                                                             NA, NA, NA), Field3_VS_Field4 = c(NA, 0.457083203173877,
                                                                                                                                                                                                               0.708927612796853, 0.151153606969186, 0.219401052468107,
                                                                                                                                                                                                               NA, 1, NA, NA, NA))))
}
)

test_that("don't allow covariate to be in nuisance formula", {
  leduc <- get_data()
  expect_error(compDistr(object = leduc, i = "peptides_log", covariate = "SampleType", nuisance = ~ -1 + SampleType:Set))
}
)
