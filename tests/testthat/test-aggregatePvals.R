test_that("compute a matrix", {
  leduc <- get_data()
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  expect_true(is(out, "matrix"))
}
)
test_that("compute the correct p values", {
  leduc <- get_data()
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = ~Set)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  expect_equal(out, structure(c(NA, NA, NA, NA, NA, NA, 1.57881004892806e-07, NA,
                                0.156430438707597, 1.7235777530268e-06), dim = c(10L, 1L), dimnames = list(
                                  NULL, "SampleTypeMelanoma_VS_SampleTypeMonocyte")))
}
)

test_that("compute the correct p values", {
  leduc <- get_data()
  leduc$Field <- as.factor(leduc$Field)
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  expect_equal(out, structure(c(NA, NA, NA, NA, NA, NA, 8.35753033439263e-09, 0.448606215296766,
                                0.0181856567593166, 0.263542908819706, NA, NA, NA, NA, NA, NA,
                                7.26889206854185e-10, 0.133874005343969, 0.793150308246409, 0.842487172361624,
                                NA, NA, NA, NA, NA, NA, 3.45305336504389e-05, 0.374491529710233,
                                0.00431153830407926, 0.0943365512455869, NA, NA, NA, NA, NA,
                                NA, 0.00131812842012296, 0.260031496990843, 0.192852792424489,
                                0.79720525514013, NA, NA, NA, NA, NA, NA, 0.0663654356091053,
                                0.533115354298812, 0.0130906858527762, 0.0215923382382018, NA,
                                NA, NA, NA, NA, NA, 0.0124066423544074, 0.517552865090733, 0.220159800679058,
                                0.7467891098581, NA, NA, NA, NA, NA, NA, 0.39990456933177, 0.591604118425169,
                                0.0191954026430765, 0.322226530052631), dim = c(10L, 7L), dimnames = list(
                                  c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2",
                                    "E9PAV3", "E9PRG8", "O00154", "O00159"), c("omnibus", "Field1_VS_Field2",
                                                                               "Field1_VS_Field3", "Field1_VS_Field4", "Field2_VS_Field3",
                                                                               "Field2_VS_Field4", "Field3_VS_Field4"))))
}
)
