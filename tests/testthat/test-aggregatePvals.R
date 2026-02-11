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
  expect_equal(out, structure(c(0.410594648663631, 0.000110430084729993, 6.39407075613526e-05,
                                0.00280535295757099, 0.00056747709848988, 0.00316630933034802,
                                8.35753033439263e-09, 0.448606215296766, 0.0181856567593166,
                                0.263542908819706, 0.524313192786383, 0.0254946450328052, 5.4615176383484e-05,
                                0.579721583895562, 0.00658743935733475, 0.0305427737976129, 7.26889206854185e-10,
                                0.133874005343969, 0.793150308246409, 0.842487172361624, 0.304517567947597,
                                0.150267375373211, 1.7361236982166e-05, 0.0031954979344023, 0.0134708090376828,
                                0.356647829570729, 3.45305336504389e-05, 0.374491529710233, 0.00431153830407926,
                                0.0943365512455869, 0.560490191880402, 0.00893256294269912, 0.00233023972782321,
                                0.231046473333615, 0.368929237971937, 0.234456681373251, 0.00131812842012296,
                                0.260031496990843, 0.192852792424489, 0.79720525514013, 0.597843675865837,
                                0.000552821607559222, 0.411613841071528, 0.000828996468356803,
                                0.486712677329631, 0.000787634319692606, 0.0663654356091053,
                                0.533115354298812, 0.0130906858527762, 0.0215923382382018, 0.296155962330694,
                                1.06601923577627e-05, 0.0720300397434854, 0.2906224486429, 0.00105349713391184,
                                0.587356466535041, 0.0124066423544074, 0.517552865090733, 0.220159800679058,
                                0.7467891098581, 0.107449917394594, 0.0461034981302682, 0.021556121883619,
                                0.016266442849401, 0.00206565992262635, 0.0239589954834463, 0.39990456933177,
                                0.591604118425169, 0.0191954026430765, 0.322226530052631), dim = c(10L,
                                                                                                   7L), dimnames = list(c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6",
                                                                                                                          "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154", "O00159"),
                                                                                                                        c("omnibus", "Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
                                                                                                                          "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"
                                                                                                                        ))))
}
)
