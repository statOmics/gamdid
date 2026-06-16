test_that("compute a matrix", {
  out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, alpha = 0.05)
  expect_true(is(out_adj$pvals, "matrix"))
}
)

test_that("compute the correct comparisons", {
  out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), "SampleTypeMelanoma_VS_SampleTypeMonocyte")
}
)

test_that("compute the correct comparisons if not adjusted", {
  out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "none", alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), c("SampleTypeMelanoma_VS_SampleTypeMonocyte_GAM", "SampleTypeMelanoma_VS_SampleTypeMonocyte_msqrob2",
                                              "SampleTypeMelanoma_VS_SampleTypeMonocyte_aggregated"))
}
)

test_that("compute the correct comparisons for multiple levels", {
  leduc_testdata$Field <- as.factor(leduc_testdata$Field)
  out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), c("omnibus", "Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
                                              "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"))
}
)

test_that("compute the correct comparisons for multiple levels if not adjusted", {
  leduc_testdata$Field <- as.factor(leduc_testdata$Field)
  out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "none", alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), c("omnibus_GAM", "Field1_VS_Field2_GAM", "Field1_VS_Field3_GAM",
                                              "Field1_VS_Field4_GAM", "Field2_VS_Field3_GAM", "Field2_VS_Field4_GAM",
                                              "Field3_VS_Field4_GAM", "omnibus_msqrob2", "Field1_VS_Field2_msqrob2",
                                              "Field1_VS_Field3_msqrob2", "Field1_VS_Field4_msqrob2", "Field2_VS_Field3_msqrob2",
                                              "Field2_VS_Field4_msqrob2", "Field3_VS_Field4_msqrob2", "omnibus_aggregated",
                                              "Field1_VS_Field2_aggregated", "Field1_VS_Field3_aggregated",
                                              "Field1_VS_Field4_aggregated", "Field2_VS_Field3_aggregated",
                                              "Field2_VS_Field4_aggregated", "Field3_VS_Field4_aggregated"))
}
)

test_that("compute the correct p values", {
  leduc_testdata$Field <- as.factor(leduc_testdata$Field)
  out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "none", alpha = 0.05)
  expect_equal(out_adj$pvals, structure(c(NA, NA, NA, NA, NA, NA, 1.02970302066109e-08, 0.645355338416144,
                                          0.0243131618874192, 0.159018157894301, NA, NA, NA, NA, NA, NA,
                                          7.62813598944185e-07, 0.349379221329827, 0.767411608329258, 0.742855803372467,
                                          NA, NA, NA, NA, NA, NA, 3.80553935626226e-05, 0.51173691635103,
                                          0.00529657905901098, 0.0535803976748089, NA, NA, NA, NA, NA,
                                          NA, 0.00159889993532359, 0.659545979678807, 0.276606130905013,
                                          0.681130986975393, NA, NA, NA, NA, NA, NA, 0.198367211473923,
                                          0.715798574445476, 0.0826809979990406, 0.0110852559978545, NA,
                                          NA, NA, NA, NA, NA, 0.0342849800069655, 0.458152518940389, 0.226413636580256,
                                          0.620595162674051, NA, NA, NA, NA, NA, NA, 0.3422370644001, 0.490889654903771,
                                          0.0100666932518853, 0.256297366966999, 0.298515338235435, 5.52934264946103e-05,
                                          3.22515236672344e-05, 0.00148842081394629, 0.000285900613288176,
                                          0.00160485886393992, 7.03285504563089e-09, 0.343793946134101,
                                          0.014525004218628, 0.769049142229933, 0.630210627806447, 0.014262367584901,
                                          2.74243207509011e-05, 0.51116311522195, 0.00355260702911753,
                                          0.0172333766913221, 3.63617850135757e-10, 0.0828006449702246,
                                          0.820675457866417, 0.972983070606951, 0.30111193452229, 0.223573507666662,
                                          8.68220838272122e-06, 0.00185846539337958, 0.00739385418798584,
                                          0.304513469775974, 3.16032957516972e-05, 0.295294863557948, 0.00363543228815316,
                                          0.394143303106639, 0.4800282189366, 0.00474794117676338, 0.0012905396056319,
                                          0.168048877732769, 0.281227440907842, 0.140275094858713, 0.00112123607830241,
                                          0.161938562078789, 0.14803070067154, 0.960967605868819, 0.491705750029814,
                                          0.000277463807031639, 0.28283722392383, 0.000463046443024833,
                                          0.416487881220139, 0.000411464965843538, 0.0398485539035542,
                                          0.424720025401639, 0.00710804267571234, 0.413987307197512, 0.188896351385011,
                                          5.33129060450405e-06, 0.0514884323402035, 0.28049956177434, 0.000545834531271173,
                                          0.434923612669141, 0.00757365310313541, 0.594650424515641, 0.214242156709305,
                                          0.937404064948075, 0.0641593116281649, 0.0281197569675531, 0.0121158622560556,
                                          0.00993108430389344, 0.00108129734098058, 0.0126325547001453,
                                          0.480944504759914, 0.744312677862637, 0.206010186690516, 0.433821455135794,
                                          NA, NA, NA, NA, NA, NA, 8.35753033439263e-09, 0.448606215296766,
                                          0.0181856567593166, 0.263542908819706, NA, NA, NA, NA, NA, NA,
                                          7.26889206854185e-10, 0.133874005343969, 0.793150308246409, 0.842487172361624,
                                          NA, NA, NA, NA, NA, NA, 3.45305336504389e-05, 0.374491529710233,
                                          0.00431153830407926, 0.0943365512455869, NA, NA, NA, NA, NA,
                                          NA, 0.00131812842012296, 0.260031496990843, 0.192852792424489,
                                          0.79720525514013, NA, NA, NA, NA, NA, NA, 0.0663654356091053,
                                          0.533115354298812, 0.0130906858527762, 0.0215923382382018, NA,
                                          NA, NA, NA, NA, NA, 0.0124066423544074, 0.517552865090733, 0.220159800679058,
                                          0.7467891098581, NA, NA, NA, NA, NA, NA, 0.39990456933177, 0.591604118425169,
                                          0.0191954026430765, 0.322226530052631), dim = c(10L, 21L), dimnames = list(
                                            c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2",
                                              "E9PAV3", "E9PRG8", "O00154", "O00159"), c("omnibus_GAM",
                                                                                         "Field1_VS_Field2_GAM", "Field1_VS_Field3_GAM", "Field1_VS_Field4_GAM",
                                                                                         "Field2_VS_Field3_GAM", "Field2_VS_Field4_GAM", "Field3_VS_Field4_GAM",
                                                                                         "omnibus_msqrob2", "Field1_VS_Field2_msqrob2", "Field1_VS_Field3_msqrob2",
                                                                                         "Field1_VS_Field4_msqrob2", "Field2_VS_Field3_msqrob2",
                                                                                         "Field2_VS_Field4_msqrob2", "Field3_VS_Field4_msqrob2",
                                                                                         "omnibus_aggregated", "Field1_VS_Field2_aggregated", "Field1_VS_Field3_aggregated",
                                                                                         "Field1_VS_Field4_aggregated", "Field2_VS_Field3_aggregated",
                                                                                         "Field2_VS_Field4_aggregated", "Field3_VS_Field4_aggregated"
                                              ))))
}
)

test_that("compute the correct p values contrastwise", {
  leduc_testdata$Field <- as.factor(leduc_testdata$Field)
  out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "contrastwise", alpha = 0.05)
  expect_equal(out_adj$pvals, structure(c(NA, NA, NA, NA, NA, NA, 3.34301213373063e-08, 0.448606215296644,
                                          0.0363713135333237, 0.351390543317299, NA, NA, NA, NA, NA, NA,
                                          2.90755682741918e-09, 0.26774801068781, 0.842487173535532, 0.842487173535532,
                                          NA, NA, NA, NA, NA, NA, 0.000138122134601815, 0.374491529710171,
                                          0.00862307661805698, 0.125782067607922, NA, NA, NA, NA, NA, NA,
                                          0.00527251368049507, 0.346708662654341, 0.346708662654341, 0.797205254641515,
                                          NA, NA, NA, NA, NA, NA, 0.0884872474788192, 0.533115354298752,
                                          0.0431846755749045, 0.0431846755749045, NA, NA, NA, NA, NA, NA,
                                          0.0496265694176787, 0.690070486787598, 0.440319600977197, 0.746789108056538,
                                          NA, NA, NA, NA, NA, NA, 0.533206092442361, 0.591604118425151,
                                          0.0767816104738106, 0.533206092442361), dim = c(10L, 7L), dimnames = list(
                                            c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2",
                                              "E9PAV3", "E9PRG8", "O00154", "O00159"), c("omnibus", "Field1_VS_Field2",
                                                                                         "Field1_VS_Field3", "Field1_VS_Field4", "Field2_VS_Field3",
                                                                                         "Field2_VS_Field4", "Field3_VS_Field4"))))
}
)
