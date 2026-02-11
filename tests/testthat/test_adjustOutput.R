test_that("compute a matrix", {
  leduc <- get_data()
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = T, alpha = 0.05)
  expect_true(is(out_adj$pvals, "matrix"))
}
)

test_that("compute the correct comparisons", {
  leduc <- get_data()
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = T, alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), "SampleTypeMelanoma_VS_SampleTypeMonocyte")
}
)

test_that("compute the correct comparisons if not adjusted", {
  leduc <- get_data()
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "SampleType")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = F, alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), c("SampleTypeMelanoma_VS_SampleTypeMonocyte_shape", "SampleTypeMelanoma_VS_SampleTypeMonocyte_location",
                                              "SampleTypeMelanoma_VS_SampleTypeMonocyte_aggregated"))
}
)

test_that("compute the correct comparisons for multiple levels", {
  leduc <- get_data()
  leduc$Field <- as.factor(leduc$Field)
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = T, alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), c("omnibus", "Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
                                              "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"))
}
)

test_that("compute the correct comparisons for multiple levels if not adjusted", {
  leduc <- get_data()
  leduc$Field <- as.factor(leduc$Field)
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = F, alpha = 0.05)
  expect_identical(colnames(out_adj$pvals), c("omnibus_shape", "Field1_VS_Field2_shape", "Field1_VS_Field3_shape",
                                              "Field1_VS_Field4_shape", "Field2_VS_Field3_shape", "Field2_VS_Field4_shape",
                                              "Field3_VS_Field4_shape", "omnibus_location", "Field1_VS_Field2_location",
                                              "Field1_VS_Field3_location", "Field1_VS_Field4_location", "Field2_VS_Field3_location",
                                              "Field2_VS_Field4_location", "Field3_VS_Field4_location", "omnibus_aggregated",
                                              "Field1_VS_Field2_aggregated", "Field1_VS_Field3_aggregated",
                                              "Field1_VS_Field4_aggregated", "Field2_VS_Field3_aggregated",
                                              "Field2_VS_Field4_aggregated", "Field3_VS_Field4_aggregated"))
}
)

test_that("compute the correct p values", {
  leduc <- get_data()
  leduc$Field <- as.factor(leduc$Field)
  out_location <- locationTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field", nuisance = NULL)
  out_shape <- shapeTest(sce = getWithColData(leduc, "proteins_norm"), covariate = "Field")
  out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
  out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = F, alpha = 0.05)
  expect_equal(out_adj$pvals, structure(c(0.657430872482492, 0.0389495794804949, 0.00366715167911284,
                                              0.024348813975198, 0.0375201768404734, 0.117062132760355, 1.02970302066109e-08,
                                              0.645355338416144, 0.0243131618874192, 0.159018157894301, 0.44888483509382,
                                              0.120001050287176, 0.00641545137790464, 0.669519159133127, 0.0451983909682338,
                                              0.134138089825946, 7.62813598944185e-07, 0.349379221329827, 0.767411608329258,
                                              0.742855803372467, 0.30800111936782, 0.113163045483905, 0.047403821001328,
                                              0.0113892369820192, 0.0756331140299931, 0.430321146939359, 3.80553935626226e-05,
                                              0.51173691635103, 0.00529657905901098, 0.0535803976748089, 0.673357950622866,
                                              0.0752881626204939, 0.0119888090109422, 0.369601396193043, 0.536120296050784,
                                              0.713515512299861, 0.00159889993532359, 0.659545979678807, 0.276606130905013,
                                              0.681130986975393, 0.762416245849508, 0.0728335774591911, 0.755674786523898,
                                              0.00395342974616231, 0.585421722558799, 0.00918198920875205,
                                              0.198367211473923, 0.715798574445476, 0.0826809979990406, 0.0110852559978545,
                                              0.685264752415757, 0.0237907586200298, 0.119841496998611, 0.301503336094088,
                                              0.0150643453164365, 0.904295763274611, 0.0342849800069655, 0.458152518940389,
                                              0.226413636580256, 0.620595162674051, 0.330346921376695, 0.127902243617102,
                                              0.0976120298042403, 0.0449265116743442, 0.0230422254827221, 0.23172811888909,
                                              0.3422370644001, 0.490889654903771, 0.0100666932518853, 0.256297366966999,
                                              0.298515338235435, 5.52934264946103e-05, 3.22515236672344e-05,
                                              0.00148842081394629, 0.000285900613288176, 0.00160485886393992,
                                              7.03285504563089e-09, 0.343793946134101, 0.014525004218628, 0.769049142229933,
                                              0.630210627806447, 0.014262367584901, 2.74243207509011e-05, 0.51116311522195,
                                              0.00355260702911753, 0.0172333766913221, 3.63617850135757e-10,
                                              0.0828006449702246, 0.820675457866417, 0.972983070606951, 0.30111193452229,
                                              0.223573507666662, 8.68220838272122e-06, 0.00185846539337958,
                                              0.00739385418798584, 0.304513469775974, 3.16032957516972e-05,
                                              0.295294863557948, 0.00363543228815316, 0.394143303106639, 0.4800282189366,
                                              0.00474794117676338, 0.0012905396056319, 0.168048877732769, 0.281227440907842,
                                              0.140275094858713, 0.00112123607830241, 0.161938562078789, 0.14803070067154,
                                              0.960967605868819, 0.491705750029814, 0.000277463807031639, 0.28283722392383,
                                              0.000463046443024833, 0.416487881220139, 0.000411464965843538,
                                              0.0398485539035542, 0.424720025401639, 0.00710804267571234, 0.413987307197512,
                                              0.188896351385011, 5.33129060450405e-06, 0.0514884323402035,
                                              0.28049956177434, 0.000545834531271173, 0.434923612669141, 0.00757365310313541,
                                              0.594650424515641, 0.214242156709305, 0.937404064948075, 0.0641593116281649,
                                              0.0281197569675531, 0.0121158622560556, 0.00993108430389344,
                                              0.00108129734098058, 0.0126325547001453, 0.480944504759914, 0.744312677862637,
                                              0.206010186690516, 0.433821455135794, 0.410594648663631, 0.000110430084729993,
                                              6.39407075613526e-05, 0.00280535295757099, 0.00056747709848988,
                                              0.00316630933034802, 8.35753033439263e-09, 0.448606215296766,
                                              0.0181856567593166, 0.263542908819706, 0.524313192786383, 0.0254946450328052,
                                              5.4615176383484e-05, 0.579721583895562, 0.00658743935733475,
                                              0.0305427737976129, 7.26889206854185e-10, 0.133874005343969,
                                              0.793150308246409, 0.842487172361624, 0.304517567947597, 0.150267375373211,
                                              1.7361236982166e-05, 0.0031954979344023, 0.0134708090376828,
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
                                                                                                                 21L), dimnames = list(c("A0AVT1", "A2A368", "A6H8Y1", "A6NCS6",
                                                                                                                                         "A6NLC8", "B0I1T2", "E9PAV3", "E9PRG8", "O00154", "O00159"),
                                                                                                                                       c("omnibus_shape", "Field1_VS_Field2_shape", "Field1_VS_Field3_shape",
                                                                                                                                         "Field1_VS_Field4_shape", "Field2_VS_Field3_shape", "Field2_VS_Field4_shape",
                                                                                                                                         "Field3_VS_Field4_shape", "omnibus_location", "Field1_VS_Field2_location",
                                                                                                                                         "Field1_VS_Field3_location", "Field1_VS_Field4_location",
                                                                                                                                         "Field2_VS_Field3_location", "Field2_VS_Field4_location",
                                                                                                                                         "Field3_VS_Field4_location", "omnibus_aggregated", "Field1_VS_Field2_aggregated",
                                                                                                                                         "Field1_VS_Field3_aggregated", "Field1_VS_Field4_aggregated",
                                                                                                                                         "Field2_VS_Field3_aggregated", "Field2_VS_Field4_aggregated",
                                                                                                                                         "Field3_VS_Field4_aggregated"))))
}
)
