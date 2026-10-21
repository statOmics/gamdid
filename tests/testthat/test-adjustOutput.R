test_that("compute a matrix", {
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, alpha = 0.05)
    expect_true(is(out_adj$pvals, "matrix"))
})

test_that("compute the correct comparisons", {
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, alpha = 0.05)
    expect_identical(colnames(out_adj$pvals), "SampleTypeMelanoma_VS_SampleTypeMonocyte")
})

test_that("compute the correct comparisons if not adjusted", {
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "SampleType")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "none", alpha = 0.05)
    expect_identical(colnames(out_adj$pvals), c(
        "SampleTypeMelanoma_VS_SampleTypeMonocyte_GAM", "SampleTypeMelanoma_VS_SampleTypeMonocyte_msqrob2",
        "SampleTypeMelanoma_VS_SampleTypeMonocyte_aggregated"
    ))
})

test_that("compute the correct comparisons for multiple levels", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, alpha = 0.05)
    expect_identical(colnames(out_adj$pvals), c(
        "omnibus", "Field1_VS_Field2", "Field1_VS_Field3", "Field1_VS_Field4",
        "Field2_VS_Field3", "Field2_VS_Field4", "Field3_VS_Field4"
    ))
})

test_that("compute the correct comparisons for multiple levels if not adjusted", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "none", alpha = 0.05)
    expect_identical(colnames(out_adj$pvals), c(
        "omnibus_GAM", "Field1_VS_Field2_GAM", "Field1_VS_Field3_GAM",
        "Field1_VS_Field4_GAM", "Field2_VS_Field3_GAM", "Field2_VS_Field4_GAM",
        "Field3_VS_Field4_GAM", "omnibus_msqrob2", "Field1_VS_Field2_msqrob2",
        "Field1_VS_Field3_msqrob2", "Field1_VS_Field4_msqrob2", "Field2_VS_Field3_msqrob2",
        "Field2_VS_Field4_msqrob2", "Field3_VS_Field4_msqrob2", "omnibus_aggregated",
        "Field1_VS_Field2_aggregated", "Field1_VS_Field3_aggregated",
        "Field1_VS_Field4_aggregated", "Field2_VS_Field3_aggregated",
        "Field2_VS_Field4_aggregated", "Field3_VS_Field4_aggregated"
    ))
})

test_that("compute the correct p values", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "none", alpha = 0.05)
    expect_equal(out_adj$pvals, structure(c(
        NA, NA, NA, NA, NA, NA, 1.0297030206397649e-08, 0.64535533841612813,
        0.024313161913675579, 0.15901815692459656, NA, NA, NA, NA, NA,
        NA, 7.6281359893748348e-07, 0.38459751652700552, 0.76997001840501134,
        0.75161608560716753, NA, NA, NA, NA, NA, NA, 6.4188030178802195e-05,
        0.57832755185543605, 0.0064326121194358387, 0.11249142417069111,
        NA, NA, NA, NA, NA, NA, 0.002365049796182672, 0.65143948581624134,
        0.25062563190762449, 0.58977768920183138, NA, NA, NA, NA, NA,
        NA, 0.16039976923272239, 0.77360604140078981, 0.09485993654556267,
        0.027790911106852702, NA, NA, NA, NA, NA, NA, 0.024470192906626866,
        0.44995899362891334, 0.20354134535371987, 0.52688131021629325,
        NA, NA, NA, NA, NA, NA, 0.34223706440009616, 0.4908896549037613,
        0.010066693238340377, 0.25629736685539473, 0.29851533823552934,
        5.5293426493848227e-05, 3.2251523667522648e-05, 0.0014884208139434043,
        0.00028590061328926509, 0.0016048588639345431, 7.0328550456368151e-09,
        0.3437939461339572, 0.01452500421863427, 0.7690491422298934,
        0.63021062780648529, 0.014262367584840041, 2.742432075111087e-05,
        0.51116311522189939, 0.0035526070291253585, 0.017233376691296542,
        3.6361785013606404e-10, 0.082800644970175696, 0.82067545786642504,
        0.97298307060694933, 0.30111193452234197, 0.22357350766637379,
        8.6822083827939771e-06, 0.0018584653933770296, 0.0073938541880000279,
        0.30451346977586496, 3.1603295751709908e-05, 0.29529486355787188,
        0.003635432288154643, 0.39414330310660384, 0.48002821893664904,
        0.004747941176736996, 0.001290539605638285, 0.16804887773271701,
        0.28122744090796603, 0.14027509485862266, 0.0011212360783026986,
        0.1619385620787224, 0.14803070067155866, 0.96096760586881635,
        0.49170575002986161, 0.00027746380702914353, 0.28283722392405936,
        0.00046304644302402518, 0.41648788122025993, 0.00041146496584225486,
        0.039848553903558716, 0.42472002540156695, 0.0071080426757148896,
        0.41398730719747812, 0.18889635138505723, 5.3312906044293969e-06,
        0.051488432340312223, 0.28049956177428137, 0.00054583453127276019,
        0.43492361266903712, 0.0075736531031367895, 0.59465042451558325,
        0.21424215670932606, 0.93740406494807083, 0.064159311628191121,
        0.02811975696745456, 0.012115862256094845, 0.0099310843038840184,
        0.0010812973409834495, 0.012632554700124884, 0.48094450475992434,
        0.74431267786259803, 0.20601018669053589, 0.43382145513575965,
        NA, NA, NA, NA, NA, NA, 8.3575303343265779e-09, 0.44860621529663919,
        0.018185656766666329, 0.26354290748796699, NA, NA, NA, NA, NA,
        NA, 7.2688920685479539e-10, 0.13626464563042967, 0.79451456257764164,
        0.84809240948143216, NA, NA, NA, NA, NA, NA, 4.2353590614914295e-05,
        0.39096330978667698, 0.0046454554329421307, 0.17502843412285474,
        NA, NA, NA, NA, NA, NA, 0.0015212631745834093, 0.25939518256584998,
        0.18612667034540964, 0.73094821669126109, NA, NA, NA, NA, NA,
        NA, 0.063837726581496307, 0.54837491507015745, 0.013225102278790028,
        0.052085340458161059, NA, NA, NA, NA, NA, NA, 0.011567197794244883,
        0.51228392531882438, 0.20875470952186209, 0.67459620961723954,
        NA, NA, NA, NA, NA, NA, 0.3999045693317706, 0.59160411842514971,
        0.019195402618452102, 0.32222652996441853
    ), dim = c(10L, 21L), dimnames = list(
        c(
            "A0AVT1", "A2A368", "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2",
            "E9PAV3", "E9PRG8", "O00154", "O00159"
        ), c(
            "omnibus_GAM",
            "Field1_VS_Field2_GAM", "Field1_VS_Field3_GAM", "Field1_VS_Field4_GAM",
            "Field2_VS_Field3_GAM", "Field2_VS_Field4_GAM", "Field3_VS_Field4_GAM",
            "omnibus_msqrob2", "Field1_VS_Field2_msqrob2", "Field1_VS_Field3_msqrob2",
            "Field1_VS_Field4_msqrob2", "Field2_VS_Field3_msqrob2", "Field2_VS_Field4_msqrob2",
            "Field3_VS_Field4_msqrob2", "omnibus_aggregated", "Field1_VS_Field2_aggregated",
            "Field1_VS_Field3_aggregated", "Field1_VS_Field4_aggregated",
            "Field2_VS_Field3_aggregated", "Field2_VS_Field4_aggregated",
            "Field3_VS_Field4_aggregated"
        )
    )))
})

test_that("compute the correct p values contrastwise", {
    leduc_testdata$Field <- as.factor(leduc_testdata$Field)
    out_location <- locationTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field", nuisance = NULL)
    out_shape <- shapeTest(sce = getWithColData(leduc_testdata, "proteins_norm"), covariate = "Field")
    out <- aggregatePvals(pval_shape = out_shape$pvals, pval_location = out_location$pvals)
    out_adj <- adjustOutput(pval = out, pval_shape = out_shape$pvals, pval_location = out_location$pvals, adjust = "contrastwise", alpha = 0.05)
    expect_equal(out_adj$pvals, structure(c(
        NA, NA, NA, NA, NA, NA, 3.3430121337306312e-08, 0.44860621529663919,
        0.036371313533332658, 0.35139054331728931, NA, NA, NA, NA, NA,
        NA, 2.9075568274191816e-09, 0.27252929126085934, 0.84809240948143216,
        0.84809240948143216, NA, NA, NA, NA, NA, NA, 0.00016941436245965718,
        0.39096330978667698, 0.0092909108658842613, 0.23337124549713964,
        NA, NA, NA, NA, NA, NA, 0.0060850526983336371, 0.34586024342113331,
        0.34586024342113331, 0.73094821669126109, NA, NA, NA, NA, NA,
        NA, 0.0851169687753284, 0.54837491507015745, 0.052900409115160113,
        0.0851169687753284, NA, NA, NA, NA, NA, NA, 0.046268791176979532,
        0.67459620961723954, 0.41750941904372418, 0.67459620961723954,
        NA, NA, NA, NA, NA, NA, 0.53320609244236072, 0.59160411842514971,
        0.076781610473808409, 0.53320609244236072
    ), dim = c(10L, 7L), dimnames = list(
        c(
            "A0AVT1", "A2A368", "A6H8Y1", "A6NCS6", "A6NLC8", "B0I1T2",
            "E9PAV3", "E9PRG8", "O00154", "O00159"
        ), c(
            "omnibus", "Field1_VS_Field2",
            "Field1_VS_Field3", "Field1_VS_Field4", "Field2_VS_Field3",
            "Field2_VS_Field4", "Field3_VS_Field4"
        )
    )))
})

test_that("adj_alpha is the BH threshold at the requested alpha", {
    p <- c(0.001, 0.008, 0.012, 0.03, 0.2, 0.6, NA, 0.9, 0.04, 0.5)
    pval <- matrix(p, ncol = 1, dimnames = list(NULL, "A_VS_B"))
    m <- sum(!is.na(p))
    ## Worked out by hand: 1, 3 and 5 BH rejections at alpha 0.01, 0.05, 0.2
    for (case in list(c(0.01, 1), c(0.05, 3), c(0.2, 5))) {
        alpha <- case[1]
        out <- adjustOutput(
            pval = pval, pval_shape = pval, pval_location = pval,
            adjust = "contrastwise", alpha = alpha
        )
        expect_equal(out$adj_alpha, case[2] / m * alpha)
        ## Raw p-values below adj_alpha are exactly the BH discoveries
        expect_identical(
            which(p <= out$adj_alpha),
            which(p.adjust(p, method = "BH") <= alpha)
        )
    }
    ## NA when nothing is significant
    out <- adjustOutput(
        pval = pval, pval_shape = pval, pval_location = pval,
        adjust = "contrastwise", alpha = 0.001
    )
    expect_true(is.na(out$adj_alpha))
})

test_that("compDistr stores the BH threshold for its alpha", {
    sce <- getWithColData(leduc_testdata, "proteins_norm")
    for (alpha in c(0.01, 0.2)) {
        out <- compDistr(
            object = sce, covariate = "SampleType", adjust = "none",
            alpha = alpha
        )
        p <- rowData(out)$SampleTypeMelanoma_VS_SampleTypeMonocyte_aggregated
        rejected <- which(p.adjust(p, method = "BH") <= alpha)
        expect_equal(
            metadata(out)$adj_alpha,
            length(rejected) / sum(!is.na(p)) * alpha
        )
        expect_identical(which(p <= metadata(out)$adj_alpha), rejected)
    }
})
