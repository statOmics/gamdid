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
        0.010066693238340377, 0.25629736685539473, 0.41410615103094128,
        0.00029595236695501947, 5.3902881659619701e-05, 0.0057881027173625772,
        0.0008735600924925902, 0.0039616169349377594, 6.5996757706715981e-08,
        0.40673204051162282, 0.031192707427080611, 0.81521672133592848,
        0.67415136251196994, 0.024487743403602848, 3.4648277110716553e-05,
        0.55221722731641065, 0.0063121275333825324, 0.026451434087710252,
        3.4637152309322939e-09, 0.10455068249207686, 0.83622632233130889,
        0.97530212821079876, 0.36615627859844402, 0.26284281392854131,
        1.4010329001446725e-05, 0.0050616268681479672, 0.012752447217665943,
        0.33819064635281948, 9.3678946734533319e-05, 0.32837571933252441,
        0.0077245133418806485, 0.43634381781653697, 0.53510944560275964,
        0.010218880289397245, 0.0015124432273819364, 0.21355192621939792,
        0.31713012860714757, 0.16941398073540945, 0.0023382539811498126,
        0.1927789381270909, 0.18967285653634788, 0.96439628008650435,
        0.54709280689752848, 0.00078790821870622417, 0.31103111389314947,
        0.00154853252192313, 0.4450387606625571, 0.00096626586544811834,
        0.053541427247293133, 0.45388317109209253, 0.013440132041621287,
        0.456084911506486, 0.24672104460340022, 3.0666277480732436e-05,
        0.064404583751970737, 0.32630085708502998, 0.0011890010910852384,
        0.46680052815952189, 0.012540871351346317, 0.61847585336356536,
        0.25900311214231325, 0.94294165404194197, 0.10252310652448822,
        0.045476203910963346, 0.018175286018934995, 0.020513109785148167,
        0.0022896403648831731, 0.019880665085225187, 0.51129430732375947,
        0.76066422588667426, 0.2486917711975577, 0.47603435214225598,
        NA, NA, NA, NA, NA, NA, 1.7814572489290558e-08, 0.49898268700186771,
        0.027326600058983585, 0.26612527104841677, NA, NA, NA, NA, NA,
        NA, 6.8961171945156473e-09, 0.16440797663485862, 0.80173161956156258,
        0.84897218875173808, NA, NA, NA, NA, NA, NA, 7.6178909328397806e-05,
        0.41889939495096046, 0.0070196168389546317, 0.17886948118175355,
        NA, NA, NA, NA, NA, NA, 0.0023515755576724059, 0.29751497661733772,
        0.21593024174642575, 0.73193789213276772, NA, NA, NA, NA, NA,
        NA, 0.080284047356523322, 0.57210566035671317, 0.023544400096216594,
        0.052389537317230409, NA, NA, NA, NA, NA, NA, 0.016583016313481178,
        0.52092792245892039, 0.22794713477546721, 0.67602472708658867,
        NA, NA, NA, NA, NA, NA, 0.41002327173880804, 0.59670175623173194,
        0.019350120793173407, 0.33319969030909974
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
        NA, NA, NA, NA, NA, NA, 7.125828995716223e-08, 0.49898268700186771,
        0.05465320011796717, 0.35483369473122234, NA, NA, NA, NA, NA,
        NA, 2.7584468778062589e-08, 0.32881595326971724, 0.84897218875173808,
        0.84897218875173808, NA, NA, NA, NA, NA, NA, 0.00030471563731359122,
        0.41889939495096046, 0.014039233677909263, 0.2384926415756714,
        NA, NA, NA, NA, NA, NA, 0.0094063022306896237, 0.39668663548978361,
        0.39668663548978361, 0.73193789213276772, NA, NA, NA, NA, NA,
        NA, 0.10704539647536443, 0.57210566035671317, 0.094177600384866375,
        0.10477907463446082, NA, NA, NA, NA, NA, NA, 0.06633206525392471,
        0.67602472708658867, 0.45589426955093443, 0.67602472708658867,
        NA, NA, NA, NA, NA, NA, 0.54669769565174398, 0.59670175623173194,
        0.077400483172693629, 0.54669769565174398
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
