test_that("a group needs at least minObs observations", {
    colDat <- data.frame(group = rep(c("A", "B"), each = 30))
    set.seed(1)
    x <- rnorm(60)
    x[c(1:5, 31:40)] <- NA # 25 observations in A, 20 in B
    expect_false(all(is.na(createCounts(x, colDat, "group", minObs = 20))))
    expect_true(all(is.na(createCounts(x, colDat, "group", minObs = 21))))
})
