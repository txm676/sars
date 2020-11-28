context("sar_monod")

test_that("sar_monod returns correct results", {
  data(niering)
  fit <- sar_monod(niering)
  expect_equal(nrow(niering), 32)
  expect_equal(round(fit$AICc, 2), 204.18)
  expect_equal(as.vector(round(fit$par[2], 2)), 0.04)
  expect_is(fit, "sars")
  expect_match(fit$normaTest[[1]], "none")
  expect_error(sar_linear(5), "data must be a matrix or dataframe")
  #checking homoCor argument
  expect_error(sar_monod(niering, homoTest = "cor.fitted",
                         homoCor = "4"))
  expect_error(sar_linear(niering, homoTest = "cor.fitted",
                         homoCor = "4"))
  fit2 <- sar_monod(niering, homoTest = "none", homoCor = "spearman")
  expect_equal(fit2$homoTest$test, "none")
  fit3 <- sar_monod(niering, homoTest = "cor.fitted", homoCor = "pearson")
  expect_equal(round(fit3$homoTest[[2]]$p.value, 2), 0.01)
  expect_equal(fit3$homoTest[[2]]$method, 
               "Pearson's product-moment correlation")
  fit4 <- sar_linear(niering, homoTest = "cor.fitted", homoCor = "pearson")
  expect_equal(round(fit4$homoTest[[2]]$p.value, 2), 0.00)
  fit5 <- sar_monod(niering, normaTest = "lillie")
  expect_match(fit5$normaTest[[1]], "lillie")
  #checking homotest matches cor.test done manually
  fit2 <- sar_monod(aegean, homoTest = "cor.fitted", homoCor = "pearson")
  rs2 <- fit2$residuals ^ 2
  pp <- cor.test(rs2, fit2$calculated)
  expect_equal(fit2$homoTest[[2]]$p.value, pp$p.value)
})


test_that("sar_monod summary returns correct results", {
  data(aegean)
  fit <- sar_monod(aegean)
  fs <- summary(fit)
  expect_equal(nrow(aegean), 90)
  expect_equal(round(sum(fs$residuals),4), 161.7902)
  expect_output(str(fs), "List of 16")
  expect_is(fs, "summary.sars")
  fit2 <- summary(sar_monod(aegean, normaTest = "lillie"))
  expect_equal(round(fit2$normaTest[[2]]$p.value, 3), 0.025)
})
