context("sar_linear")

test_that("sar_linear returns correct results", {
  fit <- sar_linear(galap)
  expect_equal(round(fit$AICc, 2), 193.77)
  expect_equal(as.vector(round(fit$par[2], 2)), 0.19)
  expect_is(fit, "sars")
  expect_match(fit$normaTest[[1]], "none")
  expect_error(sar_linear(5), "data must be a matrix or dataframe")
  fit2 <- sar_linear(galap, normaTest = "shapiro", homoTest = "cor.fitted",
                     homoCor = "pearson")
  expect_match(fit2$normaTest[[1]], "shapiro")
  expect_match(fit2$homoTest[[1]], "cor.fitted")
  expect_equal(as.vector(round(fit2$normaTest[[2]]$p.value, 2)), 0.02)
  expect_equal(as.vector(round(fit2$homoTest[[2]]$p.value, 2)), 0.44)
  expect_match(fit2$homoTest[[2]]$method, 
               "Pearson's product-moment correlation")
})


test_that("sar_linear summary returns correct results", {
  fit <- sar_linear(galap, normaTest = "lillie")
  fs <- summary(fit)
  expect_equal(sum(round(fs$residuals, 1)), 0.2)
  expect_output(str(fs), "List of 16")
  expect_is(fs, "summary.sars")
  expect_equal(round(fs$normaTest[[2]]$p.value, 3), 0.047)
  #checking homotest matches cor.test done manually
  fit2 <- sar_linear(aegean, homoTest = "cor.fitted", homoCor = "kendall")
  rs2 <- summary(fit2)$residuals^2
  pp <- cor.test(rs2, fit2$calculated, method = "kendall")
  expect_equal(fit2$homoTest[[2]]$p.value, pp$p.value)
})
