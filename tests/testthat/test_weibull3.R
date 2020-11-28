context("sar_weibull3")

test_that("sar_weibull3 returns correct results", {
  fit <- sar_weibull3(galap)
  expect_equal(round(fit$AICc, 2), 190.32)
  expect_equal(as.vector(round(fit$par[2], 2)), 0.03)
  expect_is(fit, "sars")
  expect_match(fit$normaTest[[1]], "none")
  expect_error(sar_linear(5), "data must be a matrix or dataframe")
})


test_that("sar_weibull3 summary returns correct results", {
  fit <- sar_weibull3(galap, normaTest = "lillie")
  fs <- summary(fit)
  expect_equal(round(sum(fs$residuals),1), 77.2)
  expect_output(str(fs), "List of 16")
  expect_is(fs, "summary.sars")
  expect_equal(round(fs$normaTest[[2]]$p.value, 3), 0.137)
})
