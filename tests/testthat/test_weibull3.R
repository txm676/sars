context("sar_weibull3")

test_that("sar_weibull3 returns correct results", {
  skip_on_cran()
  fit <- sar_weibull3(galap)
  expect_equal(round(fit$AICc, 2), 190.32)
  expect_equal(as.vector(round(fit$par[2], 2)), 0.03)
  expect_is(fit, "sars")
  expect_match(fit$normaTest[[1]], "none")
  expect_error(sar_linear(5), "data must be a matrix or dataframe")
})


test_that("sar_weibull3 summary returns correct results", {
  skip_on_cran()
  fit <- sar_weibull3(galap, normaTest = "lillie")
  fs <- summary(fit)
  expect_equal(round(sum(fs$residuals),1), 77.2)
  expect_output(str(fs), "List of 16")
  expect_is(fs, "summary.sars")
  expect_equal(round(fs$normaTest[[2]]$p.value, 3), 0.137)
})

test_that("gompertz summary returns correct results", {
  skip_on_cran()
  fit <- sar_gompertz(galap, normaTest = "lillie")
  fs <- summary(fit)
  expect_equal(length(capture_output_lines(fit, print = TRUE)),
               11)
  expect_equal(length(capture_output_lines(fs, print = TRUE)),
                      23)
  #one in every 1000 runs is 3 (due to stochastic nature 
  #of grid_start), so using %in% rather than =
  expect_true(round(sum(fs$residuals),0) %in% c(2,3))
  expect_output(str(fs), "List of 16")
  expect_is(fs, "summary.sars")
  expect_equal(round(fs$normaTest[[2]]$p.value, 2), 0.05)
  #start pars
  fit3 <- sar_gompertz(galap, start = c(210.88, 0.0168,38.04),
                      grid_start = "none")
  expect_equal(length(capture_output_lines(fit3, print = TRUE)),
               11)
  expect_equal(round(fit$par[1],1), 
               round(fit3$par[1]),1)
  
})
