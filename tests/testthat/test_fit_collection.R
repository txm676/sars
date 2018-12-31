context("fit_collection")
library(sars)

test_that("fit_collection returns correct results", {
  data(galap)
  fit <- sar_linear(galap)
  fit2 <- sar_power(galap)
  fitC <- fit_collection(fit, fit2)
  expect_output(str(fitC), "List of 2")
  expect_is(fitC, "sars")
})





