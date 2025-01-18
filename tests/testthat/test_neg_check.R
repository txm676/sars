context("neg_check")
library(sars)

test_that("neg_check returns correct results for individual model fits", {
  data(galap)
  data(niering)
  fit <- sar_loga(galap)
  expect_true(fit$neg_check)
  expect_true(summary(fit)$Negative_values)
  fit2 <- sar_power(galap)
  expect_false(fit2$neg_check)
  fit3 <- sar_linear(niering)
  expect_false(fit3$neg_check)
})
