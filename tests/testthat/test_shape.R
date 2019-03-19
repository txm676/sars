context("sar_shape")

test_that("observed_shape info is correct", {
  data(galap)
  fit <- sar_epm2(galap)
  fit2 <- sar_power(galap)
  x <- substr(fit$observed_shape, 1, 41)
  expect_match(x, "sigmoid - observed shape algorithm failed")
  expect_match(fit2$observed_shape, "convex up")
})
