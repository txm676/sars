context("sar_logistic")
library(sars)

test_that("sar_logistic returns correct results", {
  fit <- sar_logistic(galap)
  expect_equal(round(fit$AICc, 2), 189.21)
  expect_equal(as.vector(round(fit$par[2], 2)), 0.02)
  expect_is(fit, "sars")
  expect_match(fit$normaTest[[1]], "none")
  expect_match(fit$homoTest[[1]], "none")
  expect_error(sar_linear(5), "data must be a matrix or dataframe")
  fit2 <- sar_logistic(galap, homoTest = "cor.area", homoCor = "kendall")
  expect_equal(round(fit2$homoTest[[2]]$p.value, 2), 0.17)
  expect_match(fit2$homoTest[[2]]$method, "Kendall's rank correlation tau")
  #check it works with 0 richness values
  g2 <- galap
  g2$s[1] <- 0
  fit3 <- sar_logistic(g2)
  expect_equal(round(fit3$BIC, 2), 188.92)
})

test_that("sar_mmf returns correct results", {
  expect_warning(sar_mmf(galap, grid_start = "none"), 
                 "'sar_mmf' is deprecated.")
  s1 <- suppressWarnings(sar_mmf(niering))#deprecation warning
  s2 <- sar_heleg(niering)
  expect_equal(s1$AICc, s2$AICc)
})