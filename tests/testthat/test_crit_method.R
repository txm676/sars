context("new_information_criteria")
library(sars)

test_that("new ICs works for individual models", {
  s2 <- sar_loga(galap)
  l <- c(s2$AIC, s2$BIC, s2$AICc)
  expect_equal(round(l,4), c(187.1908, 189.5086, 189.1908))
  #match with nls
  y <- galap$s
  x <- galap$a
  n <- stats::nls(y ~ c + z*log(x), start = list("c" = 0.29, "z" = 30))
  nl <- c(stats::AIC(n), stats::BIC(n))
  expect_equal(round(nl,4), c(187.1908, 189.5086))
  #test linear model
  s2 <- sar_linear(galap)
  l <- c(s2$AIC, s2$BIC, s2$AICc)
  expect_equal(round(l,4), c(191.7655, 194.0832, 193.7655))
})


test_that("news ICs works for multi models", {
  fit3 <- sar_multi(galap, obj = c("power", "loga", "linear"))
  ff <- c(fit3$power$AIC, fit3$power$BIC, fit3$power$AICc)
  expect_equal(round(ff,4), c(187.0308, 189.3486, 189.0308))
  fit4 <- sar_average(data = galap)
  ff4 <- c(fit4$details$fits$koba$AIC, fit4$details$fits$koba$BIC, 
           fit4$details$fits$koba$AICc)
  expect_equal(round(ff4,3), c(186.073, 188.391, 188.073))
})
