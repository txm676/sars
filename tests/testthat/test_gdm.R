context("gdm")
library(sars)

test_that("gdm functions return correct results", {
  data(galap)
  galap$t <- rgamma(16, 5, scale = 2)
  g <- gdm(galap, model = "expo", mod_sel = FALSE)
  expect_equal(round(as.vector(g$m$getPars()[2]), 2), 32.94)
  expect_is(g, "gdm")
})







