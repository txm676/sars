context("lin_pow")

test_that("lin_pow returns correct results", {
   fit <- lin_pow(galap, con = 1)
   expect_equal(round(fit$Model$coefficients[2], 2), 0.34)
   expect_equal(round(fit$normaTest[[2]]$p.value, 2), 0.35)
})

test_that("lin_pow log-transformation works",{
  fit <- lin_pow(galap, con = 1)
  fit2 <- lin_pow(galap, con = 1, logT = log10)
  fit3 <- summary(lin_pow(galap, logT = log2))
  expect_match(fit$logT, "log()")
  expect_match(fit2$logT, "log10()")
  expect_match(fit3$logT, "log2()")
  expect_equal(round(fit$Model$coefficients[1], 2), 3.02)
  expect_equal(round(fit2$Model$coefficients[1], 2), 1.31)
  expect_equal(round(fit$Model$coefficients[2], 2), 
               round(fit2$Model$coefficients[2], 2), 
               round(fit3$Model$coefficients[2], 2))
})
