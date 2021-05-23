context("convergence")
library(sars)

test_that("various functions return correct convergence info", {
  fit <- sar_power(galap)
  c1 <- fit$convergence
  c2 <- fit$verge
  expect_equal(c1, 0)
  expect_true(c2)
  
  fit2 <- sar_average(data = galap, normaTest = "lillie", grid_start = "none")
  f1 <- fit2$details$convergence
  expect_equal(length(f1), length(fit2$details$fits))
  expect_identical(names(f1), names(fit2$details$fits))
  
  s2 <- summary(fit2)
  expect_equal(length(s2), 5)
  expect_identical(f1, s2$Convergence)
  expect_true(all(s2$Convergence))
  
  #created this by just randomly creating data until one of the model's did
  #not fully converge (i.e. produced a fit but with optim code != 0)
  test <- data.frame("a" =  c(0.52,  2.33,  2.59,  4.66,  4.84, 11.40, 18.39),
                     "s" = c(7.458806, 15.904833, 66.768317, 44.708306,
                             82.288296, 57.104797, 29.598247))
  #set verb to FALSE, as it produces warnings of convergence etc
  s3 <- sar_average(data = test, grid_start = "none", verb = FALSE)
  expect_false(all(s3$details$convergence))
  expect_false(all(summary(s3)$Convergence))
  expect_false(s3$details$fits$chapman$verge)
  expect_equal(s3$details$fits$chapman$convergence, 10)
})
