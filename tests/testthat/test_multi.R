context("sar_multi")
library(sars)

test_that("sar_multi returns correct results", {
  data("galap")
  fit3 <- sar_multi(galap, normaTest = "none", homoTest = "none", 
                    neg_check = FALSE)
  expect_equal(round(sum(fit3$mmi), 2), 1640.88)
  expect_output(str(fit3), "List of 2")
  expect_is(fit3, "multi")
  expect_match(fit3$details$homo_test, "none")
  expect_match(fit3$details$norm_test, "none")
  expect_match(fit3$details$ic, "AICc")
  expect_error(sar_multi(5), "argument is of length zero")
})


test_that("confidence intervals are correct", {
  skip_on_cran()
  data("galap")
  fit3 <- sar_multi(galap, normaTest = "none", homoTest = "none", 
                    neg_check = FALSE, confInt = TRUE, ciN = 20)
  ci <- fit3$details$confInt
  expect_equal(nrow(ci), nrow(galap))
  expect_output(str(ci), "data.frame")
  expect_true(all(ci[ ,1] < ci[ ,2]))
})

