context("sar_multi")

test_that("sar_multi returns correct results", {
  fitC <- sar_multi(galap, obj = c("linear", "power"))
  expect_output(str(fitC), "List of 2")
  expect_is(fitC, "sars")
  expect_no_error(plot(fitC))
  expect_no_error(plot(fitC, mfplot = TRUE))
  expect_no_error(plot(fitC, mfplot = TRUE,
                       pLeg = TRUE))
})
