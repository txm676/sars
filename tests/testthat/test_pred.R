context("sar_pred")
library(sars)

test_that("sar_pred on ind model fit returns correct results", {
  data(galap)
  fit <- sar_loga(galap)
  p <- sar_pred(fit, area = 5000)
  expect_is(p, "sars")
  expect_true(is.data.frame(p))
  expect_equal(p$Area, 5000)
  expect_equal(round(p$Prediction, 2), 258.19)
  expect_match(as.character(p$Model), "Logarithmic")
})

test_that("sar_pred on fit_collection returns correct results", {
  data(galap)
  fit2 <- sar_multi(galap, obj = c("power", "koba", "ratio"))
  p2 <- sar_pred(fit2, area = c(5000, 50000))
  expect_equal(nrow(p2), 6)
  expect_equal(round(p2$Prediction[4], 2), 382.84)
  expect_match(as.character(p2$Model[3]), "Kobayashi")
  expect_error(sar_pred(fit2, area = "a"))
})

test_that("sar_pred on multi-model curve returns correct results", {
  data(niering)
  fit3 <- sar_average(data = niering)
  p3 <- sar_pred(fit3, area = c(50, 500))
  expect_equal(nrow(p3), 2)
  expect_equal(round(p3$Prediction[2], 2), 46.82)
  expect_match(as.character(p3$Model[1]), "Multi")
})


