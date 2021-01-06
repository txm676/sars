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
  fit2 <- sar_multi(galap, obj = c("power", "koba", "ratio"), 
                    grid_start = "none")
  p2 <- sar_pred(fit2, area = c(5000, 50000))
  expect_equal(nrow(p2), 6)
  expect_equal(round(p2$Prediction[4], 2), 382.84)
  expect_match(as.character(p2$Model[3]), "Kobayashi")
  expect_error(sar_pred(fit2, area = "a"))
})

test_that("sar_pred on multi-model curve returns correct results", {
  skip_on_cran()
  #grid_start has a random component and so it is hard to use it to test,
  #as it can produce different sar_pred values each time it is run. So here
  #I have used grid_start = none except for the last.
  data(niering)
  fit3 <- sar_average(data = niering, grid_start = "none")
  p3 <- sar_pred(fit3, area = c(50, 500))
  expect_equal(nrow(p3), 2)
  expect_equal(round(p3$Prediction[2], 2), 53.76)
  expect_match(as.character(p3$Model[1]), "Multi")
  fit4 <- sar_average(data = niering, normaTest = "lillie",
                      homoTest = "cor.fitted", grid_start = "none")
  p4 <- sar_pred(fit4, area = c(50, 500))
  expect_equal(round(p4$Prediction[2], 2), 105.57)
  fit5 <- sar_average(data = niering, normaTest = "lillie", grid_start = "none")
  p5 <- sar_pred(fit5, area = c(50, 500))
  expect_equal(round(p5$Prediction[2], 2), 47.95)
  #test it using grid_start = exhaustive, which for this always seems to give the
  #same answer (in contrast to = partial). But exclude more complex models for
  #speed
  obj = c("linear","power","powerR","epm1","epm2","p1",
                "p2","loga","koba","mmf","monod","negexpo", 
                "weibull3","asymp","ratio",
        "weibull4", "heleg")
  fit5 <- sar_average(data = niering, obj = obj, normaTest = "lillie", 
                      grid_start = "exhaustive", grid_n = 1000)
  p5 <- sar_pred(fit5, area = c(50, 500))
  expect_equal(round(p5$Prediction[2], 2), 17.07)
  #nb. the expected richness for 500 here is lower than when no grid_start is
  #used and this is because grid_start improves the epm1 model which becomes
  #the best model but which starts to decrease toward the final richness point,
  #and thus so does the mmi curve.
})











