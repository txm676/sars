
context("sar_threshold")

test_that("sar_threshold returns correct results", {
  data(aegean2)
  fit <- sar_threshold(aegean2, mod = c("ContOne", "DiscOne", "ZslopeOne"),
                       non_th_models = TRUE, interval = 0.01, logAxes = "area")
  s <- summary(fit)
  expect_equal(c(s$Model_table$BIC), c(2061.57, 2066.70, 2074.14, 
                                       2314.51, 2542.14))
  expect_equal(s$Model_table$R2[1], 0.94)
  expect_equal(round(s$Model_table$Th1[3], 2), 1.37)
  expect_is(fit, "threshold")
  expect_match(fit[[5]][[1]], "area")
  expect_error(sar_threshold(aegean2, mod = c("D")), 
               "Incorrect model names provided; see help for 'mod' argument")
  expect_error(sar_threshold(aegean2, interval = 30000), 
               "interval must be smaller than max area")
  a2 <- aegean2[1:169,]
  fit2 <- sar_threshold(a2, mod = c("ContTwo", "DiscTwo", "ZslopeTwo"),
                        non_th_models = TRUE, interval = 1, logAxes = "area")
  s2 <- summary(fit2)
  expect_equal(s2$`Axes transformation`, "area")
  expect_equal(s2$Model_table$AIC, c(1974.90, 1974.87, 1975.10, 2195.69,
                                     2426.60))
  expect_equal(length(s2$Model_table$Th2[!is.na(s2$Model_table$Th2)]), 3)
})

