context("sar_average")

test_that("sar_average returns correct results", {
  skip_on_cran()

  fit2 <- sar_average(data = galap, grid_start = "none")
  expect_equal(round(sum(fit2$mmi), 1), 1647.2)
  fit3 <- sar_average(data = galap)#grid start on so round to 0
  expect_equal(round(sum(fit3$mmi), 0), 1647)
  expect_output(str(fit3), "List of 2")
  expect_is(fit3, "multi")
  expect_match(fit3$details$homo_test, "none")
  expect_match(fit3$details$norm_test, "none")
  expect_match(fit3$details$ic, "AICc")
  expect_error(sar_multi(5), "argument is of length zero")
  fit4 <- sar_average(data = galap, normaTest = "lillie", 
                      homoTest = "cor.fitted",
                      neg_check = FALSE)
  #one run in 100 gave 1638 (due to random nature of grid_start)
  expect_true(round(sum(fit4$mmi), 0) %in% c(1638, 1639))#grid_start on so round to 0
  expect_match(fit4$details$norm_test, "lillie")
  expect_match(fit4$details$homo_test, "cor.fitted")
  expect_equal(length(fit4$details$mod_names), 14)
  expect_error(sar_average(data = galap, homoTest = 4))
  expect_error(sar_average(data = galap, homoTest = "cor.fitted",
                           homoCor = "correlation"))
})

test_that("sar_average using fit_collection object works", {
  skip_on_cran()
  ff <- sar_multi(data = galap, obj = c("power", "p1", "loga", "monod",
                                        "linear"))
  expect_warning(sar_average(obj = ff, data = galap, normaTest = "none",
                    neg_check = FALSE, homoTest = "cor.fitted"))
  fit3 <- sar_average(obj = ff, data = galap)#grid start on so round to 0
  expect_equal(round(sum(fit3$mmi), 0), 1663)
  expect_output(str(fit3), "List of 2")
  expect_is(fit3, "multi")
  expect_match(fit3$details$homo_test, "none")
  expect_match(fit3$details$ic, "AICc")
  expect_error(sar_multi(5), "argument is of length zero")
})

test_that("confidence intervals are correct", {
  skip_on_cran()
  fit3 <- expect_warning(sar_average(data = galap, grid_start = "none",
                                     normaTest = "none", homoTest = "none",
                    neg_check = FALSE, confInt = TRUE, ciN = 20))
  ci <- fit3$details$confInt
  expect_equal(nrow(ci), nrow(galap))
  expect_output(str(ci), "data.frame")
  expect_true(all(ci[ ,1] < ci[ ,2]))
})



test_that("confidence intervals works with AICc", {
  skip_on_cran()
  fit3 <- expect_warning(sar_average(data = galap, grid_start = "none",
                                     normaTest = "none", homoTest = "none",
                      neg_check = FALSE, confInt = TRUE, crit = "AICc", ciN = 20))
  ci <- fit3$details$confInt
  expect_equal(nrow(ci), nrow(galap))
  expect_output(str(ci), "data.frame")
  expect_true(all(ci[ ,1] < ci[ ,2]))
})
