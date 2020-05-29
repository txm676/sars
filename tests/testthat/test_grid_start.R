context("sar_grid_start")

test_that("grid_start is working correctly", {
  data(aegean)
  expect_error(sar_power(aegean, grid_search = 1000))
  expect_error(sar_power(aegean, grid_search = TRUE, grid_n = FALSE))
  expect_warning(sar_chapman(aegean))
  s <- suppressWarnings(sar_chapman(aegean))
  expect_equal(0.316, as.vector(round(s$par[2], 3)))
  s2 <- sar_chapman(aegean, grid_start = TRUE, grid_n = 100)
  #there is a random component to grid_start so you can't check whether
  #parameter values are equal; occasionally there might be small differences.
  #expect_equal(0.044, as.vector(round(s2$par[2], 3)))
  #so instead, just check it returns a numeric parameter
  expect_true(is.numeric(s2$par[2]))
  #checking it does not error within sar_average
  df <- aegean[1:10,]
  obj <- c("linear","power","powerR",
           "chapman","gompertz")
  expect_error(sar_average(data = df, obj, grid_start = TRUE, grid_n = 5), NA)
})
