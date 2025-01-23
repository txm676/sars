context("sar_grid_start")

test_that("grid_start is working correctly", {
  skip_on_cran()
  data(aegean)
  #individual models
  expect_error(sar_power(aegean, grid_start = 1000), 
               "grid_start should be one of none, partial or exhaustive")
  expect_error(sar_power(aegean, grid_start = TRUE), 
               "grid_start should be one of none, partial or exhaustive")
  expect_error(sar_power(aegean, grid_start = "exhaustive", grid_n = NULL),
               "grid_n should be numeric if grid_start == exhaustive")
  expect_warning(sar_chapman(aegean, grid_start = "none"))#fails without GS.
  #there is a random component to grid_start so can't check par values and 
  #AIC etc. Instead check the length of the output
  expect_length(sar_chapman(aegean), 23)#same with grid_start, now works
  #check exhaustive works
  expect_length(sar_power(aegean, grid_start = "exhaustive", grid_n = 200), 23)
  #checking it does not error within sar_average
  df <- aegean[1:10,]
  obj <- c("linear","power","powerR",
           "chapman","gompertz")
  expect_error(sar_average(data = df, obj, grid_start = TRUE), 
               "grid_start should be one of 'none', 'partial' or 'exhaustive'")
  expect_error(sar_average(data = df, obj, grid_start = "exhaustive", grid_n = NULL), 
               "grid_n should be numeric if grid_start == exhaustive")
  #verb = FALSE as Gompertz & Chapman throw RSS warning
  expect_length(sar_average(data = df, obj, grid_start = "none", verb = FALSE), 2)
  expect_length(sar_average(data = df, obj, grid_start = "partial", verb = FALSE), 2)
  expect_length(sar_average(data = df, obj, 
                            grid_start = "exhaustive", grid_n = 200, verb = FALSE), 2)
  
})
