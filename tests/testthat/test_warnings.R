context("warnings")


test_that("verb warnings are correctly returned or not", {
  skip_on_cran()
  #this dataset is from the ant paper and the epm2 model does not fit
  #properly (just flat line) and should return the same min RSS warning
  ant <- data.frame("A" = c(3.210e-06, 3.740e-06, 4.050e-06, 4.210e-06,
                            4.990e-06, 5.270e-06, 5.560e-06, 5.770e-06,
                            6.910e-06, 7.970e-06, 8.000e-06, 8.560e-06,
                            3.013e-05, 3.921e-05, 4.616e-05, 6.053e-05,
                            9.337e-05),
                    "S" = c(0, 13, 13, 14, 16, 16, 17, 16, 14, 20, 20, 23,
                            38, 47, 58, 59, 77))
  expect_warning(sar_epm2(ant), paste0("Extended Power model 2: Multiple parameter", 
                                       " estimates returned the same minimum rss; one set have been", 
                                       " randomly selected"))
  #turn verb off and set regexp = NA - to check for no warning
  expect_warning(sar_epm2(ant, verb = FALSE),  regexp = NA)
  expect_error(sar_epm2(ant, verb = "T"))
  
  ##without residual checks
  mods <- c("power", "linear", "monod", "weibull3", "p1", "epm2")
  #checking verb
  expect_warning(sar_average(data = ant, obj = mods, grid_start = "partial", 
              verb = TRUE, display = FALSE),
              paste0("Extended Power model 2: Multiple parameter", 
                     " estimates returned the same minimum rss; one set have been", 
                     " randomly selected"))
  expect_warning(sar_average(data = ant, obj = mods, grid_start = "partial", 
                             verb = FALSE, display = FALSE),  regexp = NA)
  ##with residual checks
  #checking verb
  expect_warning(sar_average(data = ant, obj = mods, grid_start = "partial", 
                             verb = TRUE, display = TRUE, normaTest = "lillie"),
                 paste0("Extended Power model 2: Multiple parameter", 
                        " estimates returned the same minimum rss; one set have been", 
                        " randomly selected"))
  expect_message(sar_average(data = galap, obj = mods, grid_start = "none", 
                             verb = TRUE, display = TRUE, normaTest = "lillie"),
                 "1 models failed the residuals normality test")
  expect_message(sar_average(data = galap, obj = mods, grid_start = "none", 
                             verb = TRUE, display = FALSE, homoTest = "cor.area"),
                 "3 models failed the residuals homogeneity")
  expect_message(sar_average(data = galap, obj = mods, grid_start = "none", 
                             verb = FALSE, display = TRUE, normaTest = "lillie"),
                 regexp = NA)
  expect_error(sar_average(data = galap, verb = 3))
  expect_error(sar_average(data = galap, display = 3))
  expect_error(sar_multi(data = galap, verb = 3))
  expect_error(sar_multi(data = galap, display = "TRUE"))
})