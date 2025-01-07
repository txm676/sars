context("gdm")
library(sars)

test_that("gdm functions return correct results", {
  data(galap)
  galap$t <- c(4, 1, 13, 16, 15, 2, 6, 4, 5, 11, 3, 9, 8, 10, 12, 7)
  data(aegean)
  
  ###matches with sar_power_
  g2 <- sars::gdm(galap, 
                  model = "power_area", 
                  mod_sel = TRUE)
  expect_equal(length(capture_output_lines(g2, print = TRUE)),
               22)
  pow <- g2[[3]]
  pow_sars <- sar_power(galap)
  
  expect_equal(as.vector(round(pow_sars$par[2], 2)), 
               as.vector(round(pow$m$getPars()[2], 2)))
  expect_equal(as.vector(round(pow_sars$par[1], 1)), 
               as.vector(round(exp(pow$m$getPars()[1]), 1)))
  expect_equal(pow_sars$AIC, AIC(pow))
  expect_error(sars::gdm(aegean), "Not enough columns/variables to fit GDM")
  expect_error(sars::gdm(galap, start_vals = c(5,4,3,6)), 
               "start_vals should be a dataframe with 1 row and four columns")
  expect_error(sars::gdm(galap, model = "betap"), "provided model name not available")
  expect_error(sars::gdm(galap, mod_sel = "d"), 
               "mod_sel argument should be TRUE or FALSE")
  #matches with BAT values
 # BTP <- BAT::gdm(galap$s, area = galap$a, time = galap$t)
 # expect_equal(as.vector(round(BTP[3,3], 2)),
 #             as.vector(round(g2[[1]]$m$getAllPars()[3], 2)))
 # expect_equal(as.vector(round(BTP[3,4], 2)),
 #             as.vector(round(g2[[1]]$m$getAllPars()[4], 2)))
 # expect_equal(as.vector(round(BTP[3,2], 2)),
 #              as.vector(round(g2[[1]]$m$getAllPars()[2], 2)))
 #  
  #power_area_time matches
 g2b <- sars::gdm(galap, model = "power_area_time", mod_sel = TRUE)
 pow <- g2b[[3]]
 pow_sars <- sar_power(galap)
 
 expect_equal(as.vector(round(pow_sars$par[2], 2)), 
              as.vector(round(pow$m$getPars()[2], 2)))
 expect_equal(as.vector(round(pow_sars$par[1], 1)), 
              as.vector(round(exp(pow$m$getPars()[1]), 1)))
 expect_equal(pow_sars$AIC, AIC(pow))
 #matches with BAT values
 # expect_equal(as.vector(round(BTP[4,3], 2)),
 #              as.vector(round(g2b[[1]]$m$getAllPars()[3], 2)))
 # expect_equal(as.vector(round(BTP[4,4], 2)),
 #              as.vector(round(g2b[[1]]$m$getAllPars()[4], 2)))
 # expect_equal(as.vector(round(BTP[4,2], 2)),
 #              as.vector(round(g2b[[1]]$m$getAllPars()[2], 2)))
 # expect_equal(round(as.vector(round(BTP[4,5], 2)), 2), #R2
 #              as.vector(round(sars:::R2_gdm(g2b[[1]]), 2)))
 
 ###matches with sar_loga
  g3 <- sars::gdm(galap, model = "loga", mod_sel = TRUE)
  loga <- g3[[3]]
  loga_sars <- sar_loga(galap)
  
  expect_equal(as.vector(round(loga_sars$par[2], 2)), 
               as.vector(round(loga$m$getPars()[2], 2)))
  expect_equal(as.vector(round(loga_sars$par[1], 1)), 
               as.vector(round(loga$m$getPars()[1], 1)))
  expect_equal(loga_sars$AIC, AIC(loga))
  #matches with BAT values
# expect_equal(as.vector(round(BTP[2,3], 2)),
#             as.vector(round(g3[[1]]$m$getAllPars()[3], 2)))
# expect_equal(as.vector(round(BTP[2,4], 2)),
#              as.vector(round(g3[[1]]$m$getAllPars()[4], 2)))
#  expect_equal(as.vector(round(BTP[2,2], 2)),
#               as.vector(round(g3[[1]]$m$getAllPars()[2], 2)))
#  expect_equal(round(as.vector(round(BTP[2,5], 2)), 2), #R2
#               as.vector(round(sars:::R2_gdm(g3[[1]]), 2)))
 
  ###matches with sar_linear
  g4 <- sars::gdm(galap, model = "linear", mod_sel = TRUE)
  linear <- g4[[3]]
  lin_sars <- sar_linear(galap)
  
  expect_equal(as.vector(round(lin_sars$par[2], 2)), 
               as.vector(round(linear$m$getPars()[2], 2)))
  expect_equal(as.vector(round(lin_sars$par[1], 1)), 
               as.vector(round(linear$m$getPars()[1], 1)))
  expect_equal(lin_sars$AIC, AIC(linear))
  #matches with BAT values
# expect_equal(as.vector(round(BTP[1,3], 2)),
#              as.vector(round(g4[[1]]$m$getAllPars()[3], 2)))
#  expect_equal(as.vector(round(BTP[1,4], 2)),
#               as.vector(round(g4[[1]]$m$getAllPars()[4], 2)))
#  expect_equal(as.vector(round(BTP[1,2], 2)),
#               as.vector(round(g4[[1]]$m$getAllPars()[2], 2)))
  
  ###all model comparison matches with individual model fits
  g5 <- sars::gdm(galap, model = "all", mod_sel = FALSE)
  expect_equal(length(capture_output_lines(g5, print = TRUE)),
               8)
  expect_equal(AIC(g4[[1]]), AIC(g5[[2]]))
  expect_equal(AIC(g2[[1]]), AIC(g5[[3]]))
  expect_equal(AICcmodavg::AICc(g3[[1]]), AICcmodavg::AICc(g5[[1]]))
  expect_equal(AICcmodavg::AICc(g2b[[1]]), AICcmodavg::AICc(g5[[4]]))
  #matches with BAT values (Delta AIC of linear (and _area_time) vs. power)
 #  expect_equal(round(AIC(g5[[2]]) - AIC(g5[[3]]), 2), round(BTP[1,7], 2))
 # expect_equal(round(AICcmodavg::AICc(g5[[2]]) - AICcmodavg::AICc(g5[[3]]), 2),
 #             round(BTP[1,9], 2))
 # expect_equal(round(AIC(g5[[4]]) - AIC(g5[[3]]), 2), round(BTP[4,7], 2))
 # expect_equal(round(AICcmodavg::AICc(g5[[4]]) - AICcmodavg::AICc(g5[[3]]), 2),
 #              round(BTP[4,9], 2))

  g5b <- sars::gdm(galap, model = "all", mod_sel = TRUE)
  expect_equal(length(g5b), 4)
  expect_equal(length(g5b[[4]]), 4)
  expect_equal(length(g5b[[2]]), 4)
  
  ###linear power version
  g6 <- sars::gdm(galap, model = "ATT2", mod_sel = TRUE)
  expect_equal(length(capture_output_lines(g6, print = TRUE)),
               19)
  linP <- g6[[3]]
  lin_sars <- sar_loga(galap)
  
  expect_equal(as.vector(round(linP$coefficients[2], 2)), 
               as.vector(round(lin_sars$par[2], 2)))
  expect_equal(as.vector(round(linP$coefficients[1], 2)), 
               as.vector(round(lin_sars$par[1], 2)))
  
  g6b <- sars::gdm(galap, model = "ATT2", mod_sel = FALSE)
  expect_equal(class(g6b), c("gdm", "lm"))
  
  expect_error(sars::gdm(galap, model = "lin_pow", mod_sel = FALSE))
  
  g7 <- sars::gdm(galap, model = "ATT2", mod_sel = FALSE)
  expect_equal(round(g7$coefficients, 2), round(g6[[1]]$coefficients, 2))
  expect_equal(as.vector(round(g7$coefficients, 2)),
               c(-104.00, 14.35, 49.97, -2.84))
  
})
