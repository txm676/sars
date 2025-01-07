context("sar_habitat")
library(sars)

test_that("sar_habitat errors where it should", {
  skip_on_cran()
  data(habitat)
  expect_error(sar_habitat(habitat, modType = "expo"),
               "modType should be one of 'power', 'logarithmic' or 'power_log'")
  expect_error(sar_habitat(habitat, logT = "expo"))
  habitat2 <- habitat
  habitat2$Species[1] <- 0
  expect_error(sar_habitat(habitat2, con = NULL),
               "The dataset has richness values of zero, con should be a numeric vector of length 1")
  expect_message(sar_habitat(habitat2, con = 1),
                 "The dataset has zero richness values, 1 has been added to all richness values.")
})

test_that("sar_habitat log_log returns correct values", {
  skip_on_cran()
  data(habitat)
  s <- sar_habitat(data = habitat, modType = "power_log", 
  con = NULL, logT = log)
  expect_equal(length(s), 4)
  expect_equal(class(s), c("habitat", "sars","list"))
  expect_equal(attributes(s)$modType, "power_log")
  s2 <- summary(s)
  expect_equal(s2$modType, "power_log")
  #jigsaw
  wjig <- which(s2$Model_table$Model == "jigsaw")
  expect_equal(round(c(s2$Model_table$AIC[wjig],
                 s2$Model_table$z[wjig],
                 s2$Model_table$`d-z`[wjig],
                 s2$Model_table$R2[wjig]), 3),
               c(-3.100, 0.086, 0.414, 0.940))
  #choros
  wcho <- which(s2$Model_table$Model == "choros")
  expect_equal(round(c(s2$Model_table$AICc[wcho],
                       s2$Model_table$z[wcho],
                       s2$Model_table$R2[wcho]), 3),
               c(2.228, 0.147, 0.914))
  #Kalli
  wkal <- which(s2$Model_table$Model == "Kallimanis")
  expect_equal(round(c(s2$Model_table$AICc[wkal],
                       s2$Model_table$z[wkal],
                       s2$Model_table$d[wkal]), 3),
               c(11.150, 0.179, -0.000))
  #power
  wpow <- which(s2$Model_table$Model == "power")
  expect_equal(round(c(s2$Model_table$AICc[wpow],
                       s2$Model_table$z[wpow]), 3),
               c(6.444, 0.177))
  
  #Check IC equations match up
  ll_cho <- logLik(s$choros)
  n <- length(s$choros$residuals)
  K <- attributes(ll_cho)$df
  AICcm <- -2*ll_cho+2*K*(n/(n-K-1))
  IC_cho <- round(c("AIC" = AIC(s$choros), 
                            "BIC" = BIC(s$choros), 
                    "AICc" = AICcm),3)
  expect_equal(unlist(round(s2$Model_table[wcho,
          c("AIC", "BIC", "AICc")], 3)),
          IC_cho)
  
  ##Check plots
  expect_no_error(plot(s))
  expect_no_error(plot(s, col = "red"))
  })

test_that("sar_habitat logarithmic returns correct values", {
  skip_on_cran()
  data(habitat)
  s3 <- sar_habitat(data = habitat, 
                    modType = "logarithmic", 
                   con = NULL, logT = log)
  expect_equal(length(s3), 4)
  expect_equal(class(s3), c("habitat", "sars","list"))
  expect_equal(attributes(s3)$modType, "logarithmic")
  s4 <- summary(s3)
  expect_equal(s4$modType, "logarithmic")
  
  #jigsaw
  wjig <- which(s4$Model_table$Model == "jigsaw")
  expect_equal(round(c(s4$Model_table$AIC[wjig],
                       s4$Model_table$z[wjig],
                       s4$Model_table$d[wjig],
                       s4$Model_table$R2[wjig]), 3),
               c(66.053, 0.615, 0.454, 0.923))
  
  #logarithmic model
  wlog <- which(s4$Model_table$Model == "logarithmic")
  s5 <- sar_loga(habitat[,c("Area", "Species")])
  expect_equal(round(s5$AIC, 3), 
               round(s4$Model_table$AIC[wlog], 3))
  expect_equal(round(as.vector(s5$par[2]), 3), 
               round(s4$Model_table$z[wlog], 3))
  
  #compare minpack.lm::nlsLM with nls
  s3_nls <- nls(Species ~ (Het^d) * 
                  log(c1 * (Area / Het)^z),
                 start = list("c1" = 5,
                              "z" = 1,
                              "d" = 0.6),
                 data = habitat)
  expect_equal(AIC(s3_nls), AIC(s3$jigsaw))
  expect_equal(round(sum(s3$jigsaw$m$resid()), 3), 
               round(sum(s3_nls$m$resid())),3)
  
  ##hashed out as requires external dataset
  # dat <- read.csv("Data - Table 1 - Reed 1981.csv")
  # s9 <- sar_habitat(data = dat, 
  #                   modType = "logarithmic", 
  #                   con = NULL, logT = log)
  # s9_nls <- nls(Species ~ (Het^d) * 
  #                 logT(c1 * (Area / Het)^z),
  #               start = list("c1" = 1.73525,
  #                            "z" = 0.08795,
  #                            "d" = 1.28740),
  #               data = dat)
  # expect_equal(AIC(s9_nls), AIC(s9$jigsaw))
  # expect_equal(round(sum(s9$jigsaw$m$resid()), 3), 
  #              round(sum(s9_nls$m$resid())),3)
  # expect_equal(s9_nls$convInfo$stopMessage, 
  #              "converged")
  # expect_true(s9$jigsaw$convInfo$isConv)
  ##Check plots
  expect_no_error(plot(s3))
  expect_no_error(plot(s3, col = "red"))
}) 

test_that("sar_habitat untransformed returns correct values", {
  skip_on_cran()
  data(habitat)
  s6 <- sar_habitat(data = habitat, 
                    modType = "power", 
                    con = NULL, logT = log)
  expect_no_error(plot(s6))
  expect_equal(length(s6), 4)
  expect_equal(class(s6), c("habitat", "sars","list"))
  expect_equal(attributes(s6)$modType, "power")
  s7 <- summary(s6)
  expect_equal(s7$modType, "power")
  
  #kallimanis
  wkal <- which(s7$Model_table$Model == "Kallimanis")
  expect_equal(round(c(s7$Model_table$AICc[wkal],
                       s7$Model_table$z[wkal],
                       s7$Model_table$d[wkal],
                       s7$Model_table$R2[wkal]), 3),
               c(59.758, 0.172, 0.000, 0.972))
  
  #power model
  wpow <- which(s7$Model_table$Model == "power")
  s8 <- sar_power(habitat[,c("Area", "Species")])
  expect_equal(round(s8$AIC, 3), 
               round(s7$Model_table$AIC[ wpow], 3))
  expect_equal(round(as.vector(s8$par[2]), 3), 
               round(s7$Model_table$z[ wpow], 3))
}) 
