
context("sar_countryside")
library(sars)

test_that("sar_countryside errors where it should", {
  skip_on_cran()
  data(countryside)
  expect_error(sar_countryside(countryside, modType = "expo"),
               "modType should be one of power or logarithmic")
  expect_error(sar_countryside(countryside, zLower = 1:4))
  c2 <- countryside[,1:6]
  expect_error(sar_countryside(c2, ubiSp = TRUE))
  expect_error(sar_countryside(countryside, habNam = letters[1:5]))
  c2[1,1:3] <- 0
  expect_error(sar_countryside(c2, modType = "logarithmic"))
  c2 <- c2[1:10,]
  expect_warning(sar_countryside(c2, modType = "power"),
                 "Some sites have total area equal to zero")
})

test_that("sar_countryside power returns correct values", {
  skip_on_cran()
  data(countryside)
  expect_equal(colnames(countryside), 
               c("Area_AG", "Area_SH", "Area_QF", "Spcs_AG",
                "Spcs_SH", "Spcs_QF", "Spcs_UB"))
  s <- sar_countryside(data = countryside, ubiSp = TRUE)
  expect_equal(length(s), 8)
  expect_equal(class(s), c("habitat", "sars","list"))
  expect_equal(attributes(s)$modType, "power")

  expect_equal(as.vector(round(s$c, 0)),
               c(11, 4, 6, 1))
  
  
  c(round(s$affinity$Spcs_AG[1], 1), 
    round(s$affinity$Spcs_SH[3], 7),
    round(s$affinity$Spcs_QF[2], 7),
    round(s$affinity$Spcs_UB[1], 2))
  
  c(1, 2.12e-05, )
  
  
  
  
  
  
})

test_that("sar_countryside logarithmic returns correct values", {
  skip_on_cran()
  data(countryside)
  s3 <- sar_countryside(data = countryside, 
                    modType = "logarithmic", 
                    con = NULL, logT = log)
  expect_equal(length(s3), 4)
  expect_equal(class(s3), c("countryside", "sars","list"))
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
  s5 <- sar_loga(countryside[,c("Area", "Species")])
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
                data = countryside)
  expect_equal(AIC(s3_nls), AIC(s3$jigsaw))
  expect_equal(round(sum(s3$jigsaw$m$resid()), 3), 
               round(sum(s3_nls$m$resid())),3)
  
  ##hashed out as requires external dataset
  # dat <- read.csv("Data - Table 1 - Reed 1981.csv")
  # s9 <- sar_countryside(data = dat, 
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
}) 

test_that("sar_countryside untransformed returns correct values", {
  skip_on_cran()
  data(countryside)
  s6 <- sar_countryside(data = countryside, 
                    modType = "power", 
                    con = NULL, logT = log)
  expect_equal(length(s6), 4)
  expect_equal(class(s6), c("countryside", "sars","list"))
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
  s8 <- sar_power(countryside[,c("Area", "Species")])
  expect_equal(round(s8$AIC, 3), 
               round(s7$Model_table$AIC[ wpow], 3))
  expect_equal(round(as.vector(s8$par[2]), 3), 
               round(s7$Model_table$z[ wpow], 3))
}) 