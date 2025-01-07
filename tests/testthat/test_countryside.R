
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
  expect_no_error(capture_output_lines(s))
  expect_equal(length(s), 8)
  expect_equal(class(s), c("habitat", "sars","list"))
  expect_equal(attributes(s)$modType, "power")

  expect_equal(as.vector(round(s$c, 0)),
               c(11, 4, 6, 1))
  
  expect_equal(as.vector(c(round(s$affinity$Spcs_AG[1], 1), 
    round(s$affinity$Spcs_SH[3], 7),
    round(s$affinity$Spcs_QF[2], 8),
    round(s$affinity$Spcs_UB[1], 2))),
    c(1, 2.12e-05, 2.3e-07, 0.42))
  
  expect_equal(as.vector(round(s$rss, 0)),
               c(5256, 10573))
  expect_equal(round(sum(s$fits$Spcs_AG$m$resid()^2),0),
               2590)
  
  expect_no_error(plot(s, type = 1))
  expect_no_error(plot(s, type = 2,
                  lcol = c("black", "aquamarine4",
                  "#CC661AB3", "darkblue"), 
                  pLeg = TRUE, lwd = 1.5, 
                  legPos = "topright"))
  expect_warning(plot(s, type = 2,
                       lcol = c("black", "aquamarine4"), 
                       pLeg = TRUE, lwd = 1.5, 
                       legPos = "topright"))
  
  #Check countryside_extrap
  expect_error(countryside_extrap(s, area = 1:5))
  b <- countryside_extrap(s, area = 1:3)
  expect_equal(round(b$Total,2), 23.64)
  expect_false(b$Failed_mods)
  
  #check plots
  expect_no_error(plot(s))
  expect_no_error(plot(s, type = 1))
  expect_no_error(plot(s, type = 2))
  expect_error(plot(s, type = 3))
  
  #Check provision of starting pars
  M2 <- matrix(c(3.061e+08, 2.105e-01, 1.075e+00, 1.224e-01,
  3.354e-08, 5.770e+05, 1.225e+01, 1.090e-01,
  6.848e-01, 1.054e-01, 4.628e+05, 1.378e-01,
  0.20747, 0.05259, 0.49393, 0.18725), nrow = 4,
  byrow = TRUE)
  s4 <- sar_countryside(data = countryside,
                      modType = "power",
                     startPar = M2, ubiSp = TRUE)
  expect_equal(as.vector(c(round(s4$affinity$Spcs_AG[1], 1), 
                           round(s4$affinity$Spcs_SH[3], 7),
                           round(s4$affinity$Spcs_QF[2], 8),
                           round(s4$affinity$Spcs_UB[1], 2))),
               c(1, 2.12e-05, 2.3e-07, 0.42))
  expect_equal(round(sum(s4$fits$Spcs_AG$m$resid()^2),0),
               2590)
})

test_that("sar_countryside logarithmic returns correct values", {
  skip_on_cran()
  data(countryside)
  s2 <- sar_countryside(data = countryside, 
                        ubiSp = TRUE, 
                        modType = "logarithmic")
  expect_equal(length(s2), 8)
  expect_equal(class(s2), c("habitat", "sars","list"))
  expect_equal(attributes(s2)$modType, "logarithmic")
  
  expect_equal(as.vector(round(s2$c, 0)),
               c(12, 5, 7, 1))
  
  expect_equal(as.vector(c(round(s2$affinity$Spcs_AG[1], 1), 
                           round(s2$affinity$Spcs_SH[3], 4),
                           round(s2$affinity$Spcs_QF[2], 4),
                           round(s2$affinity$Spcs_UB[1], 2))),
               c(1, 3e-04, 7e-04, 0.49))
  
  expect_equal(as.vector(round(s2$rss, 0)),
               c(8615, 15763))
  expect_equal(round(sum(s2$fits$Spcs_AG$m$resid()^2),0),
               5716)
  
  expect_no_error(plot(s2, type = 1))
  expect_no_error(plot(s2, type = 2,
                       lcol = c("black", "aquamarine4",
                                "#CC661AB3", "darkblue"), 
                       pLeg = TRUE, lwd = 1.5, 
                       legPos = "topright"))
  expect_warning(plot(s2, type = 2,
                      lcol = c("black", "aquamarine4"), 
                      pLeg = TRUE, lwd = 1.5, 
                      legPos = "topright"))
  
  #Check countryside_extrap
  expect_error(countryside_extrap(s2, area = 1:5))
  b <- countryside_extrap(s2, area = 1:3)
  expect_equal(round(b$Total,2), 26.19)
  expect_false(b$Failed_mods)
}) 