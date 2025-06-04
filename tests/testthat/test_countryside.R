
context("sar_countryside")
library(sars)

test_that("sar_countryside errors where it should", {
  skip_on_cran()
  data(countryside)
  expect_error(sar_countryside(countryside, modType = "expo"),
               "modType should be one of power or logarithmic")
  expect_error(sar_countryside(countryside, zLower = 1:4))
  c2 <- countryside[,1:6]
  expect_error(sar_countryside(c2, habNam = NULL))
  expect_error(sar_countryside(c2, habNam = 1:3,spNam = NULL))
  expect_error(sar_countryside(c2, habNam = 1:3, spNam = TRUE))
  expect_error(sar_countryside(c2,
                               habNam = letters[1:5], spNam = 6:7))
  expect_error(sar_countryside(countryside,
                               habNam = letters[1:5], spNam = 1:2))
  c2[1,1:3] <- 0
  expect_error(sar_countryside(c2, modType = "logarithmic"))
  c2 <- c2[1:10,]
  expect_warning(sar_countryside(c2, modType = "power",
                    habNam = letters[1:3], spNam = 4:6),
                 "Some sites have total area equal to zero")
})

test_that("sar_countryside power returns correct values", {
  skip_on_cran()
  data(countryside)
  expect_equal(colnames(countryside),
               c("Area_AG", "Area_SH", "Area_QF", "Spcs_AG",
                "Spcs_SH", "Spcs_QF", "Spcs_UB"))
  s <- sar_countryside(data = countryside,
                       habNam = 1:3, spNam = 4:7)
  expect_equal(length(capture_output_lines(s, print = TRUE)),
               91)
  expect_equal(length(s), 8)
  expect_equal(class(s), c("habitat", "sars","list"))
  expect_equal(attributes(s)$modType, "power")

  expect_equal(as.vector(round(s$c, 0)),
               c(11, 4, 6, 1))

  expect_equal(as.vector(c(round(s$affinity$Sp_grp1[1], 1),
    round(s$affinity$Sp_grp2[3], 7),
    round(s$affinity$Sp_grp3[2], 8),
    round(s$affinity$Sp_grp4[1], 2))),
    c(1, 2.12e-05, 2.3e-07, 0.42))

  expect_equal(as.vector(round(s$rss, 0)),
               c(5256, 10573))
  expect_equal(round(sum(s$fits$Sp_grp1$m$resid()^2),0),
               2590)

  expect_equal(round(sum(s$fits$Sp_grp3$m$resid()^2),
                    0),1084)

  #Calculate AICc using Proenca approach (our old approach)
  nc1 <- nrow(countryside)
  kc1 <- 5
  RSSc <- sum(s$fits$Sp_grp3$m$resid()^2)
  AICC <- (nc1 * log(RSSc/nc1)) + (2*kc1)*(nc1 / (nc1 - kc1 - 1))
  expect_equal(round(AICC,0), 408)
  aa <- rowSums(countryside[,1:3])
  df <- data.frame(aa, countryside[,6])
  sss <- sar_power(df)
  aiccP <- (nc1 * log(sss$value/nc1)) +  (2*3)*(nc1 / (nc1 - 3 - 1))
  expect_equal(round(aiccP,0), 1156)

  expect_no_error(plot(s, type = 1))
  #Circle CI doesn't work with the Enter Return plots, so
  #have to use which
  # expect_no_error(plot(s, type = 2,
  #                 lcol = c("black", "aquamarine4",
  #                 "#CC661AB3", "darkblue"),
  #                 pLeg = TRUE, lwd = 1.5,
  #                 legPos = "topright"))
  expect_no_error(plot(s, type = 2,
                       lcol = c("black", "aquamarine4",
                                "#CC661AB3", "darkblue"),
                       pLeg = TRUE, lwd = 1.5,
                       legPos = "topright",
                       which = 2, ModTitle = "S"))
  expect_no_error(plot(s, type = 2,
                       lcol = c("black", "aquamarine4",
                                "#CC661AB3", "darkblue"),
                       pLeg = TRUE, lwd = 1.5,
                       legPos = "topright",
                       which = 2, ModTitle = letters[1:3]))
  expect_warning(plot(s, type = 2,
                       lcol = c("black", "aquamarine4"),
                       pLeg = TRUE, lwd = 1.5,
                       legPos = "topright", which = 1))
  #Circle CI doesn't work with the Enter Return plots, so
  #have to use which
  # expect_no_error(plot(s, type = 3,
  #                      lcol = c("black", "aquamarine4",
  #                               "#CC661AB3", "darkblue"),
  #                      pLeg = TRUE, lwd = 1.5,
  #                      legPos = "topright"))
  expect_no_error(plot(s, type = 3,
                       lcol = c("black", "aquamarine4",
                                "#CC661AB3", "darkblue"),
                       pLeg = TRUE, lwd = 1.5,
                       legPos = "topright",
                       which = 3, ModTitle = "S"))
  
  expect_error(plot(s, type = 4, which = 1:2))
  expect_error(plot(s, type = 5))
  
  expect_no_error(plot(s, type = 4, which = 1))
  
  expect_no_error(plot(s, type = 3, which = 1, 
                       ModTitle = c("Agricultural land")))
  
  expect_no_error(plot(s, type = 3, which = 1))
  
  expect_no_error(plot(s, type = 3, which = 1, 
                       ModTitle = "none"))
  expect_error(plot(s, type = 3, which = 1, 
                       ModTitle = c("R", "E")))
  expect_error(plot(s, type = 4, 
                    ModTitle = c("R", "E")))
  #Check countryside_extrap
  expect_error(countryside_extrap(s, area = 1:5))
  b <- countryside_extrap(s, area = 1:3)
  expect_equal(round(b$Total,2), 23.64)
  expect_false(b$Failed_mods)

  #Check provision of starting pars
  M2 <- matrix(c(3.061e+08, 2.105e-01, 1.075e+00, 1.224e-01,
  3.354e-08, 5.770e+05, 1.225e+01, 1.090e-01,
  6.848e-01, 1.054e-01, 4.628e+05, 1.378e-01,
  0.20747, 0.05259, 0.49393, 0.18725), nrow = 4,
  byrow = TRUE)
  s4 <- sar_countryside(data = countryside,
                      modType = "power",
                     startPar = M2,
                     habNam = 1:3, spNam = 4:7)
  expect_equal(as.vector(c(round(s4$affinity$Sp_grp1[1], 1),
                           round(s4$affinity$Sp_grp2[3], 7),
                           round(s4$affinity$Sp_grp3[2], 8),
                           round(s4$affinity$Sp_grp4[1], 2))),
               c(1, 2.12e-05, 2.3e-07, 0.42))
  expect_equal(round(sum(s4$fits$Sp_grp1$m$resid()^2),0),
               2590)

  ##Test output still works if you mix up columns and remove
  #a species group
  Ff <- countryside[,c(1,3,2,7,4,6,5)]
  s5 <- sar_countryside(data = Ff, modType = "power",
                        gridStart = "none",
                        habNam = c("AG", "F", "SH"),
                        spNam = c("UB_Sp","AG_Sp",  "F_Sp",
                                  "SH_Sp"))
  expect_equal(as.vector(round(s5$c, 0)),
               c(1, 11, 6, 4))

  expect_equal(as.vector(c(round(s5$affinity$AG_Sp[1], 1),
                           round(s5$affinity$SH_Sp[2], 7),
                           round(s5$affinity$F_Sp[3], 8),
                           round(s5$affinity$UB_Sp[1], 2))),
               c(1, 2.12e-05, 2.3e-07, 0.42))

  expect_equal(as.vector(round(s5$rss, 0)),
               c(5256, 10573))
  expect_equal(round(sum(s5$fits$AG_Sp$m$resid()^2),0),
               2590)
  ##To try solve weird fails I was getting on Circle-CI
  Ff <- NULL
  s5 <- NULL

 ##Test function works with fewer richness cols than area
  c4 <- countryside[1:50,c(1,2,3,4,5)]
  s6 <- sar_countryside(data = c4, modType = "power",
                        gridStart = "partial",
                        habNam = c("AG", "SH", "F"),
                        spNam = c("AG_Sp", "SH_Sp"))
  expect_equal(length(s6), 8)
  expect_equal(names(s6$fits), c("AG_Sp", "SH_Sp"))

  ##Test using gridStart = "none"
  Jj2 <- sar_countryside(data = countryside,
                         modType = "power",
                        gridStart = "none",
                        habNam = c("AG", "F", "SH"),
                        spNam = c("SH_Sp","AG_Sp",  "F_Sp",
                                  "UB_Sp"))
  expect_equal(as.vector(round(Jj2$rss, 0)),
               c(5256, 10573))
  #check works with tibble
  cp2 <- tibble::as_tibble(countryside)
  expect_no_error(sar_countryside(data = cp2,
                                  modType = "power",
                                  gridStart = "none",
                                  habNam = c("AG", "F", "SH"),
                                  spNam = c("SH_Sp","AG_Sp",  "F_Sp",
                                            "UB_Sp")))
})

 test_that("sar_countryside logarithmic returns correct values", {
  skip_on_cran()
  data(countryside)
  s2 <- sar_countryside(data = countryside,
                        modType = "logarithmic",
                        habNam =  c("Area_AG", "Area_SH",
                                    "Area_QF"),
                        spNam = c("Spcs_AG",
                                  "Spcs_SH", "Spcs_QF",
                                  "Spcs_UB"))
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
                       legPos = "topright", which = 1))
  expect_no_error(plot(s2, type = 3,
                       lcol = c("black", "aquamarine4",
                                "#CC661AB3", "darkblue"),
                       pLeg = TRUE, lwd = 1.5,
                       legPos = "topright", which = 2))
  expect_warning(plot(s2, type = 2,
                      lcol = c("black", "aquamarine4"),
                      pLeg = TRUE, lwd = 1.5,
                      legPos = "topright", which = 3))
  expect_no_error(plot(s2, type = 4, which = 2))
  expect_error(plot(s, type = 4, which = 5))

  #Check countryside_extrap
  expect_error(countryside_extrap(s2, area = 1:5))
  b <- countryside_extrap(s2, area = 1:3)
  expect_equal(round(b$Total,2), 26.19)
  expect_false(b$Failed_mods)
})

# #Tested on second dataset: hashed out for speed
# test_that("sar_countryside power works with 2nd dataset", {
#   skip_on_cran()
# 
#   #In Henrique_tests drive
#   guillerme <- read.csv("guilherme.csv")
# 
#   #This version works with gridStart = "none"
#   sg <- sar_countryside(data = guillerme,
#                        habNam = c("AG", "SH","F"),
#                        spNam = c("AG_Sp", "SH_Sp",
#                                  "F_Sp", "UB_Sp"),
#                        gridStart = "none")
#   expect_equal(length(capture_output_lines(sg, print = TRUE)),
#                90)
#   expect_equal(length(sg), 8)
#   expect_equal(class(sg), c("habitat", "sars","list"))
#   expect_equal(attributes(sg)$modType, "power")
# 
#   expect_equal(as.vector(round(sg$c, 1)),
#                c(0.1, 0.4, 0.6, 0.7))
# 
#   expect_equal(as.vector(c(round(sg$affinity$AG_Sp[1], 1),
#                            round(sg$affinity$SH_Sp[3], 7),
#                            round(sg$affinity$F_Sp[2], 8),
#                            round(sg$affinity$UB_Sp[1], 2))),
#                c(1, 0.001008, 7.85e-06 , 0.91 ))
# 
#   # #Calculate AICc using Proenca approach (our old approach)
#  nc1 <- 125 #HP confirmed this was the actual N used
#  kc1 <- 5
#  RSSc <- sum(sg$fits$AG_Sp$m$resid()^2)
#  AICC <- (nc1 * log(RSSc/nc1)) +
#    (2*kc1)*(nc1 / (nc1 - kc1 - 1))
#  expect_equal(round(AICC,0), -8)
# 
#   expect_no_error(plot(sg, type = 1))
# 
#   expect_no_error(plot(sg, type = 2,
#                        lcol = c("black", "aquamarine4",
#                                 "#CC661AB3", "darkblue"),
#                        pLeg = TRUE, lwd = 1.5,
#                        legPos = "topright",
#                        which = 2, ModTitle = letters[1:3]))
#   expect_warning(plot(sg, type = 2,
#                       lcol = c("black", "aquamarine4"),
#                       pLeg = TRUE, lwd = 1.5,
#                       legPos = "topright", which = 1))
# 
#   expect_no_error(plot(sg, type = 3,
#                        lcol = c("black", "aquamarine4",
#                                 "#CC661AB3", "darkblue"),
#                        pLeg = TRUE, lwd = 1.5,
#                        legPos = "topright",
#                        which = 3,
#                        ModTitle = "F"))
# })
