
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
  expect_equal(s2$Model_table$AIC, c(1974.90, 1974.87, 1977.10, 2195.69,
                                     2426.60))
  expect_equal(length(s2$Model_table$Th2[!is.na(s2$Model_table$Th2)]), 3)
  
  #table 1 in paper
  fit4 <- sar_threshold(data = aegean2, mod = "All", interval = 0.1, 
                        non_th_models = TRUE, logAxes = "area", 
                        logT = log10, parallel = TRUE, cores = 2)
  
  s4 <- summary(fit4, order = "BIC") 
  expect_equal(c(s4$Model_table$AIC), c(2020.25, 2019.48, 2022.42, 
                                        2045.94, 2047.78, 2061.70,
                                        2305.05, 2535.84))
  #repeat but without parallel processing
  fit5 <- sar_threshold(data = aegean2, mod = "All", interval = 0.1, 
                        non_th_models = TRUE, logAxes = "area", 
                        logT = log10, parallel = F)
  
  s5 <- summary(fit5, order = "BIC") 
  expect_equal(c(s5$Model_table$AIC), c(2020.25, 2019.48, 2022.42, 
                                        2045.94, 2047.78, 2061.70,
                                        2305.05, 2535.84))
  
  #check AIC and BIC returns same as normal AIC/BIC functions
  data(aegean)
  xd2 <- sar_threshold(aegean, mod = c("DiscTwo"), non_th_models = F)
  xs <- summary(xd2)
  obj <- xd2[[1]][[1]]
  n <- length(obj$residuals)
  val <- logLik(obj)
  expect_equal(attributes(val)$df, 7)#lm doesn't account for two threshold parameters
  attr(val, 'df') <- 9 #so need to add these two in for this model
  P <- 9
  lAIC2 <- (2 * P) - (2 * val)
  lBIC2 <- (-2 * val) + (P * log(n))
  expect_equal(as.vector(lAIC2), AIC(val))
  expect_equal(round(xs$Model_table$AIC, 2), round(AIC(val), 2))
  expect_equal(as.vector(lBIC2), BIC(val))
  expect_equal(round(xs$Model_table$BIC, 2), round(BIC(val), 2))
  #AICc 
  LL <- xs$Model_table$LL
  K <- 9
  n <- nrow(aegean)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  expect_equal(round(xs$Model_table$AICc, 1), round(AICcm, 1))
})


test_that("nisl argument returns correct results", {
  data(aegean2)
  a2 <- aegean2[1:168,]
  fitT <- sar_threshold(data = a2, mod = "All",
                        interval = 0.1, nisl = 75, non_th_models = TRUE, 
                        logAxes = "both", logT = log2)
  
  s <- summary(fitT)
  expect_false(any(na.omit(s$Model_table$seg1) < 75))
  expect_false(any(na.omit(s$Model_table$seg3) < 75))
  expect_error(sar_threshold(data = a2, mod = "All",
                             interval = 0.1, nisl = 750), 
          "nisl is equal or larger than half of the total number of islands")
  #check with exactly half, i.e. nisl half the number of islands (Should error)
  expect_error(sar_threshold(data = a2, mod = "All",
                             interval = 0.1, nisl = 84), 
       "nisl is equal or larger than half of the total number of islands")
})


test_that("threshold_ci returns correct results", {
  data(aegean2)
  a2 <- aegean2[1:168,]
  #fit with just one model (boot and F methods)
  fitT <- sar_threshold(data = a2, mod = "ContOne", interval = 0.1, 
                        non_th_models = TRUE, logAxes = "area", logT = log10)
  CI1 <- threshold_ci(fitT, method = "boot", interval = NULL, Nboot = 3)
  
  fitT <- sar_threshold(data = a2, mod = "ContOne", interval = 0.1, 
                        non_th_models = TRUE, logAxes = "area", logT = log10)
  CI2 <- threshold_ci(fitT, method = "F", interval = NULL, Nboot = 3)
  
  #fit with both models (boot and F methods)
  fitT <- sar_threshold(data = a2, mod = c("ContOne","ZslopeOne"),
                        interval = 0.1, non_th_models = TRUE, 
                        logAxes = "area", logT = log10)
  CI3 <- threshold_ci(fitT, method = "boot", interval = NULL, Nboot = 3)
  
  fitT <- sar_threshold(data = a2, mod = c("ContOne","ZslopeOne"),
                        interval = 0.1, non_th_models = TRUE, 
                        logAxes = "area", logT = log10)
  CI4 <- threshold_ci(fitT, method = "F", interval = NULL, Nboot = 3)

  #can't check actual boot answers as are random draws
  expect_equal(round(CI2[[1]], 1), c(0.4, 1.2))
  expect_equal(length(CI1[[1]]$ContOne), 3)
  expect_equal(length(CI1[[2]]$ContOne), 2)
  expect_true(is.numeric(CI3[[2]]$ZslopeOne[1]))
  expect_equal(length(CI3[[2]]), 2)
  expect_equal(round(CI4[[2]], 1), c(0.1, 0.7))
  expect_match(CI1$Method, "boot")
  expect_match(CI2[[2]], "F")
  expect_error(threshold_ci(fitT, method = "T"), 
               "method should be one of 'boot' or 'F'")
  expect_is(CI4, "sars")
})


test_that("get_coef returns correct results", {
  data(aegean2)
  a2 <- aegean2[1:168,]
  fitT <- sar_threshold(data = a2, mod = c("ContOne", "DiscOne", "ZslopeOne"),
                        interval = 0.1, non_th_models = TRUE, logAxes = "area", logT = log10)
  coefs <- get_coef(fitT)
  expect_equal(coefs[[1]], c(134.13,  54.31, 125.28))
  expect_equal(coefs[[3]], c(NA, NA, -208.65))
  expect_is(coefs, "sars")
  expect_error(get_coef(5), "fit object should be of class 'threshold'")
  fitT2 <- sar_threshold(data = a2, mod = c("ContTwo"),
                         interval = 1, non_th_models = FALSE, logAxes = "area")
  ContTwo <- fitT2[[1]][[1]]
  coefs2 <- get_coef(fitT2)
  p1 <- (predict(ContTwo, data.frame("x" = -2)) - predict(ContTwo, data.frame("x" = -4)))/2
  expect_equal(as.vector(round(p1,2)), coefs2[[2]])
  p2 <- (predict(ContTwo, data.frame("x" = 2)) - predict(ContTwo, data.frame("x" = 0)))/2
  expect_equal(as.vector(round(p2,2)), coefs2[[4]])
  p3 <- (predict(ContTwo, data.frame("x" = 6)) - predict(ContTwo, data.frame("x" = 4)))/2
  expect_equal(as.vector(round(p3,2)), coefs2[[6]])
})

