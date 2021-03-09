context("new_information_criteria")
library(sars)

test_that("new ICs works for individual models", {
  skip_on_cran()
  #loga
  s2 <- sar_loga(galap)
  l <- c(s2$AIC, s2$BIC, s2$AICc)
  expect_equal(round(l,4), c(187.1908, 189.5086, 189.1908))
  #match with nls
  y <- galap$s
  x <- galap$a
  n <- stats::nls(y ~ c + z*log(x), start = list("c" = 0.29, "z" = 30))
  nl <- c(stats::AIC(n), stats::BIC(n))
  expect_equal(round(nl,4), c(187.1908, 189.5086))
  #check residuals
  expect_equal(as.vector(round(residuals(n), 2)), round(s2$residuals, 2))
  #try with different model (p1)
  s3 <- sar_p1(galap)
  n2 <- stats::nls(y ~ c * x^z *(exp(-d*x)), 
                  start = list("c" = 8, "z" = 0.67, "d" = 0.002))
  nl <- c(stats::AIC(n2), stats::BIC(n2))
  LL <- logLik(n2)
  K <- 3 + 1# 1 for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(nl,4), round(c(s3$AIC, s3$BIC, s3$AICc), 4))

  #test linear model
  s2 <- sar_linear(galap)
  l <- c(s2$AIC, s2$BIC, s2$AICc)
  ll <- lm(s ~ a, galap)
  llic <- c(AIC(ll), BIC(ll))
  LL <- logLik(ll)
  K <- 2 + 1 #for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  llic <- c(llic, AICcm)
  expect_equal(round(l,4), round(llic, 4))
  #check residuals
  expect_equal(as.vector(round(residuals(ll), 4)), round(s2$residuals, 4))
  
  #test chapman model
  s3 <- sar_chapman(galap)
  n2 <- stats::nls(y ~  d * (1 - exp(-z*x))^c, 
  start = list("d" = 2.080189e+02, "z" = 9.975764e-03, "c" = 6.794128e-01))
  nl <- c(stats::AIC(n2), stats::BIC(n2))
  LL <- logLik(n2)
  K <- 3 + 1# 1 for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(nl,4), round(c(s3$AIC, s3$BIC, s3$AICc), 4))
  
  #test logistic model
  s3 <- sar_logistic(galap)
  n2 <- stats::nls(y ~  d/(1 + exp(-z*x + c)), 
                   start = list("d" = 212.84164599, 
                                "z" = 0.02274824, "c" = 1.59457723))
  nl <- c(stats::AIC(n2), stats::BIC(n2))
  LL <- logLik(n2)
  K <- 3 + 1# 1 for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(nl,4), round(c(s3$AIC, s3$BIC, s3$AICc), 4))
  
  #test negexpo model (no grid_start)
  s3 <- sar_negexpo(galap, grid_start = "none")
  n2 <- stats::nls(y ~ d * (1 - exp(-z * x)), 
                   start = list("d" = 207.91232765, "z" = 0.01432358))
  nl <- c(stats::AIC(n2), stats::BIC(n2))
  LL <- logLik(n2)
  K <- 2 + 1# 1 for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(nl,4), round(c(s3$AIC, s3$BIC, s3$AICc), 4))
  
  #test asymp model
  s3 <- sar_asymp(galap)
  n2 <- stats::nls(y ~ d - c * z^x, 
          start = list("d" = 209.1301838, "c" = 184.9069161, "z" = 0.9886514 ))
  nl <- c(stats::AIC(n2), stats::BIC(n2))
  LL <- logLik(n2)
  K <- 3 + 1# 1 for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(nl,4), round(c(s3$AIC, s3$BIC, s3$AICc), 4))
})


test_that("news ICs works for multi models", {
  skip_on_cran()
  y <- galap$s
  x <- galap$a
  fit3 <- sar_multi(galap, obj = c("power", "loga", "linear"))
  #power
  ff <- c(fit3$power$AIC, fit3$power$BIC, fit3$power$AICc)
  n3 <- nls(y ~ c * x^z, start = list("c" = 30, "z" = 0.25))
  nl <- c(stats::AIC(n3), stats::BIC(n3))
  LL <- logLik(n3)
  K <- 2 + 1 #for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(ff,4), round(nl,4))
  #koba model
  fit4 <- sar_average(data = galap, obj = c("koba", "gompertz", "heleg"))
  ff4 <- c(fit4$details$fits$koba$AIC, fit4$details$fits$koba$BIC, 
           fit4$details$fits$koba$AICc)
  n4 <- nls(y ~ c * log(1 + x/z), start = list("c" = 40, "z" = 3.5))
  nl <- c(stats::AIC(n4), stats::BIC(n4))
  LL <- logLik(n4)
  K <- 2 + 1 #for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(ff4,4), round(nl,4))
  expect_equal(round(ff4,3), c(186.073, 188.391, 188.073))
  #check residuals (round to 0 as grid start can make v small changes)
  expect_equal(as.vector(round(residuals(n4), 0)), 
               round(fit4$details$fits$koba$residuals, 0))
  #Gompertz model
  ff4 <- c(fit4$details$fits$gompertz$AIC, fit4$details$fits$gompertz$BIC, 
           fit4$details$fits$gompertz$AICc)
  n4 <- nls(y ~ d * exp(-exp(-z * (x - c))), 
            start = list("d" = 210.88694840, "z" = 0.01680535, "c" = 38.04431280))
  nl <- c(stats::AIC(n4), stats::BIC(n4))
  LL <- logLik(n4)
  K <- 3 + 1 #for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(ff4,4), round(nl,4))
  #heleg model
  ff4 <- c(fit4$details$fits$heleg$AIC, fit4$details$fits$heleg$BIC, 
           fit4$details$fits$heleg$AICc)
  n4 <- nls(y ~ c/(f + x^(-z)), 
            start = list("c" = 4.95780952, "f" = 0.02223783, "z" = 0.98822915))
  nl <- c(stats::AIC(n4), stats::BIC(n4))
  LL <- logLik(n4)
  K <- 3 + 1 #for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(ff4,4), round(nl,4))
  ##check it works with grid_start turned off as well
  fit8 <- sar_average(data = galap, obj = c("koba", "gompertz", "ratio"),
                      grid_start = "none")
  ff8 <- c(fit8$details$fits$ratio$AIC, fit8$details$fits$ratio$BIC, 
           fit8$details$fits$ratio$AICc)
  n4 <- nls(y ~ (c + z * x)/(1 + d * x), 
            start = list("c" = 20.08886448, "z" = 3.48524535, "d" = 0.01529))
  nl <- c(stats::AIC(n4), stats::BIC(n4))
  LL <- logLik(n4)
  K <- 3 + 1 #for variance
  n <- nrow(galap)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(ff8,4), round(nl,4))
  
  #weibull4 - with this and the betap model I can't get nls to converge at all,
  #even for sigmoid shaped fits, so I've used the warnOnly option and the port
  #algorithm to get it to return a non-converged object just for purposes of
  #this IC test. The parameters it returns differ from sars only at fourth
  #decimal point, so round by 3 instead.
  test <- data.frame("A" = c( 1,  2,  3,  5,  6,  7,  8, 12, 15, 20, 
                              23, 26, 27, 30, 34, 37, 40, 50),
                     "R" =  c(1,  1,  1,  2,  1,  2,  1,  4,  7,  9, 12, 14, 
                              15, 16, 16, 17, 16, 16))
  y2 <- test$R
  x2 <- test$A
  fit5 <- sar_average(data = test, obj = c("power", "loga", "weibull4"))
  ff4 <- c(fit5$details$fits$weibull4$AIC, fit5$details$fits$weibull4$BIC, 
           fit5$details$fits$weibull4$AICc)
  n4 <- suppressWarnings(nls(y2 ~ d * (1 - exp(-c * x2^z))^f, 
            start = list("d" = 1.626089e+01, "c" = 1.665468e-16, 
                         "z" = 1.092198e+01, "f" = 1.456715e-01),
            control = list("warnOnly" = T), algorithm = "port"))
  
  nl <- c(stats::AIC(n4), stats::BIC(n4))
  LL <- logLik(n4)
  K <- 4 + 1 #for variance
  n <- nrow(test)
  AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
  nl <- c(nl, AICcm)
  expect_equal(round(ff4,3), round(nl,3))
}) 
  
#hashed out for speed

# test_that("news ICs works for betap", {
#   skip_on_cran()
#   #betap model - using sars without grid_start gives suboptimal
#   #parameter estimates so have to use it with that turned on for it
#   #to match with nls
#   y <- galap$s
#   x <- galap$a
#   fit6 <- sar_betap(galap, grid_start = "exhaustive", grid_n = 100)
#   ff4 <- c(fit6$AIC, fit6$BIC,
#            fit6$AICc)
#   n4 <- suppressWarnings(nls(y ~ d * (1 - (1 + (x/c)^z)^-f),
#             start = list("d" = 207.993, "c" = 3455203.376,
#                          "z" = 0.832, "f" = 8109.855),
#             control = list("warnOnly" = T), algorithm = "port"))
#   nl <- c(stats::AIC(n4), stats::BIC(n4))
#   LL <- logLik(n4)
#   K <- 4 + 1 #for variance
#   n <- nrow(galap)
#   AICcm <- -2*LL+2*K*(n/(n-K-1))#- function taken directly from AICcmodavg
#   nl <- c(nl, AICcm)
#   expect_equal(round(ff4,3), round(nl,3))
# })
