context("indiv_par_CIs")
library(sars)

test_that("new ICs works for individual model parameter estimates", {
  skip_on_cran()
  s2 <- sar_loga(galap)
  s3 <- sar_power(galap)
  #check power model sigConf has the extra columns for nls fits
  expect_equal(ncol(s2$sigConf), 6)
  expect_equal(ncol(s3$sigConf), 9)
  #check the c and z-value returned from our fit and nls are the same
  expect_equal(round(s3$sigConf[,"Estimate"], 2), round(s3$sigConf[,"nls.Est."], 2))
  
  #check par values in sigConf match those stored in $par
  expect_equal(s3$sigConf[,"Estimate"],s3$par) 
  expect_equal(s2$sigConf[,"Estimate"],as.vector(s2$par)) 
  y <- galap$s
  x <- galap$a
  n3 <- nls(y ~ c * x^z, start = list("c" = s3$par[1], "z" = s3$par[2]))
  #check nls fit pars match the nls pars stored in sigCof
  expect_equal(round(coef(n3), 2), round(s3$sigConf[,"nls.Est."], 2))
  sn3 <- summary(n3)$coefficients
  #check that the other coefficient values returned from the nls fit match those
  #in our table; they wont be exactly the same though as often small differences
  #in the par estimates (hence the round to 0) from ours and nls (but should
  #only be at decimal places)
  expect_equal(round(sn3[,"Std. Error"], 0), round(s3$sigConf[,"Std. Error"], 0)) 
  expect_equal(round(sn3[,"t value"], 0), round(s3$sigConf[,"t value"], 0))
  expect_equal(round(sn3[,"Pr(>|t|)"], 0), round(s3$sigConf[,"Pr(>|t|)"], 0))
  
  #Repeat for heleg
  #check par values in sigConf match those stored in $par
  s4 <- sar_heleg(niering)
  expect_equal(s4$sigConf[,"Estimate"], as.vector(s4$par))
  y <- niering$s
  x <- niering$a
  n4 <- nls(y ~ c/(f + x^(-z)), 
                    start = list("c" = s4$par[1], "f" = s4$par[2], 
                                 "z" = s4$par[3]))
  sn4 <- summary(n4)$coefficients
  #check that the other coefficient values returned from the nls fit match those
  #in our table; they wont be exactly the same though as often small differences
  #in the par estimates (hence the round to -1) from ours and nls (but should
  #only be at decimal places)
  #first par v large and so can differ by up to 1, rather than 1 decimal place;
  #thus exclude this from standard error calculation
  expect_equal(as.vector(round(sn4[2:3,"Std. Error"], 0)), 
               round(s4$sigConf[2:3,"Std. Error"], 0)) 
  expect_equal(as.vector(round(sn4[,"t value"], 0)), 
               round(s4$sigConf[,"t value"], 0))
  expect_equal(as.vector(round(sn4[,"Pr(>|t|)"], 0)), 
               round(s4$sigConf[,"Pr(>|t|)"], 0))
  
  #Repeat for logistic
  #check par values in sigConf match those stored in $par
  s5 <- sar_logistic(aegean)
  expect_equal(s5$sigConf[,"Estimate"], as.vector(s5$par))
  y <- aegean$s
  x <- aegean$a
  n5 <- nls(y ~ d/(1 + exp(-z * x + c)), 
            start = list("d" = s5$par[1], "z" = s5$par[2], 
                         "c" = s5$par[3]))
  sn5 <- summary(n5)$coefficients
  #check that the other coefficient values returned from the nls fit match those
  #in our table; they wont be exactly the same though as often small differences
  #in the par estimates (hence the round to 0) from ours and nls (but should
  #only be at decimal places)
  #first par large and so can differ by up to 1, rather than 1 decimal place;
  #thus exclude this from standard error calculation
  expect_equal(as.vector(round(sn5[2:3,"Std. Error"], 0)), 
               round(s5$sigConf[2:3,"Std. Error"], 0)) 
  expect_equal(as.vector(round(sn5[,"t value"], 0)), 
               round(s5$sigConf[,"t value"], 0))
  expect_equal(as.vector(round(sn5[,"Pr(>|t|)"], 0)), 
               round(s5$sigConf[,"Pr(>|t|)"], 0))
})