context("lin_pow")

test_that("lin_pow returns correct results", {
   fit <- lin_pow(galap, con = 1, normaTest = "lillie",
                  homoTest = "cor.fitted", homoCor = "spearman")
   expect_equal(round(fit$Model$coefficients[2], 2), 0.34)
   expect_equal(round(fit$normaTest[[2]]$p.value, 2), 0.35)
   expect_equal(round(fit$homoTest[[2]]$p.value, 2), 0.13)
   #checking homotest matches cor.test done manually
   fit2 <- lin_pow(aegean, homoTest = "cor.area", homoCor = "pearson")
   AO <- aegean[order(aegean$a),]#is ordered inside function
   rs2 <-  fit2$Model$residuals^2
   rs3 <-  summary(fit2)$Model$residuals^2
   pp <- cor.test(rs2, log(AO$a), method = "pearson")
   pp2 <- cor.test(rs3, log(AO$a), method = "pearson")
   expect_equal(fit2$homoTest[[2]]$p.value, pp$p.value)
   expect_equal(fit2$homoTest[[2]]$p.value, pp2$p.value)
   expect_no_error(plot(fit))
   expect_equal(length(capture_output_lines(fit, print = TRUE)),
                11)
   sfit2 <- summary(fit)
   expect_equal(length(capture_output_lines(sfit2, print = TRUE)),
                22)
})

test_that("lin_pow log-transformation works",{
  fit <- lin_pow(galap, con = 1)
  fit2 <- lin_pow(galap, con = 1, logT = log10)
  fit3 <- summary(lin_pow(galap, logT = log2))
  expect_match(fit$logT, "log()")
  expect_match(fit2$logT, "log10()")
  expect_match(fit3$logT, "log2()")
  expect_equal(round(fit$Model$coefficients[1], 2), 3.02)
  expect_equal(round(fit2$Model$coefficients[1], 2), 1.31)
  expect_equal(round(fit$Model$coefficients[2], 2), 
               round(fit2$Model$coefficients[2], 2), 
               round(fit3$Model$coefficients[2], 2))
})


test_that("lin_pow returns warning for all identical species", {
  d <- data.frame("A" = 1:4, "S" = 0)
  #Note in R-devel 4.5.0, there is a change inside grepl, which 
  #made this warning matching no longer work unless the long string, 
  #which is split across lines, is put inside paste0()
  expect_warning(lin_pow(d, compare = TRUE), paste0("All richness",
                 " values are zero: ",
                 "parameter estimates of non-linear models should be ",
                 "interpreted with caution"))
  expect_warning(lin_pow(d, compare = FALSE), "All richness values identical")
  d$S <- 1
  expect_warning(lin_pow(d), "All richness values identical")
})