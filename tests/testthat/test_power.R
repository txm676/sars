context("sar_power")
library(sars)

test_that("sar_power returns correct results", {
  fit <- sar_power(galap)
  expect_equal(round(fit$AICc, 2), 189.03)
  expect_equal(as.vector(round(fit$par[2], 2)), 0.28)
  expect_is(fit, "sars")
  expect_match(fit$normaTest[[1]], "none")
  expect_error(sar_linear(5), "data must be a matrix or dataframe")
  fit2 <- sar_power(galap, homoTest = "cor.fitted")
  expect_equal(as.vector(round(fit2$homoTest[[2]]$p.value, 2)), 0.04)
  expect_match(fit2$homoTest[[2]]$method, "Spearman's rank correlation rho")
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, lcol = "black"))
  expect_error(plot(fit, col = "black"))
  #check works with tibble
  gp2 <- tibble::as_tibble(galap)
  expect_no_error(sar_power(data = gp2, grid_start = "none"))
})


#changed neg_expo function allowing z to be > 1, and also the change to 
#asymp setting z to Rplus rather than R
test_that("neg_expo and asymp returns correct results", {
  fit <- sar_negexpo(niering, grid_start = "none")
  expect_equal(round(fit$AICc, 2), 207.39)
  expect_equal(as.vector(round(fit$par[2], 2)), 26.14)
  fit2 <- sar_asymp(niering)
  expect_equal(round(fit2$BIC, 2), 204.32)
  expect_equal(as.vector(round(fit2$par[3], 5)), 0)
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, lcol = "black"))
  expect_error(plot(fit, col = "black"))
})


test_that("sar_power summary returns correct results", {
  fit <- sar_power(galap, normaTest = "lillie")
  fs <- summary(fit)
  expect_equal(round(sum(fs$residuals), 0), -31)#grid start so round to 0
  expect_output(str(fs), "List of 16")
  expect_is(fs, "summary.sars")
  expect_equal(round(fs$normaTest[[2]]$p.value, 3), 0.056)
})

test_that("sar_power returns warning for all identical species", {
  d <- data.frame("A" = 1:4, "S" = 0)
  #Note in R-devel 4.5.0, there is a change inside grepl, which 
  #made this warning matching no longer work unless the long string, 
  #which is split across lines, is put inside paste0()
  expect_warning(sar_power(d), paste0("All richness values are zero: ",
                 "parameter estimates of non-linear models should be ",
                  "interpreted with caution"))
  d$S <- 1
  expect_warning(sar_power(d), "All richness values identical")
})