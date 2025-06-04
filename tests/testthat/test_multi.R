context("sar_average")

test_that("sar_average returns correct results", {
  skip_on_cran()

  fit2 <- sar_average(data = galap, grid_start = "none")
  expect_equal(length(capture_output_lines(fit2, print = TRUE)),
               6)
  sfit2 <- summary(fit2)
  expect_equal(length(capture_output_lines(sfit2, print = T)),
               30)
  expect_equal(round(sum(fit2$mmi), 1), 1647.2)
  expect_no_error(plot(fit2))
  expect_no_error(plot(fit2, type = "bar"))
  expect_no_error(plot(fit2, type = "bar",
                       subset_weights = 0.1))
  expect_no_error(plot(fit2, allCurves = FALSE))
  expect_no_error(plot(fit2, mmSep = TRUE,
                       col.Sep = "orange",
                       lwd.Sep = 4))
  expect_no_error(plot(fit2, pcol = "red",
                       ModTitle = "A",
                       TiAdj = 0,
                       cex.main = 5,
                       cex.axis = 4,
                       pLeg = FALSE,
                       yRange = c(0,500)))

  fit3 <- sar_average(data = galap)#grid start on so round to 0
  expect_equal(round(sum(fit3$mmi), 0), 1647)
  expect_output(str(fit3), "List of 2")
  expect_is(fit3, "multi")
  expect_match(fit3$details$homo_test, "none")
  expect_match(fit3$details$norm_test, "none")
  expect_match(fit3$details$ic, "AICc")
  expect_error(sar_multi(5), "data must be a matrix or dataframe")
  #check works with tibble
  gp2 <- tibble::as_tibble(galap)
  expect_no_error(sar_average(data = gp2, grid_start = "none"))
  fit4 <- sar_average(data = galap, normaTest = "lillie",
                      homoTest = "cor.fitted",
                      neg_check = FALSE)
  expect_equal(length(capture_output_lines(fit4, print = TRUE)),
               8)
  #one run in 100 gave 1638 (due to random nature of grid_start)
  expect_true(round(sum(fit4$mmi), 0) %in% c(1638, 1639))#grid_start on so round to 0
  expect_match(fit4$details$norm_test, "lillie")
  expect_match(fit4$details$homo_test, "cor.fitted")
  #every 200th run or so is 13 rather than 14
  expect_true(length(fit4$details$mod_names) %in% c(13,14))
  expect_error(sar_average(data = galap, homoTest = 4))
  expect_error(sar_average(data = galap, homoTest = "cor.fitted",
                           homoCor = "correlation"))

  expect_length(sars_models(), 21)
  expect_length(display_sars_models(), 6)
})

test_that("sar_average using fit_collection object works", {
  skip_on_cran()
  ff <- sar_multi(data = galap, obj = c("power", "p1", "loga", "monod",
                                        "linear"))
  expect_warning(sar_average(obj = ff, data = galap, normaTest = "none",
                    neg_check = FALSE, homoTest = "cor.fitted"))
  fit3 <- sar_average(obj = ff, data = galap)#grid start on so round to 0
  expect_equal(round(sum(fit3$mmi), 0), 1663)
  expect_output(str(fit3), "List of 2")
  expect_is(fit3, "multi")
  expect_match(fit3$details$homo_test, "none")
  expect_match(fit3$details$ic, "AICc")
  expect_error(sar_multi(5), "data must be a matrix or dataframe")
})

 test_that("sar_average correctly deals with only 1 or 2 good fits", {
  skip_on_cran()
  expect_message(sar_average(obj = c("power", "powerR"),
                             data = galap, grid_start = "none",
                             normaTest = "shapiro",
                             homoTest = "cor.fitted",
                             homoCor = "spearman",
                             verb = FALSE,
                             display = FALSE))
  uusr <- sar_average(obj = c("power", "powerR"),
              data = galap, grid_start = "none",
              normaTest = "shapiro",
              homoTest = "cor.fitted",
              homoCor = "spearman",
              verb = FALSE,
              display = FALSE)
  expect_is(uusr, "sars")
  expect_equal(length(uusr), 23)
  expect_error(sar_average(obj = c("power", "linear"),
                             data = galap, grid_start = "none",
                             normaTest = "shapiro",
                             homoTest = "cor.fitted",
                             homoCor = "spearman",
                             verb = FALSE,
                             display = FALSE))
})

test_that("confidence intervals are correct", {
  skip_on_cran()
  fit3 <- expect_warning(sar_average(data = galap, grid_start = "none",
                                     normaTest = "none", homoTest = "none",
                    neg_check = FALSE, confInt = TRUE, ciN = 20))
  ci <- fit3$details$confInt
  expect_equal(nrow(ci), nrow(galap))
  expect_output(str(ci), "data.frame")
  expect_true(all(ci[ ,1] < ci[ ,2]))
})



test_that("confidence intervals works with AICc", {
  skip_on_cran()
  fit3 <- expect_warning(sar_average(data = galap, grid_start = "none",
                                     normaTest = "none", homoTest = "none",
                      neg_check = FALSE, confInt = TRUE, crit = "AICc", ciN = 20))
  ci <- fit3$details$confInt
  expect_equal(nrow(ci), nrow(galap))
  expect_output(str(ci), "data.frame")
  expect_true(all(ci[ ,1] < ci[ ,2]))
})

test_that("asymptotes are correctly calculated", {
  skip_on_cran()
  x <- 1:10
  y <- 10*x^0.23
  y <- c(y, rep(17,5))
  x <- c(x, 11:15)
  df <- data.frame(x, y)
  fit5 <- sar_average(data = df)
  s2 <- summary(fit5)
  #6 or 7 as gompertz doesn't tend to fit properly in this
  #use case with partial grid_start
  expect_true(length(which(s2$Model_table$Asymptote)) %in%
                c(6,7))
})
