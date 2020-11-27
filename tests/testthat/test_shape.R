context("sar_shape")

test_that("observed_shape info is correct", {
  data(galap)
  fit <- sar_epm2(galap)
  fit2 <- sar_power(galap)
  fit3 <- sar_average(data = galap)
  s3 <- summary(fit3)
  x <- fit$observed_shape
  x3 <- s3$Model_table$Shape[17]#epm2 again
  #in R 3.6.3 it returns the s3..[17] as a factor but later
  #versions of R return it as a character vector. The as.vector()
  #is thus here to pass Travis test for past versions.
  expect_match(as.vector(s3$Model_table$Model[17]), "epm2")
  expect_match(x, "sigmoid")
  expect_match(fit2$observed_shape, "convex up")
  expect_match(x3, "sigmoid")
  #convex down test
  test <- data.frame("A" = c( 1,  2,  3,  5,  6,  7,  8, 12, 15, 20, 
             23, 26, 27, 30, 34, 37, 40, 50),
             "R" =  c(30, 25, 22, 19, 18, 17, 17, 15, 14, 13,
             12, 12, 12, 12, 11, 11, 11, 10))
  fit4 <- sar_power(test)
  expect_match(fit4$observed_shape, "convex down")
  fit5 <- sar_p1(test)
  expect_match(fit4$observed_shape, "convex down")
  #sigmoid test
  test$R <- c(1,  1,  1,  2,  1,  2,  1,  4,  7,  9, 12, 14, 15, 16, 16, 17, 16,
              16)
  fit6 <- sar_weibull4(test)
  fit7 <- sar_p1(test)
  fit8 <- sar_loga(test)
  fit9 <- sar_linear(test)
  expect_match(fit6$observed_shape, "sigmoid")
  expect_match(fit7$observed_shape, "sigmoid")
  expect_match(fit8$observed_shape, "convex up")
  expect_match(fit9$observed_shape, "linear")
})
