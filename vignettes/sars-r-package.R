## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  out.width = "100%"
)
options(cli.unicode = FALSE)

## ----include = FALSE----------------------------------------------------------
library(sars)

## ---- fig.width=6, fig.height=6-----------------------------------------------
#load an example dataset (Preston, 1962), fit the logarithmic SAR model using
#the grid_start method of selecting starting parameter values, return a model 
#fit summary and plot the model fit. 
data(galap) 
fit <- sar_loga(data = galap, grid_start = TRUE, grid_n = 100) 
summary(fit) 
plot(fit)

## ---- fig.width=16, fig.height=12---------------------------------------------
#Create a fit_collection object containing multiple SAR model fits, and 
#plot all fits. 
fitC <- sar_multi(data = galap, obj = c("power", "loga", "monod"))
plot(fitC) #see Fig.1

## -----------------------------------------------------------------------------
#load an example dataset, fit the linear SAR model whilst running residual
#normality and homogeneity tests, and return the results of the residual
#normality test 
data(galap) 
fit <- sar_linear(data = galap, normaTest ="lillie", homoTest = "cor.fitted") 
summary(fit) #a warning is provided  indicating the normality test failed 
fit$normaTest

## ---- fig.width=7, fig.height=19----------------------------------------------
#load an example dataset (Niering, 1963), run the ‘sar_average’ function
#using a vector of model names and with no model validation tests, and
#produce the plots in Figure 2 of the paper 
data(niering) 

#run the ‘sar_average’ function using a vector of model names 
fit <- sar_average(data= niering, obj =c("power","loga","koba","mmf","monod",
                                         "negexpo","chapman","weibull3","asymp"),
normaTest = "none", homoTest = "none", neg_check = FALSE, confInt = TRUE, ciN
= 50) #a message is provided indicating that one model (asymp) could not be
#fitted

par(mfrow = c(3,1)) #plot all model fits and the multimodel SAR curve as a separate curve on top
plot(fit, ModTitle = "a) Multimodel SAR", mmSep = TRUE)

#plot the multimodel SAR curve (with confidence intervals; see explanation
#in the main text, above) on its own 
plot(fit, allCurves = FALSE, ModTitle =
      "c) Multimodel SAR with confidence intervals", confInt = TRUE)

#Barplot of the information criterion weights of each model 
plot(fit, type = "bar", ModTitle = "b) Model weights", cex.lab = 1.3)

## ---- fig.width=6, fig.height=6-----------------------------------------------
#load an example dataset, fit the log-log power model, return a model fit
#summary and plot the model fit. When ‘compare’ == TRUE, the non-linear
#power model is also fitted and the resultant parameter values compared. 
#If any islands have zero species, a constant (‘con’) is added to all 
#species richness values. 
data(galap) 
fit <- lin_pow(dat = galap, compare = TRUE, con = 1) 
summary(fit) 
plot(fit)

#load an example dataset, fit the random placement model and plot the 
#model fit and standard deviation. The ‘data’ argument requires a species-
#site abundance matrix: rows are species and columns are sites. The area 
#argument requires a vector of site (island) area values. 
data(cole_sim) 
fit <- coleman(data = cole_sim[[1]], area = cole_sim[[2]]) 
plot(fit, ModTitle = "Hetfield")

#load an example dataset, fit the GDM using the logarithmic SAR model, and
#compare the GDM with three alternative (nested) models: area and time 
#(age of each island), area only, and intercept only. 
data(galap) 
galap$t <- rgamma(16, 5, scale = 2)#add a random time variable 
gdm(data = galap, model = "loga", mod_sel = TRUE)

## -----------------------------------------------------------------------------
#fit the power model and predict richness on an island of area = 5000
data(galap)
p <- sar_power(data = galap)
sar_pred(p, area = 5000)

#fit three SAR models and predict richness on islands of area = 5000 & 10000
p2 <- sar_multi(galap, obj = c("power", "loga", "koba"))
sar_pred(p2, area = c(5000, 10000))

#calculate a multi-model curve and predict richness on islands of area = 5000 & 10000
p3 <- sar_average(data = galap)
sar_pred(p3, area = c(5000, 10000))

## ---- fig.width=6, fig.height=6-----------------------------------------------

#load an example dataset, and fit the continuous two-threshold model 
#to the data (with area transformed using log to the base 10), using an 
#interval of 0.1 (for speed) and parallel processing.
data(aegean2)
fit <- sar_threshold(data = aegean2, mod = c("ContTwo"), interval = 0.1, 
                     non_th_models = FALSE, logAxes = "area", con = 1,
                     logT = log10, nisl = NULL, parallel = TRUE, cores = 3)

#generate model fitting summary table (generally more useful when fitting multiple models)
summary(fit)

#Plot the resultant model fit 
plot(fit, cex = 0.8, cex.main = 1, cex.lab = 0.8, pcol = "grey") 

