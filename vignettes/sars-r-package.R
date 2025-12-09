## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  out.width = "100%"
)
options(cli.unicode = FALSE)

## ----include = FALSE----------------------------------------------------------
library(sars)

## ----fig.width=6, fig.height=6------------------------------------------------
#load an example dataset (Preston, 1962), fit the logarithmic SAR model using
#the grid_start method of selecting starting parameter values, return a model 
#fit summary and plot the model fit. 
data(galap) 
fit <- sar_loga(data = galap, grid_start = "partial") 
summary(fit) 
plot(fit)

## ----fig.width=16, fig.height=12----------------------------------------------
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

## ----fig.width=7, fig.height=19-----------------------------------------------
#load an example dataset (Niering, 1963), run the ‘sar_average’ function
#using a vector of model names and with no model validation tests, and
#produce the plots in Figure 2 of the paper 
data(niering) 

#run the ‘sar_average’ function using a vector of model names, and with simply
#using the default model starting parameter estimates (grid_start = "none")
fit <- sar_average(data= niering, obj =c("power","loga","koba","logistic","monod",
                                         "negexpo","chapman","weibull3","asymp"),
grid_start = "none", normaTest = "none", homoTest = "none", neg_check = FALSE, 
confInt = TRUE, ciN = 50) #a message is provided indicating that one model 
#(asymp) could not be used in the confidence interval calculation

par(mfrow = c(3,1)) #plot all model fits and the multimodel SAR curve as a separate curve on top
plot(fit, ModTitle = "a) Multimodel SAR", mmSep = TRUE)

#plot the multimodel SAR curve (with confidence intervals; see explanation
#in the main text, above) on its own 
plot(fit, allCurves = FALSE, ModTitle =
      "c) Multimodel SAR with confidence intervals", confInt = TRUE)

#Barplot of the information criterion weights of each model 
plot(fit, type = "bar", ModTitle = "b) Model weights", cex.lab = 1.3)

## ----fig.width=6, fig.height=6------------------------------------------------
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
#grid_start set to "none" for speed
p3 <- sar_average(data = galap, grid_start = "none")
sar_pred(p3, area = c(5000, 10000))

## ----fig.width=6, fig.height=6------------------------------------------------

#load an example dataset, and fit the continuous two-threshold model 
#to the data (with area transformed using log to the base 10), using an 
#interval of 0.1 (for speed)
data(aegean2)
fit <- sar_threshold(data = aegean2, mod = c("ContTwo"), interval = 0.1, 
                     non_th_models = FALSE, logAxes = "area", con = 1,
                     logT = log10, nisl = NULL)

#generate model fitting summary table (generally more useful when fitting multiple models)
summary(fit)

#Plot the resultant model fit 
plot(fit, cex = 0.8, cex.main = 1, cex.lab = 0.8, pcol = "grey") 

## ----fig.width=6, fig.height=6------------------------------------------------

#Load the dataset and fit the three habitat-heterogeneity models
#(power model form) alongside the Arrhenius power model, and view
#the raw model fit objects
data(habitat)
s <- sar_habitat(data = habitat, modType = "power")
s

#Fit the three habitat-heterogeneity models in log–log power
#model form alongside the linear (log–log) power model, and
#compare the model fits using information criteria
s2 <- sar_habitat(data = habitat, modType = "power_log")
summary(s2)

#Generate barplot of model AICc weights
plot(s2) 

## ----fig.width=7.5, fig.height=6----------------------------------------------
# Fit the countryside SAR model (power form) to the data, and use the function’s
# starting parameter value selection procedure (setting 'gridStart' here to
# "none" just for speed). Abbreviations: AG = agricultural land, SH = shrubland,
# F = oak forest, UB = ubiquitous species.data(countryside)
data(countryside)
s3 <- sar_countryside(data = countryside, modType = "power",
                      gridStart = "none", 
                      habNam = c("AG", "SH", "F"), 
                      spNam = c("AG_Sp", "SH_Sp", "F_Sp", "UB_Sp"))

# Use the fitted models to predict the richness of a site with given amounts of
# agricultural land, shrub land and forest. The ‘Indiv_mods’ element shows the
# predicted number of species (based on the given model) in a given group in the
# site, while the ‘Total’ element is the predicted total number of species
# across all groups.
countryside_extrap(s3, area = c(1000, 1000, 1000))

#Plot the fitted individual SAR curves for each species group (Type 2 plot),
#including a total species richness curve, modifying various aspects of the
#plot (e.g., legend position). Use the ‘which’ argument to only generate the
#second plot (i.e., the plot for shrubland)
par(mar=c(5.1, 4.1, 4.1, 7.5), xpd=TRUE)
plot(s3, type = 2, totSp = TRUE,
     lcol = c("black", "aquamarine4", "#CC661AB3",
              "darkblue", "darkgrey"), pLeg = TRUE,
     legPos ="topright", legInset = c(-0.27,0.3), 
     lwd = 1.5, ModTitle = c("Agricultural land", 
                             "Shrubland", "Forest"),
     which = 2)

# For a given habitat and a fixed site area, plot the number of species in each
# species group as a function of the proportion of the given habitat in the site
# (Type 3 plot).  Use the ‘which’ argument to only generate the first plot
# (i.e., the plot for agricultural land)
plot(s3, type = 3, totSp = TRUE, 
     lcol = c("black", "aquamarine4","#CC661AB3",
              "darkblue", "darkgrey"), pLeg = TRUE, 
     legPos ="topright", legInset = c(-0.27,0.3), 
     lwd = 1.5, ModTitle = c("Agricultural land", 
                             "Shrubland", "Forest"), 
     which = 1)

## ----fig.width=12, fig.height=6-----------------------------------------------
# Generate an effective area plot (Type 4 plot), customising point type and size
# etc 
par(mar=c(5.1, 4.1, 4.1, 2.1), xpd=FALSE)
plot(s3, type = 4, pch = 16, lwd = 2, cex  = 0.75, cex.axis = 0.75, 
     cex.lab = 0.9, ModTitle = c("Agricultural species"), which = 1)

