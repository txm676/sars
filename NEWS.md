
##Version 1.3.4
  * For the power model only we now return (in the sigConf object) the parameter
    estimates from a nls fit (using as starting par estimates the par values from
    our fitting process). Using the nls fit, we also return the confidence intervals
    around the par estimates generated using the confint function - these should be
    more accurate than the default sars CIs (2 * standard error)
  * Corrected bug in find_one_threshold_cont and find_one_threshold_disc which didn't
    properly report warning that multiple par estimates returned same minimum RSS

##Version 1.3.3
  *It was realised that the mmf and heleg models are (basically) identical and 
   almost always provide exactly the same fit. As such, the mmf model has been
   deprecated and replaced with the sar_logistic() model, i.e., the standard logistic
   model given in Tjorve (2003).

##version 1.3.2
  * Changed mmi confidence interval function to fit all the models the user
    originally selects to the boostrapped samples. Also so it takes all of the original           arguments (e.g.        normatest, IC crit etc)  provided by the user.
  * Now return residuals from the non-linear models as observed - fitted rather than vice        versa, to match lm and nls etc. 
  * Added power_area_time variant to gdm function, and changed the linear power GDM function     to fit the        orginal
    GDM model of Whittaker et al. (2008) - the ATT2 model (i.e. semi-log SAR)
  * Removed the check which assigned a model with R2 < 0 as having not converged.
  * Now return optim model convergence info for all models in a sar_average() fit (and its       summary             function).
  * Changed grid_start to have three options: none, partial (default) and exhaustive - key       change to           previous versions as now grid_start is implemented as default.

##version 1.3.1
  * Corrected bug in obs_shape function which meant it was not recognising sigmoid fits
  * this was linked to a bug in the function used to calculate 1st and 2nd derivatives (also     corrected)
  * Added new model shape category (convex/sigmoid) for epm1, asymp
    and p1 models (allowing them to be tested for observed sigmoid shape)
  * Corrected bug in plot function for sar_threshold meaning dataframes with more
    than 2 columns could not be plotted.
  * Changed gdm power function to return c rather than log c. And corrected a bug
    for the intercept only power gdm model.
  * Updated gdm function to return AICc and R2, and fit the GDM using the linear version of      the power and logarithmic models. Also implemented a start_vals argument to provide          starting values
  * Corrected bug in homogeneity of variance test - it was previously using the raw residuals     rather than the squared residuals. And lin_pow homog test was also using untransformed       area rather than           log(area).
  * Changed the defaults for residual normality and homogeneity tests to be "none" - then up     to 
    the user to turn them on and select which test they want.
  * Correcting number of parameters for DiscTwo in sar_threshold (from 8 to 9)
  * Corrected bug in Chapman model function - the ^c term was in the wrong place.

##version 1.3.0
  * Added a set of functions for fitting, evaluating and plotting a range of commonly 
    used piecewise SAR models (see help page and accompanying paper)

##version 1.2.3
  * return shape algorithm fail info to model summary table
  * added grid_start argument option to sar_average and sar_multi
  * edited grid_start to ensure very small starting par values are always included
  * edited how grid_start works and have added a grid_n argument
  * changed the negative exponential model fitting process so that the z parameter can
    be any number rather than constrained between 0 and 1.
  * changed the asymptotic model fitting process so that a negative z-value cannot be returned
  * bug corrections in the confidence interval function, and adapting it to work with
    grid_start
  * Changing model plotting to plot smooth curves by creating 1000 fitted values
    using fitted parameters

## Version 1.2.2
  * bug fix - no AICc option in confidence interval function
  * changed AIC, BIC and AICc equations to be calculated using
    the same approach as the nls and lm functions.

## Version 1.2.1
  * added a warning for when all richness values are identical
  * adding functionality to plot.multi for plotting the multimodel curve on top of the other model fits

## Version 1.2.0
  * bug fixes (bar plot of weight x-axis labels order was incorrect in some cases)
  * logT argument added to lin_pow function (can choose log-transformation function)

## Version 1.1.2
  * bug fixes
  * additional tests added
  * sar_pred() function added for SAR extrapolation
  * negative fitted value check (fit$neg_check) returned for each individual model fit
  * knitr error correction (for CRAN)
  
## Version 1.1.1
 * added vignette
 * bug fixes
 * internal updates
 * Version archived on Zenodo

## Version 1.1.0  
  * sar_multi() split into two functions: sar_multi() and sar_average()
  * New sar_multi() replaces fit_collection()
  * unit testing added
  * functions added to bring up a list of all models with information
  * bug fixes
  
  
## Version 1.0.0

  * Initial version on CRAN
  
