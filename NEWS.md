
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
  
