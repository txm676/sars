
<!-- README.md is generated from README.Rmd. Please edit that file -->
sars
====

> *fit and compare species-area relationship models using multimodel inference*

**sars** provides functionality to fit twenty SAR model using non-linear regression, and to calculate multi-model averaged curves using various information criteria. The software also provides easy to use functionality to plot multi-model SAR curves and to generate confidence intervals using bootstrapping.

Installation
------------

You can install the released version of sars from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sars")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("txm676/sars")
```

Example
-------

Basic usage of **sars** will result in using two types of functions:

Fitting the power sar model (Arrhenius 1921) to the 'galapagos' data set

``` r
fit_pow <- sar_power(data = galap)
```

Attempting to fit all 20 sar models to the 'galapagos' data set and generate a multimodel inference

``` r
mm_galap <- multi_sars(data = galap)
#> -- fitting model:  sar_power 
#> -- fitting model:  sar_powerR 
#> -- fitting model:  sar_epm1 
#> -- fitting model:  sar_epm2 
#> -- fitting model:  sar_p1 
#> -- fitting model:  sar_p2 
#> -- fitting model:  sar_expo 
#> -- fitting model:  sar_koba 
#> -- fitting model:  sar_mmf 
#> -- fitting model:  sar_monod 
#> -- fitting model:  sar_negexpo 
#> -- fitting model:  sar_chapman 
#> -- fitting model:  sar_weibull3 
#> -- fitting model:  sar_asymp 
#> -- fitting model:  sar_ratio 
#> -- fitting model:  sar_gompertz 
#> -- fitting model:  sar_weibull4 
#> -- fitting model:  sar_betap 
#> -- fitting model:  sar_heleg
```

Most of 'fitting' functions have corresponding plot methods.

Fitting the exponential sar model (Gleason 1922) to the 'galapagos' data set and plot it

``` r
fit_expo <- sar_expo(data = galap)

plot(fit_expo)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

Troubleshoutting :bomb:
-----------------------

If, despite the :heart: brought during the programming of this R :package: and writing of this documentation, you have difficulties to install or run sars, if you have questions about the procedures or calculations, or if you want to report bugs :beetle:, do not hesitate to connect with us on [GitHub](https://github.com/txm676/sars).

References
==========

Arrhenius, Olof. 1921. “Species and Area.” *The Journal of Ecology* 9 (1). British Ecological Society: 95. doi:[10.2307/2255763](https://doi.org/10.2307/2255763).

Gleason, Henry Allan. 1922. “On the Relation Between Species and Area.” *Ecology* 3 (2). Ecological Society of America: 158–62. doi:[10.2307/1929150](https://doi.org/10.2307/1929150).
