
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The **‘sars’** R Package <img src="man/figures/sars_logo.png" align="right" width="10%"/>

[![Build
Status](https://travis-ci.org/txm676/sars.svg?branch=master)](https://travis-ci.org/txm676/sars)
[![Downloads](https://cranlogs.r-pkg.org/badges/sars?color=brightgreen)](https://cran.r-project.org/package=sars)
[![CRAN](https://www.r-pkg.org/badges/version/sars)](https://cran.r-project.org/package=sars)
[![codecov.io](https://codecov.io/github/txm676/sars/coverage.svg?branch=master)](https://codecov.io/github/txm676/sars?branch=master)
[![DOI](https://zenodo.org/badge/95295704.svg)](https://zenodo.org/badge/latestdoi/95295704)

> *fit and compare **Species-Area Relationship (SAR)** models using
> multi-model inference*

**sars** provides functionality to fit twenty SAR model using non-linear
regression, and to calculate multi-model averaged curves using various
information criteria. The software also provides easy to use
functionality to plot multi-model SAR curves and to generate confidence
intervals using bootstrapping. Additional SAR related functions include
fitting the linear version of the power model and comparing parameters
with the non-linear version, fitting the general dynamic model of island
biogeography, and fitting the random placement model to a species
abundance - site matrix.

As this is version 1.1.1 of the package, it is possible that there are
some bugs in places. Please report any issues to us via GitHub.

The package has an associated vignette that provides examples of how to
use the package.

A website for the package can be found here:
<https://txm676.github.io/sars/>

Version 1.1.1 of the package has been archived on the Zenodo research
data repository (DOI: 10.5281/zenodo.2573067).

## Table of Contents

  - [Installation](#installation)
  - [Example](#example-usage)
  - [References](#references)

## Installation

You can install the released version of sars from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("sars")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("txm676/sars")
```

## Example usage

Basic usage of **sars** will result in using two types of functions:

To fit the power sar model (Arrhenius 1921) to the ‘galapagos’ (Preston
1962) data set:

``` r
fit_pow <- sar_power(data = galap)
```

Attempting to fit all 20 sar models to the ‘galapagos’ (Preston 1962)
data set and get a multi-model SAR:

``` r
mm_galap <- sar_average(data = galap)
#> 
#>  Now attempting to fit the 20 SAR models: 
#> 
#> --  multi_sars ----------------------------------------------------------------- multi-model SAR --
#> > power    : v
#> > powerR   : v
#> > epm1     : v
#> observed shape algorithm failed: observed shape set to
#>                 theoretical shape (sigmoid)
#> > epm2     : v
#> > p1       : v
#> > p2       : v
#> > loga     : v
#> > koba     : v
#> > mmf      : v
#> > monod    : v
#> > negexpo  : v
#> > chapman  : Warning: could not compute parameters statistics
#> > weibull3 : v
#> > asymp    : v
#> > ratio    : v
#> > gompertz : v
#> > weibull4 : v
#> > betap    : v
#> > heleg    : v
#> > linear   : v
#> 
#> Model fitting completed - all models succesfully fitted.  Now undertaking model validation checks.
#> Additional models  will be excluded if necessary:
#> 
#> 4 models failed the residuals normality test and  have been excluded  from the multi SAR:
#> Extended Power model 1, Asymptotic regression, Cumulative Weibull 4 par., Linear model
#> 16 remaining models used to construct the multi  SAR:
#>  Power, PowerR, Extended Power model 2, Persistence function 1, Persistence function 2, Logarithmic, Kobayashi, MMF, Monod, Negative exponential, Chapman Richards, Cumulative Weibull 3 par., Rational function, Gompertz, Beta-P cumulative, Heleg(Logistic) 
#> ---------------------------------------------------------------------------------------------------
```

Each of the ‘fitted’ objects have corresponding plot methods:

to fit the logarithmic SAR model (Gleason 1922) to the ‘galapagos’ data
set and plot it

``` r
fit_loga <- sar_loga(data = galap)

plot(fit_loga)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

to fit a multimodel SAR curve to the ‘galapagos’ data set and plot it

``` r
mm_galap <- suppressMessages(sar_average(data = galap, verb = FALSE))

plot(mm_galap, pLeg = FALSE)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

## References

<div id="refs" class="references">

<div id="ref-Arrhenius1921">

Arrhenius, Olof. 1921. “Species and Area.” *The Journal of Ecology* 9
(1). British Ecological Society:95. <https://doi.org/10.2307/2255763>.

</div>

<div id="ref-Gleason1922">

Gleason, Henry Allan. 1922. “On the Relation Between Species and Area.”
*Ecology* 3 (2). Ecological Society of America:158–62.
<https://doi.org/10.2307/1929150>.

</div>

<div id="ref-Preston1962">

Preston, F. W. 1962. “The Canonical Distribution of Commonness and
Rarity: Part I.” *Ecology* 43 (2). Ecological Society of America:185.
<https://doi.org/10.2307/1931976>.

</div>

</div>
