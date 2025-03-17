<!-- README.md, last updated 2023-01-29, F.Osorio -->

# L1pack: Routines for L1 estimation

<!-- badges: start -->
[![CRAN status](http://www.r-pkg.org/badges/version/L1pack)](https://cran.r-project.org/package=L1pack)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/L1pack)](https://cran.r-project.org/package=L1pack)
<!-- badges: end -->

Provides routines to perform L1 estimation in linear regression models. Estimation of mean and covariance matrix using the multivariate Laplace distribution, and multivariate random number generation for the Laplace distribution. A basic set of methods for printing the results is also available.

## Reference Manual

* [L1pack.pdf](https://cran.r-project.org/web/packages/L1pack/L1pack.pdf)

## Resources

Latest binaries and sources can be found at the [CRAN package repository](https://cran.r-project.org/package=L1pack):

* [L1pack_0.50.tar.gz](https://cran.r-project.org/src/contrib/L1pack_0.50.tar.gz) - Package sources
* [L1pack_0.50.zip](https://cran.r-project.org/bin/windows/contrib/4.4/L1pack_0.50.zip) - Windows binaries (R-release)
* [L1pack_0.50.tgz](https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.4/L1pack_0.50.tgz) - MacOS binaries (R-release, arm64)
* [L1pack_0.50.tgz](https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.4/L1pack_0.50.tgz) - MacOS binaries (R-release, x86_64)

## Installation

Install `L1pack` from CRAN using.

``` r
install.packages("L1pack")
```

You can install the latest development version from github with:

``` r
# install.packages("devtools")
devtools::install_github("faosorios/L1pack")
```

## Methods

The methods implemented in `L1pack` include:

-   Barrodale and Roberts (1974) procedure for L1 estimation in linear regression.
-   EM algorithm for LAD estimation in linear regression (Phillips, 2002).
-   Estimation of center and Scatter matrix using the multivariate Laplace distribution.
-   Density, distribution function, quantile function and random number generation for univariate Laplace distribution.
-   Density and random number generation for the multivariate Laplace distribution (Gomez et al., 1998).
-   Computation of the generalized spatial median estimator as defined by Rao (1988)

## Citation

To cite package `L1pack` in publications use:

``` r
citation("L1pack")
#> 
#> To cite package L1pack in publications use:
#> 
#>   Osorio, F., Wolodzko, T. (2024). Routines for L1 estimation. R
#>   package version 0.50. URL: https://github.com/faosorios/L1pack
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>    title = {Routines for L1 estimation},
#>    author = {F. Osorio and T. Wolodzko},
#>    year = {2024},
#>    note = {R package version 0.50},
#>    url = {https://github.com/faosorios/L1pack},
#>   }
```
## Providing Feedback

Please report any bugs/suggestions/improvements to [Felipe Osorio](https://faosorios.github.io/). 
If you find these routines useful or not then please let me know. Also, acknowledgement 
of the use of the routines is appreciated.

## About the Authors

Felipe Osorio is an applied statistician and creator of several R packages
* Webpage: [faosorios.github.io](https://faosorios.github.io/)

Tymoteusz Wolodzko is Software Engineer and Machine Learning Engineer. Also is elected moderator 
for [CrossValidated.com](https://stats.stackexchange.com/)
* Webpage: [twolodzko.github.io](https://twolodzko.github.io/)
