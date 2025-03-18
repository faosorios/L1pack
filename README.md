<!-- README.md, last updated 2023-01-29, F.Osorio -->

# L1pack: Routines for L1 estimation

<!-- badges: start -->
[![CRAN status](http://www.r-pkg.org/badges/version/L1pack)](https://cran.r-project.org/package=L1pack)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/L1pack)](https://cran.r-project.org/package=L1pack)
<!-- badges: end -->

**L1pack** provides routines to perform L1 estimation in linear regression models, estimation of mean and covariance matrix using the multivariate Laplace distribution, and multivariate random number generation for the Laplace distribution. A basic set of methods for printing the results is also available.

## Reference Manual

* [L1pack.pdf](https://cran.r-project.org/web/packages/L1pack/L1pack.pdf)

## Resources

Latest binaries and sources can be found at the [CRAN package repository](https://cran.r-project.org/package=L1pack):

* [L1pack_0.52.tar.gz](https://cran.r-project.org/src/contrib/L1pack_0.52.tar.gz) - Package sources
* [L1pack_0.50.zip](https://cran.r-project.org/bin/windows/contrib/4.4/L1pack_0.50.zip) - Windows binaries (R-release)
* [L1pack_0.52.tgz](https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.4/L1pack_0.52.tgz) - MacOS binaries (R-release, arm64)
* [L1pack_0.52.tgz](https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.4/L1pack_0.52.tgz) - MacOS binaries (R-release, x86_64)

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
Alternatively, you can download the source as a tarball or as a zip file. Unpack the tarball or zipfile (thereby creating a directory named, L1pack) and install the package source by executing (at the console prompt)

``` r
R CMD INSTALL L1pack
```
Next, you can load the package by using the command `library(L1pack)`

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
#>   Osorio, F., Wolodzko, T. (2025). Routines for L1 estimation. R
#>   package version 0.52. URL: https://github.com/faosorios/L1pack
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>    title = {Routines for L1 estimation},
#>    author = {F. Osorio and T. Wolodzko},
#>    year = {2025},
#>    note = {R package version 0.52},
#>    url = {https://github.com/faosorios/L1pack},
#>   }
```
## Some papers using L1pack
- Plate, M., Bernstein, R., Hoppe, A., Bienefeld, K. (2019). Comparison of infinitesimal and finite locus models for long-term breeding simulations with direct and maternal effects at the example of honeybess. [PLOS ONE](https://doi.org/10.1371/journal.pone.0213270) 14, e0213270.
- Wang, W., Yu, P., Lin, L., Tong, T. (2019). Robust estimation of derivatives using locally weighted least absolute deviation regression. [Journal of Machine Learning Research](http://jmlr.org/papers/v20/17-340.html) 20 (60), 1-49.
- VandenHeuvel, D., Wu, J., Wang, Y.G. (2023). Robust regression for electricity demand forecast against cyberattacks. [International Journal of Forecasting](https://doi.org/10.1016/j.ijforecast.2022.10.004) 39, 1573-1592

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
