
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCOR

<!-- badges: start -->

<!-- badges: end -->

The goal of SCOR is to solve optimization problem where parameters
belong to the surface of unit hyper-sphere

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("synx21/SCOR")
```

## Example

This package consists function to run SCOR and functions to evaluate
EHUM, SHUM and ULBA estimates.

``` r
library(SCOR)
optimized_EHUM(rep(1,12),colnames(AL),AL)
#> [1] 0.8403805
```
