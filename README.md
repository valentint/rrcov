
<!-- README.md is generated from README.Rmd. Please edit that file -->

R package providing scalable robust estimators with high breakdown
point.

## Installation

You can install ‘rrcov’ from github with:

``` r
# install.packages("remotes")
remotes::install_github("valentint/rrcov")
```

## Example

This is a basic example which shows you if the package is properly
installed:

``` r

library(rrcov)
#> Loading required package: robustbase
#> Scalable Robust Estimators with High Breakdown Point (version 1.5-5)
data(hbk)
(out <- CovMcd(hbk))
#> 
#> Call:
#> CovMcd(x = hbk)
#> -> Method:  Fast MCD(alpha=0.5 ==> h=40); nsamp = 500; (n,k)mini = (300,5) 
#> 
#> Robust Estimate of Location: 
#>       X1        X2        X3         Y  
#>  1.50345   1.85345   1.68276  -0.06552  
#> 
#> Robust Estimate of Covariance: 
#>     X1        X2        X3        Y       
#> X1   1.56742   0.15447   0.28699   0.16560
#> X2   0.15447   1.60912   0.22130  -0.01917
#> X3   0.28699   0.22130   1.55468  -0.21853
#> Y    0.16560  -0.01917  -0.21853   0.45091
```
