# hermiter v2.1.0

## New features

* A method has been added for estimating the Kendall rank correlation 
coefficient in the bivariate setting.
* The univariate quantile estimation method has been significantly enhanced in
accuracy using series acceleration techniques. Series acceleration is enabled
by default.
* The univariate pdf and cdf methods have been significantly enhanced in
accuracy using series acceleration techniques. Series acceleration is enabled
by default.
* The new default method for the univariate quantile estimation method, 
'interpolate' is much faster than the alternate method, 'bisection' with nearly
the same accuracy.
* Added print and summary methods for both the univariate and bivariate 
hermite_estimator objects.
* Convenience function added to calculate sums of Hermite functions.

## Documentation improvements

* The vignette `hermiter`, namely `vignette("hermiter")` has been extended to 
included examples pertaining to estimation of the Kendall Tau nonparametric
correlation coefficient in the bivariate setting.


## Test environments
* local R installation, Windows 10, R 4.1.2
* local R installation, Ubuntu Linux 21.10, R 4.1.2
* Ubuntu Linux 20.04.1 LTS (r-hub)
* Debian Linux, R-release, GCC (r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (r-hub)
* Fedora Linux (r-hub)
* Apple Silicon (M1), macOS 11.6 Big Sur, R-release
* macOS 10.13.6 High Sierra, R-release, CRAN's setup
* Oracle Solaris 10, x86, 32 bit, R-release (r-hub)
* win-builder (devel and release)

## R CMD check results

There was 1 NOTE:

Found the following (possibly) invalid URLs:
  URL: https://projecteuclid.org/euclid.ejs/1488531636
    From: inst/doc/hermiter.html
    Status: 500
    Message: Internal Server Error
  URL: https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-11/issue-1/Sequential-quantiles-via-Hermite-series-density-estimation/10.1214/17-EJS1245.full
    From: README.md
    Status: 500
    Message: Internal Server Error

Found the following (possibly) invalid DOIs:
  DOI: 10.1214/17-EJS1245
    From: DESCRIPTION
    Status: Internal Server Error
    Message: 500
    
    
* The URL and DOI listed above appear to be valid and work when tested directly.

