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
* local R installation, Windows 10, R 4.0.3
* local Linux Mint 20 installation, R 4.0.3 with valgrind
* rocker/r-devel-ubsan-clang docker image 
* Ubuntu 16.04 (on travis-ci and r-hub), R 4.0.3
* Debian Linux (r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (r-hub)
* Fedora Linux (r-hub)
* Oracle Solaris 10, x86, 32 bit, R-release (r-hub)
* Oracle Solaris 10, x86, 32 bit, R-release, Oracle Developer Studio 12.6 (r-hub)
* macOS 10.13.6 High Sierra, R-release, brew (r-hub)
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (r-hub)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTES.

