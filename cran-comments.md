# hermiter v2.0.1

Fixed minor issues in c++ source code causing the build to fail on Solaris 
along with ASAN/valgrind issues. Updated the test environments as listed below. 
Added appropriate test cases. This is a re-submission.

# hermiter v2.0.0

This is a major update to include bivariate estimators.  

## Breaking changes

* `get_coefficients` has been removed as it is redundant.
* `combine_hermite` has been renamed to `merge_hermite` for clarity.
* `combine_pair` has been renamed to `merge_pair` for clarity.
* `hermite_integral_val_quantile_adap` has been renamed to 
`hermite_integral_val_upper` for clarity.

## New features

* Bivariate Hermite estimators have been added with methods for estimating 
bivariate probability density functions and cumulative distribution functions 
along with Spearman's rank correlation coefficients.
* The bivariate estimators include methods to batch update or sequentially 
update.
* Methods are also provided to consistently merge bivariate hermite_estimators.
* Convenience methods have been added for calculating normalized Hermite 
functions, along with upper, lower and full domain integrals of the 
normalized Hermite functions. 

## Documentation improvements

* The vignette `hermiter`, namely `vignette("hermiter")` has been extended to 
included examples pertaining to the bivariate Hermite series based estimators.

## Minor improvements and bug fixes
  
* The method for merging univariate Hermite series based estimators has been
improved, yielding greater accuracy when the hermite_estimators are 
standardized.
* The method for estimating quantiles with the univariate Hermite series based
estimator has been improved and is now consistent with the estimator in the
literature.
* Added further error trapping and other minor enhancements.

## Test environments
* local R installation, Windows 10, R 4.0.3
* local Linux Mint 20 installation, R 4.0.3 with valgrind
* Ubuntu 16.04 (on travis-ci and r-hub), R 4.0.3
* Debian Linux (r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (r-hub)
* Fedora Linux (r-hub)
* Oracle Solaris 10, x86, 32 bit, R-release (r-hub)
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (r-hub)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

Days since last update: 1

Explanation: Fixed minor issues in c++ source code causing the build to fail on 
Solaris along with ASAN/valgrind issues as brought to the package maintainer's 
attention by Prof Ripley.

