# hermiter v2.2.0

## Breaking changes

* The interface of `hermiter` has been simplified. The `update_batch` method 
has been removed in favor of providing the ability to initialize the 
`hermite_estimator` with an initial batch of observations. Several internal 
methods are no longer exported in the interests of simplicity.
* The default values of N have been optimized for different settings. For 
univariate, non-exponentially weighted estimators, the default is now N = 50. 
For univariate, exponentially weighted estimators, the default is now N = 20. 
For bivariate, non-exponentially weighted estimators, the default is now N = 30.
Finally, For bivariate, exponentially weighted estimators, the default is now 
N = 20.

## Major enhancements

* Parallel implementation of batch updating using RcppParallel provides 
significant performance improvements on multicore systems.

## Minor improvements and bug fixes

* Updated citation information.
* Additional test cases have been added.
* Bug fixes for series acceleration algorithm.


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

GNU make is a SystemRequirements.
    
    
* This appears to be related to using RcppParallel and does not seem 
problematic. GNU make has been added to SystemRequirements in the Description
file.

