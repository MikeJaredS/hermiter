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
significant performance improvements on multicore systems. Note that this 
can be disabled by using options(hermiter.parallel = FALSE).

## Minor improvements and bug fixes

* Updated citation information.
* Additional test cases have been added.
* Bug fixes for series acceleration algorithm.


## Test environments
* local R installation, Windows 10, R 4.2.2
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (r-hub)
* Debian Linux, R-devel, GCC ASAN/UBSAN (r-hub)
* Fedora Linux, R-devel, clang, gfortran (r-hub)
* macOS 10.13.6 High Sierra, R-release, brew
* Windows Server 2022, R-devel, 64 bit (r-hub)
* win-builder (devel, release. oldrel)

## R CMD check results

There was 1 NOTE across environments:

GNU make is a SystemRequirements.
    
* This appears to be related to using RcppParallel and does not seem 
problematic. GNU make has been added to SystemRequirements in the Description
file.

* In addition, on Ubuntu 20.04.1 LTS, the package 
appears to be slightly larger than 5MB, around 5.9 MB. This does not
appear to occur on any of the other environments however (i.e. Fedora Linux, 
Windows, MacOS)

