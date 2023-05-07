# hermiter v2.3.0

## Additions and improvements

* Enhanced the update_sequential method to incorporate a single or multiple
new observations.
* Added density generic function which outputs an object with associated print
and plot generics.
* Added quantile generic function for convenience.
* Added hcdf function which outputs an object with associated print, 
plot and summary generics.
* Added median and IQR convenience functions.
* Added a wrapper around the stats::cor function with two new methods, namely
"hermite.spearman" and "hermite.kendall".

## Minor improvements and bug fixes

* Additional test cases have been added.
* Updated the vignette.
* Fixed the SystemRequirements: C++11 note.


## Test environments
* local R installation, Windows 10, R 4.3.0
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

