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
* Debian Linux, R-release, GCC (r-hub)
* Fedora Linux, R-devel, clang, gfortran (r-hub)
* macOS 10.13.6 High Sierra, R-release, brew (r-hub)
* Windows Server 2022, R-devel, 64 bit (r-hub)
* win-builder (devel, release. oldrel)

## R CMD check results

There was 1 NOTE across all environments:

GNU make is a SystemRequirements.
    
* This appears to be related to using RcppParallel and does not seem 
problematic. GNU make has been added to SystemRequirements in the Description
file.

* In addition, on Debian and Ubuntu Linux, the package 
appears to be slightly larger than 5MB, around 5.2- 6 MB. This does not
appear to occur on any of the other environments however (i.e. Fedora Linux, 
Windows, MacOS)

* Finally, there is the following note on Ubuntu and Fedora Linux:

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable

This seems build system specific and not directly tied to the package.

* There is also a PREPERROR for Debian Linux, R-devel, GCC ASAN/UBSAN on rhub 
- we believe this is an issue on the build server and not package specific.
