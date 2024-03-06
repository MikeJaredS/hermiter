# hermiter v2.3.1

## Minor improvements and bug fixes

* Updated parallelized C++ methods to be more thread-safe and included 
improved bound checking.
* Updated citation for hermiter article.

## Test environments
* local R installation, Mac OS 14.0, R 4.3.2
* Ubuntu Linux 20.04.1 LTS, R-release, GCC (r-hub)
* Fedora Linux, R-devel, clang, gfortran (r-hub)
* win-builder (devel, release, oldrel)

## R CMD check results

There was 1 NOTE across all environments:

GNU make is a SystemRequirements.
    
* This appears to be related to using RcppParallel and does not seem 
problematic. GNU make has been added to SystemRequirements in the Description
file.

* In addition, on Ubuntu Linux, the package 
appears to be slightly larger than 5MB, around 6.3 MB. This does not
appear to occur on any of the other environments however (i.e. Fedora Linux, 
Windows)

* Finally, there is the following note on Ubuntu and Fedora Linux:

* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
Skipping checking math rendering: package 'V8' unavailable

This seems build system specific and not directly tied to the package.

* There is also a PREPERROR for Debian Linux, R-devel, GCC ASAN/UBSAN on rhub 

We believe this is an issue on the build server and not package specific.
