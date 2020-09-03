# hermiter v1.0.0

## INITIAL RELEASE

## Test environments
* local R installation, Windows 10, R 4.0.2
* Ubuntu 16.04 (on travis-ci and r-hub), R 4.0.2
* Debian Linux (r-hub)
* Fedora Linux (r-hub)
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTES:

New submission

Possibly mis-spelled words in DESCRIPTION:
  Iain (16:55)
  Macdonald (16:60)
  Metrika (17:123)
  Stephanou (16:13, 17:3)
  Varughese (16:33, 17:26)

Explanation: The above are not mis-spelled words but are in fact proper nouns.
This note does not appear in the local R installation checks or the r-hub checks
on Debian.

* checking for future file timestamps ... NOTE
  unable to verify current time
  
Explanation: It appears the worldclockapi service is down. This note is 
only present in the local R installation and Ubuntu 16.04 checks. It is not 
present in win-builder and certain r-hub checks (Debian Linux and Fedora Linux).
