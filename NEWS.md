# hermiter v2.0.2

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
* Vectorized the univariate quantile estimation method.
* Added further error trapping and other minor enhancements (also for C++
routines).

# hermiter v1.0.0

## INITIAL RELEASE