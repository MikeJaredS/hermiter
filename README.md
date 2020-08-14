# hermiter

## What does hermiter do?

`hermiter` facilitates the estimation of the full probability density function, cumulative distribution function and quantile function using Hermite series based estimators. The package is applicable to streaming, batch and grouped data.

These estimators are particularly useful in the sequential setting (both stationary and nonstationary data streams). In addition, they are useful in efficient, one-pass batch estimation which is particularly relevant in the context of large data sets. Finally, the Hermite series based estimators are applicable in decentralized (distributed) settings in that estimators formed on subsets of the data can be consistently combined. The Hermite series based estimators have the distinct advantage of being able to estimate the full density function, distribution function and quantile function in an online manner.The theoretical and empirical properties of these estimators for distribution function and quantile estimation have been studied in,

* [Stephanou, Michael, Melvin Varughese, and Iain Macdonald. "Sequential quantiles via Hermite series density estimation." Electronic Journal of Statistics 11.1 (2017): 570-607.](https://projecteuclid.org/euclid.ejs/1488531636) 
* [Stephanou, Michael, and Melvin Varughese. "On the properties of hermite series based distribution function estimators." Metrika (2020): 1-25.](https://link.springer.com/article/10.1007/s00184-020-00785-z)

## Features

* fast batch estimation of pdf, cdf and quantile function
* consistent combining of estimates on different subsets of a larger data set
* fast sequential estimation of pdf, cdf and quantile function on streaming data
* adaptive sequential estimation on non-stationary streams via exponential weighting

## Installation

```r
install.packages("hermiter")
```

## Usage

```r
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
hermite_est <- hermite_est %>% update_batch(observations)
x <- seq(-15,15,0.1)
pdf_est <- dens(hermite_est,x)
cdf_est <- cum_prob(hermite_est,x)
p <- seq(0.05,1,0.05)
quantile_est <- quant(hermite_est,p)
```
