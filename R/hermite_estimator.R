#' A class to sequentially estimate pdfs, cdfs and quantile functions
#'
#' This method constructs an S3 object with associated methods for univariate
#' nonparametric estimation of pdfs, cdfs and quantiles.
#'
#' The hermite_estimator class allows the sequential or one-pass batch
#' estimation of the full probability density function, cumulative distribution
#' function and quantile function. It is well suited to streaming data (both
#' stationary and non-stationary) and to efficient estimation in the context of
#' massive or distributed data sets. Indeed, estimators constructed on different
#' subsets of a distributed data set can be consistently combined.
#'
#' @author Michael Stephanou <michael.stephanou@gmail.com>
#'
#' @param N An integer between 0 and 100. The Hermite series based estimator
#' is truncated at N+1 terms.
#' @param standardize A boolean value. Determines whether the observations are
#' standardized, a transformation which often improves performance.
#' @param exp_weight_lambda A numerical value between 0 and 1. This parameter
#' controls the exponential weighting of the Hermite series based estimator.
#' If this parameter is NA, no exponential weighting is applied.
#' @return An S3 object of class hermite_estimator, with methods for
#' density function, distribution function and quantile function estimation.
#' @export
#' @examples
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
hermite_estimator <-
  function(N = 10,
           standardize = FALSE,
           exp_weight_lambda = NA) {
    if (!is.numeric(N)) {
      stop("N must be numeric.")
    }
    if (N < 0 | N > 100) {
      stop("N must be >= 0 and N <= 100.")
    }
    if (!(standardize == TRUE | standardize == FALSE)) {
      stop("standardize can only take on values TRUE or FALSE.")
    }
    if (!is.na(exp_weight_lambda)) {
      if (!is.numeric(exp_weight_lambda)) {
        stop("exp_weight_lambda must be numeric.")
      }
      if (exp_weight_lambda <= 0 | exp_weight_lambda > 1) {
        stop("exp_weight_lambda must be a real number > 0 and <= 1.")
      }
    }
    this <-
      list(
        N_param = N,
        coeff_vec = rep(0, N + 1),
        num_obs = 0,
        standardize_obs = standardize,
        running_mean = 0,
        running_variance = 0,
        exp_weight = exp_weight_lambda,
        normalization_hermite_vec = c()
      )
    this$normalization_hermite_vec <-
      hermite_normalization(this$N_param)
    class(this) <- append(class(this), "hermite_estimator")
    return(this)
  }

#' Returns Hermite series expansion coefficients
#'
#' @param this A hermite_estimator object. 
#' @return Numeric vector of length N+1
#' @export
#' @examples
#' hermite_est_1 <- hermite_estimator(N = 10, standardize = FALSE)
#' hermite_est_1 <- update_batch(hermite_est_1, rnorm(30))
#' coefficient_vec <- get_coefficients(hermite_est_1)
get_coefficients <- function(this) {
  UseMethod("get_coefficients", this)
}

#' @export
get_coefficients.hermite_estimator <-
  function(this) {
    return(this$coeff_vec)
}


#' Internal method to consistently combine the number of observations, means and 
#' variances of two Hermite estimators
#' 
#' The algorithm to combine the variances consistently comes from
#' Schubert, Erich, and Michael Gertz. "Numerically stable parallel computation 
#' of (co-) variance." Proceedings of the 30th International Conference on 
#' Scientific and Statistical Database Management. 2018.
#'
#' @param hermite_estimator1 A hermite_estimator object.
#' @param hermite_estimator2 A hermite_estimator object.
#' @return An object of class hermite_estimator.
combine_moments_and_count <- function(hermite_estimator1,hermite_estimator2){
  num_obs_1 <- hermite_estimator1$num_obs
  num_obs_2 <- hermite_estimator2$num_obs
  hermite_combined <- hermite_estimator(hermite_estimator1$N_param, 
                                        hermite_estimator1$standardize_obs)
  hermite_combined$num_obs <- num_obs_1 + num_obs_2
  hermite_combined$running_mean <- (num_obs_1*hermite_estimator1$running_mean +
                num_obs_2*hermite_estimator2$running_mean)/(num_obs_1+num_obs_2)
  hermite_combined$running_variance <- (hermite_estimator1$running_variance +
      hermite_estimator2$running_variance) + ((num_obs_1 * num_obs_2) / 
     (num_obs_1 + num_obs_2)) *(hermite_estimator1$running_mean 
                                - hermite_estimator2$running_mean)^2
  return(hermite_combined)
}

#' Internal method to combine a list of standardized Hermite estimators with 
#' improved accuracy
#'
#'
#' @param hermite_estimators A list of hermite_estimator objects.
#' @return An object of class hermite_estimator.
combine_standardized_helper <- function(hermite_estimators) {
  N <- hermite_estimators[[1]]$N_param
  hermite_estimator_combined <- base::Reduce(f=combine_moments_and_count, 
                                             x = hermite_estimators)
  integrand <- function(t,hermite_est_current, hermite_estimator_combined,
                        current_k){
    normalization_hermite <- hermite_est_current$normalization_hermite_vec
    t <- sqrt(2) * t
    original_sd <- sqrt(hermite_est_current$running_variance / 
                          (hermite_est_current$num_obs-1))
    original_mean <- hermite_est_current$running_mean
    new_sd <- sqrt(hermite_estimator_combined$running_variance / 
                     (hermite_estimator_combined$num_obs-1))
    new_mean <- hermite_estimator_combined$running_mean
    herm_mod <- hermite_polynomial(hermite_est_current$N_param, t) *
      normalization_hermite
    sqrt(2) * hermite_function(current_k,((t*original_sd +original_mean) - 
                                            new_mean)/new_sd, 
              normalization_hermite[1:current_k])[current_k, ] * 
              as.vector(crossprod(herm_mod, hermite_est_current$coeff_vec))
  }
  hermite_estimator_combined$coeff_vec <-
    sapply(1:(N+1),FUN=function(k){sum(sapply(hermite_estimators, 
          FUN=function(x){(x$num_obs / hermite_estimator_combined$num_obs) *
      gauss_hermite_quad_100(function(t){integrand(t, x, 
                                       hermite_estimator_combined, k)})}))})
  return(hermite_estimator_combined)
}

#' Combines two Hermite estimators
#'
#' This method allows a pair of Hermite based estimators of class
#' hermite_estimator to be consistently combined.
#'
#' Note that the N and standardize arguments must be the same for the two
#' estimators in order to combine them. In addition, note that exponentially
#' weighted estimators cannot be combined. If the Hermite estimators are not
#' standardized, the combined estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param this A hermite_estimator object. The first Hermite series based
#' estimator.
#' @param hermite_estimator_other A hermite_estimator object. The second Hermite
#' series based estimator.
#' @return An object of class hermite_estimator.
#' @export
#' @examples
#' hermite_est_1 <- hermite_estimator(N = 10, standardize = FALSE)
#' hermite_est_1 <- update_batch(hermite_est_1, rnorm(30))
#' hermite_est_2 <- hermite_estimator(N = 10, standardize = FALSE)
#' hermite_est_2 <- update_batch(hermite_est_2, rnorm(30))
#' hermite_combined <- combine_pair(hermite_est_1, hermite_est_2)
combine_pair <- function(this, hermite_estimator_other) {
  UseMethod("combine_pair", this)
}

#' @export
combine_pair.hermite_estimator <-
  function(this, hermite_estimator_other) {
    if (!is(hermite_estimator_other, "hermite_estimator")) {
      stop("combine_pair can only be applied to hermite_estimator objects.")
    }
    if (this$N_param != hermite_estimator_other$N_param) {
      stop("N must be equal to combine estimators.")
    }
    if (this$standardize_obs != hermite_estimator_other$standardize_obs) {
      stop("Standardization setting must be the same to combine estimators.")
    }
    if (!is.na(this$exp_weight) |
        !is.na(hermite_estimator_other$exp_weight)) {
      stop("Cannot combine exponentially weighted estimators.")
    }
    if (this$standardize_obs == FALSE) {
      hermite_estimator_combined <- combine_moments_and_count(this, 
                                                      hermite_estimator_other)
      hermite_estimator_combined$coeff_vec <-
        (
          this$coeff_vec * this$num_obs + hermite_estimator_other$coeff_vec 
          * hermite_estimator_other$num_obs
        ) / (this$num_obs + hermite_estimator_other$num_obs)
    } else {
      hermite_estimator_combined <- combine_standardized_helper(list(this, 
                                                     hermite_estimator_other))
    }
    return(hermite_estimator_combined)
  }

#' Combines a list of Hermite estimators
#'
#' This method allows a list of Hermite based estimators of class
#' hermite_estimator to be consistently combined.
#'
#' Note that the N and standardize arguments must be the same for all estimators
#' in order to combine them. In addition, note that exponentially weighted
#' estimators cannot be combined. If the Hermite estimators are not
#' standardized, the combined estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param hermite_estimators A list of hermite_estimator objects.
#' @return An object of class hermite_estimator.
#' @export
#' @examples
#' hermite_est_1 <- hermite_estimator(N = 10, standardize = FALSE)
#' hermite_est_1 <- update_batch(hermite_est_1, rnorm(30))
#' hermite_est_2 <- hermite_estimator(N = 10, standardize = FALSE)
#' hermite_est_2 <- update_batch(hermite_est_2, rnorm(30))
#' hermite_combined <- combine_hermite(list(hermite_est_1, hermite_est_2))
combine_hermite <- function(hermite_estimators) {
  UseMethod("combine_hermite", hermite_estimators)
}

#' @export
combine_hermite.list <- function(hermite_estimators) {
  if (length(hermite_estimators) == 1) {
    return(hermite_estimators[[1]])
  }
  if (hermite_estimators[[1]]$standardize_obs==FALSE){
   hermite_estimator_combined <- base::Reduce(combine_pair,hermite_estimators)
  } else {
   hermite_estimator_combined <- combine_standardized_helper(hermite_estimators)
  }
  return(hermite_estimator_combined)
}

#' Updates the Hermite series based estimator sequentially
#'
#' This method can be applied in sequential estimation settings.
#' 
#'
#' @param this A hermite_estimator object.
#' @param x A numeric value. An observation to be incorporated into the
#' estimator.
#' @return An object of class hermite_estimator.
#' @export
#' @examples
#' hermite_estimator <- hermite_estimator(N = 10, standardize = TRUE)
#' hermite_estimator <- update_sequential(hermite_estimator, x = 2)
update_sequential <- function(this, x) {
  UseMethod("update_sequential", this)
}

#' @export
update_sequential.hermite_estimator <- function(this, x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x) != 1) {
    stop("The sequential method is only 
         applicable to one observation at a time.")
  }
  this$num_obs <- this$num_obs + 1
  if (this$standardize_obs == TRUE) {
    if (is.na(this$exp_weight)) {
      processed_vec <-
        standardizeInputs(x,
                          this$num_obs,
                          this$running_mean,
                          this$running_variance)
      this$running_mean <- processed_vec[1]
      this$running_variance <- processed_vec[2]
      x <- processed_vec[3]
      if (this$num_obs < 2) {
        return(this)
      }
    } else {
      processed_vec <-
        standardizeInputsEW(x,
                            this$num_obs,
                            this$exp_weight,
                            this$running_mean,
                            this$running_variance)
      this$running_mean <- processed_vec[1]
      this$running_variance <- processed_vec[2]
      x <- processed_vec[3]
      if (this$num_obs < 2) {
        return(this)
      }
    }
  }
  h_k <-
    as.vector(hermite_function(this$N_param, x, this$normalization_hermite_vec))
  if (is.na(this$exp_weight)) {
    this$coeff_vec <-
      (this$coeff_vec * (this$num_obs - 1) + h_k) / this$num_obs
  } else {
    this$coeff_vec <-
      this$coeff_vec * (1 - this$exp_weight) + h_k * this$exp_weight
  }
  return(this)
}

#' Updates the Hermite series based estimator with a batch of data
#'
#' This method can be applied in one-pass batch estimation settings. This
#' method cannot be used with an exponentially weighted estimator.
#'
#' @param this A hermite_estimator object.
#' @param x A numeric vector. A vector of observations to be incorporated
#' into the estimator.
#' @return An object of class hermite_estimator.
#' @export
#' @examples
#' hermite_estimator <- hermite_estimator(N = 10, standardize = TRUE)
#' hermite_estimator <- update_batch(hermite_estimator, x = c(1, 2))
update_batch <- function(this, x) {
  UseMethod("update_batch", this)
}

#' @export
update_batch.hermite_estimator <- function(this, x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (!is.na(this$exp_weight)) {
    stop("The Hermite estimator cannot be exponentially weighted.")
  }
  this$num_obs <- this$num_obs + length(x)
  if (this$standardize_obs == TRUE) {
    this$running_mean <-
      (this$running_mean * (this$num_obs - length(x)) + length(x) * mean(x)) / 
      this$num_obs
    this$running_variance <-
      (
        this$running_variance * (this$num_obs - length(x)) 
        + length(x) * stats::var(x) * (length(x) - 1)
      ) / this$num_obs
    x <-
      (x - this$running_mean) / sqrt(this$running_variance / (this$num_obs - 1))
  }
  h_k <-
    hermite_function(this$N_param, x, this$normalization_hermite_vec)
  this$coeff_vec <-
    (this$coeff_vec * (this$num_obs - length(x)) + rowSums(h_k)) / this$num_obs
  return(this)
}

#' Standardize a vector of observations x
#'
#' This helper method standardizes the observations x using the online mean and
#' online standard deviation contained in the hermite_estimator object (this).
#'
#' @param this A hermite_estimator object.
#' @param x A numeric vector of observations.
#' @return An object of class hermite_estimator.
#' @export
standardize_value <- function(this, x) {
  UseMethod("standardize_value", this)
}

#' @export
standardize_value.hermite_estimator <- function(this, x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (is.na(this$exp_weight)) {
    running_std <- sqrt(this$running_variance / (this$num_obs - 1))
  } else {
    running_std <- sqrt(this$running_variance)
  }
  return((x - this$running_mean) / running_std)
}

#' Estimates the cumulative probability for a vector of x values
#'
#' This method calculates the cumulative probability values at a vector of
#' x values using the hermite_estimator object (this).
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param this A hermite_estimator object.
#' @param x A numeric vector. Values at which to estimate the cumulative
#' probability
#' @param clipped A boolean value. This value determines whether cumulative
#' probabilities are clipped to lie within the range [0,1].
#' @return A numeric vector of cumulative probability values.
#' @export
#' @examples
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, rnorm(30))
#' cdf_est <- cum_prob(hermite_est, c(0, 0.5, 1))
cum_prob <- function(this, x, clipped) {
  UseMethod("cum_prob", this)
}

#' @export
cum_prob.hermite_estimator <- function(this, x, clipped = FALSE) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (this$num_obs < 2) {
    return(NA)
  }
  if (this$standardize_obs == TRUE) {
    x <- standardize_value(this, x)
  }
  h_k <-
    hermite_function(this$N_param, x, this$normalization_hermite_vec)
  integrals_hermite <- hermite_integral_val(this$N_param, x, h_k)
  cdf_val <- crossprod(integrals_hermite, this$coeff_vec)
  if (clipped == TRUE) {
    cdf_val <- pmin(pmax(cdf_val, 1e-08), 1)
  }
  return(as.vector(cdf_val))
}

#' Estimates the probability density for a vector of x values
#'
#' This method calculates the probability density values at a vector of
#' x values using the hermite_estimator object (this).
#'
#' The object must be updated with observations prior to the use of the method.
#'
#' @param this A hermite_estimator object.
#' @param x A numeric vector. Values at which to estimate the probability
#' density.
#' @param clipped A boolean value. This value determines whether
#' probability densities are clipped to be bigger than zero.
#' @return A numeric vector of probability density values.
#' @export
#' @examples
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, rnorm(30))
#' pdf_est <- dens(hermite_est, c(0, 0.5, 1))
dens <- function(this, x, clipped) {
  UseMethod("dens", this)
}

#' @export
dens.hermite_estimator <- function(this, x, clipped = FALSE) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (this$num_obs < 2) {
    return(NA)
  }
  factor <- 1
  if (this$standardize_obs == TRUE) {
    if (is.na(this$exp_weight)) {
      running_std <- sqrt(this$running_variance / (this$num_obs - 1))
    } else {
      running_std <- sqrt(this$running_variance)
    }
    x <- (x - this$running_mean) / running_std
    factor <- 1 / running_std
  }
  h_k <-
    hermite_function(this$N_param, x, this$normalization_hermite_vec)
  pdf_val <- crossprod(h_k, this$coeff_vec) * factor
  if (clipped == TRUE) {
    pdf_val <- pmax(pdf_val, 1e-08)
  }
  return(as.vector(pdf_val))
}

#' Estimates the quantile at a single probability value
#'
#' This helper method is intended for internal use by the hermite_estimator
#' class.
#'
#' @param this A hermite_estimator object.
#' @param p A numeric value. The probability at which to calculate the
#' associated quantile.
#' @return A numeric value. The value of the quantile associated with p.
#' @export
quantile_helper <- function(this, p) {
  UseMethod("quantile_helper", this)
}

#' @export
quantile_helper.hermite_estimator <- function(this, p) {
  if (!is.numeric(p)) {
    stop("p must be numeric.")
  }
  if (p < 0 | p > 1) {
    return(NA)
  }
  h_k <-
    hermite_function(this$N_param, x=0, this$normalization_hermite_vec)
  p_lower <- this$coeff_vec %*% hermite_integral_val(this$N_param,x=0,h_k)
  p_upper <- 1-as.numeric(this$coeff_vec %*% 
                      hermite_integral_val_quantile_adap(this$N_param,x=0,h_k))
  if (is.na(p_lower) | is.na(p_upper)){
    return(NA)
  }
  if (p_upper < p_lower){
    x_lower <- tryCatch({
      stats::uniroot(
        f = function(x) {
          this$coeff_vec %*% hermite_int_lower(this$N_param,x) - p_upper
        },
        interval = c(-100, 100)
      )$root
    },
    error = function(e) {NA})
    x_upper <- tryCatch({
      stats::uniroot(
        f = function(x) {
          1-as.numeric(this$coeff_vec %*% 
                         hermite_int_upper(this$N_param,x)) - p_lower
        },
        interval = c(-100, 100)
      )$root
    },
    error = function(e) {NA})
  } else if (p_upper > p_lower) {
    x_lower <- -1e-6
    x_upper <- 1e-6
  } else if (p_upper == p_lower){
    x_lower <- 0
    x_upper <- 0
  }
  if (is.na(x_lower) | is.na(x_upper)){
    return(NA)
  }
  est <- tryCatch({
    stats::uniroot(
      f = function(x) {
        lower_idx <- which(x < x_lower)
        upper_idx <- which(x > x_upper)
        ambig_idx <- which(x >= x_lower & x <= x_upper)
        res <- rep(NA,length(x))
        if (length(lower_idx)>0){
          res[lower_idx] <- crossprod(hermite_int_lower(this$N_param, 
                                        x[lower_idx]), this$coeff_vec) - p
        }
        if (length(upper_idx)>0){
          res[upper_idx] <- 1 - 
            crossprod(hermite_int_upper(this$N_param, 
                                        x[upper_idx]), this$coeff_vec) - p
        }
        if (length(ambig_idx)>0){
          if (p_upper < p_lower){
            res[ambig_idx] <- (p_upper + (x[ambig_idx]-x_lower) * 
                           as.numeric((p_lower-p_upper)/(x_upper-x_lower))) - p
          }
          else if (p_upper > p_lower){
            res[ambig_idx] <- (p_lower + (x[ambig_idx]-x_lower) * 
                                 (p_upper-p_lower)/(x_upper-x_lower)) - p
          } else if (p_upper == p_lower){
            res[ambig_idx] <- p_upper
          }
        }
        return(res)
      },
      interval = c(-100, 100)
    )$root
  }, error = function(e) {NA})
  if (is.na(this$exp_weight)) {
    est <-
      est * sqrt(this$running_variance / (this$num_obs - 1)) + this$running_mean
  } else {
    est <- est * sqrt(this$running_variance) + this$running_mean
  }
  return(est)
}

#' Estimates the quantiles at a vector of probability values
#' 
#' This method utilizes the estimator (13) in paper Stephanou, Michael, 
#' Varughese, Melvin and Iain Macdonald. "Sequential quantiles via Hermite 
#' series density estimation." Electronic Journal of Statistics 11.1 (2017): 
#' 570-607 <doi:10.1214/17-EJS1245>, with some modifications to improve the 
#' stability of numerical root finding. 
#'
#' @param this A hermite_estimator object.
#' @param p A numeric vector. A vector of probability values.
#' @return A numeric vector. The vector of quantile values associated with the
#' probabilities p.
#' @export
#' @examples
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, rnorm(30))
#' quant_est <- quant(hermite_est, c(0.25, 0.5, 0.75))
quant <- function(this, p) {
  UseMethod("quant", this)
}

#' @export
quant.hermite_estimator <- function(this, p) {
  if (!is.numeric(p)) {
    stop("p must be numeric.")
  }
  if (this$standardize_obs != TRUE) {
    stop("Quantile estimation requires standardization to be true.")
  }
  if (length(p) < 1) {
    stop("At least one quantile must be specified.")
  }
  if (this$num_obs < 2) {
    return(NA)
  }
  result <- rep(0, length(p))
  for (idx in seq_along(p)) {
    result[idx] <- quantile_helper(this, p[idx])
  }
  return(result)
}
