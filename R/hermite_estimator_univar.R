#' A class to sequentially estimate univariate pdfs, cdfs and quantile functions
#'
#' This method constructs an S3 object with associated methods for univariate
#' nonparametric estimation of pdfs, cdfs and quantiles.
#'
#' The hermite_estimator_univar class allows the sequential or one-pass batch
#' estimation of the full probability density function, cumulative distribution
#' function and quantile function. It is well suited to streaming data (both
#' stationary and non-stationary) and to efficient estimation in the context of
#' massive or distributed data sets. Indeed, estimators constructed on different
#' subsets of a distributed data set can be consistently merged.
#'
#' @author Michael Stephanou <michael.stephanou@gmail.com>
#'
#' @param N An integer between 0 and 75. The upper bound has been chosen
#' as a value that yields an estimator that is reasonably fast and that remains 
#' robust to numerical issues. The Hermite series based estimator is truncated 
#' at N+1 terms.
#' @param standardize A boolean value. Determines whether the observations are
#' standardized, a transformation which often improves performance.
#' @param exp_weight_lambda A numerical value between 0 and 1. This parameter
#' controls the exponential weighting of the Hermite series based estimator.
#' If this parameter is NA, no exponential weighting is applied.
#' @param observations A numeric vector. A vector of observations to be 
#' incorporated into the estimator.
#' @return An S3 object of class hermite_estimator_univar, with methods for
#' density function, distribution function and quantile function estimation.
hermite_estimator_univar <-
  function(N = 50,
           standardize = TRUE,
           exp_weight_lambda = NA,
           observations = c()) {
    h_est_obj <-
      list(
        N_param = N,
        coeff_vec = rep(0, N + 1),
        num_obs = 0,
        standardize_obs = standardize,
        running_mean = 0,
        running_variance = 0,
        exp_weight = exp_weight_lambda,
        normalization_hermite_vec = c(),
        h_int_full_domain_serialized
      )
    h_est_obj$normalization_hermite_vec <- 
      h_norm_serialized[1:(h_est_obj$N_param+1)]
    h_est_obj$h_int_full_domain_serialized <- 
      h_int_full_domain_serialized[1:(h_est_obj$N_param+1),,drop = FALSE]
    class(h_est_obj) <- c("hermite_estimator_univar", "list")
    if (length(observations) > 0) {
      if (!is.numeric(observations)) {
        stop("observations must be numeric.")
      }
      if (any(is.na(observations)) | any(!is.finite(observations))){
        stop("The batch initialization method is only 
         applicable to finite, non NaN, non NA values.")
      }
      h_est_obj <- initialize_batch_univar(h_est_obj, observations)
    }
    return(h_est_obj)
  }

#' Internal method to consistently merge the number of observations, means and 
#' variances of two Hermite estimators
#' 
#' The algorithm to merge the variances consistently comes from
#' Schubert, Erich, and Michael Gertz. "Numerically stable parallel computation 
#' of (co-) variance." Proceedings of the 30th International Conference on 
#' Scientific and Statistical Database Management. 2018.
#'
#' @param hermite_estimator1 A hermite_estimator_univar object.
#' @param hermite_estimator2 A hermite_estimator_univar object.
#' @return An object of class hermite_estimator_univar.
merge_moments_and_count_univar <- function(hermite_estimator1, 
                                             hermite_estimator2){
  num_obs_1 <- hermite_estimator1$num_obs
  num_obs_2 <- hermite_estimator2$num_obs
  hermite_merged <- hermite_estimator_univar(hermite_estimator1$N_param, 
                                        hermite_estimator1$standardize_obs)
  hermite_merged$num_obs <- num_obs_1 + num_obs_2
  hermite_merged$running_mean <- (num_obs_1*hermite_estimator1$running_mean +
                num_obs_2*hermite_estimator2$running_mean)/(num_obs_1+num_obs_2)
  hermite_merged$running_variance <- (hermite_estimator1$running_variance +
      hermite_estimator2$running_variance) + ((num_obs_1 * num_obs_2) / 
     (num_obs_1 + num_obs_2)) *(hermite_estimator1$running_mean 
                                - hermite_estimator2$running_mean)^2
  return(hermite_merged)
}

#' Internal method to merge a list of standardized Hermite estimators
#'
#'
#' @param hermite_estimators A list of hermite_estimator_univar objects.
#' @return An object of class hermite_estimator_univar.
merge_standardized_helper_univar <- function(hermite_estimators) {
  all_N <- lapply(hermite_estimators, FUN =
                          function(x){return(x$N_param)})
  if (length(unique(all_N)) >1) {
    stop("List must contain Hermite estimators with a consistent N")
  }
  N <- hermite_estimators[[1]]$N_param
  hermite_estimator_merged <- base::Reduce(f=merge_moments_and_count_univar, 
                                             x = hermite_estimators)
  hermite_estimator_merged$coeff_vec <-
    vapply(1:(N+1),FUN=function(k){sum(vapply(hermite_estimators, 
          FUN=function(x){(x$num_obs / hermite_estimator_merged$num_obs) *
      gauss_hermite_quad_100(function(t){integrand_coeff_univar(t, x, 
                           hermite_estimator_merged, k)})},
      FUN.VALUE=numeric(1)))}, FUN.VALUE=numeric(1))
  return(hermite_estimator_merged)
}

#' Merges two Hermite estimators
#'
#' This method allows a pair of Hermite based estimators of class
#' hermite_estimator_univar to be consistently merged.
#'
#' Note that the N and standardize arguments must be the same for the two
#' estimators in order to merge them. In addition, note that exponentially
#' weighted estimators cannot be merged. If the Hermite estimators are not
#' standardized, the merged estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator_univar inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param h_est_obj A hermite_estimator_univar object. The first Hermite series 
#' based estimator.
#' @param hermite_estimator_other A hermite_estimator_univar object. The 
#' second Hermite series based estimator.
#' @return An object of class hermite_estimator_univar.
merge_pair.hermite_estimator_univar <-
  function(h_est_obj, hermite_estimator_other) {
    if (!is(hermite_estimator_other, "hermite_estimator_univar")) {
      stop("merge_pair.hermite_estimator_univar can only be applied to 
           hermite_estimator_univar objects.")
    }
    if (h_est_obj$N_param != hermite_estimator_other$N_param) {
      stop("N must be equal to merge estimators.")
    }
    if (h_est_obj$standardize_obs != hermite_estimator_other$standardize_obs) {
      stop("Standardization setting must be the same to merge estimators.")
    }
    if (!is.na(h_est_obj$exp_weight) |
        !is.na(hermite_estimator_other$exp_weight)) {
      stop("Cannot merge exponentially weighted estimators.")
    }
    if (h_est_obj$standardize_obs == FALSE) {
      hermite_estimator_merged <- merge_moments_and_count_univar(h_est_obj, 
                                                      hermite_estimator_other)
      hermite_estimator_merged$coeff_vec <-
        (
          h_est_obj$coeff_vec * h_est_obj$num_obs + 
            hermite_estimator_other$coeff_vec * hermite_estimator_other$num_obs
        ) / (h_est_obj$num_obs + hermite_estimator_other$num_obs)
    } else {
      hermite_estimator_merged <- 
        merge_standardized_helper_univar(list(h_est_obj, 
                                              hermite_estimator_other))
    }
    return(hermite_estimator_merged)
  }

#' Merges a list of Hermite estimators
#'
#' This method allows a list of Hermite based estimators of class
#' hermite_estimator_univar to be consistently merged.
#'
#' Note that the N and standardize arguments must be the same for all estimators
#' in order to merge them. In addition, note that exponentially weighted
#' estimators cannot be merged. If the Hermite estimators are not
#' standardized, the merged estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator_univar inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param hermite_estimators A list of hermite_estimator_univar objects.
#' @return An object of class hermite_estimator_univar.
merge_hermite_univar <- function(hermite_estimators) {
  if (length(hermite_estimators) == 1) {
    return(hermite_estimators[[1]])
  }
  if (hermite_estimators[[1]]$standardize_obs==FALSE){
   hermite_estimator_merged <- base::Reduce(merge_pair,hermite_estimators)
  } else {
   hermite_estimator_merged <- 
     merge_standardized_helper_univar(hermite_estimators)
  }
  return(hermite_estimator_merged)
}

# Helper function for univariate update_sequential function.
update_sequential_univar_helper <- function(h_est_obj,x){
  h_est_obj$num_obs <- h_est_obj$num_obs + 1
  if (h_est_obj$standardize_obs == TRUE) {
    if (is.na(h_est_obj$exp_weight)) {
      prev_running_mean <- h_est_obj$running_mean
      h_est_obj$running_mean <-  (h_est_obj$running_mean * 
                                    (h_est_obj$num_obs - 1) + x) / 
        h_est_obj$num_obs
      if (h_est_obj$num_obs < 2) {
        return(h_est_obj)
      }
      h_est_obj$running_variance <- h_est_obj$running_variance + 
        (x - prev_running_mean) * (x - h_est_obj$running_mean)
      x <- (x - h_est_obj$running_mean) /  sqrt(h_est_obj$running_variance /
                                                  (h_est_obj$num_obs - 1))
    } else {
      if (h_est_obj$num_obs < 2){
        h_est_obj$running_mean <- x
        h_est_obj$running_variance <- 1
        return(h_est_obj)
      }
      h_est_obj$running_mean <-  (1 - h_est_obj$exp_weight) * 
        h_est_obj$running_mean + h_est_obj$exp_weight * x
      h_est_obj$running_variance <- (1 - h_est_obj$exp_weight) * 
        h_est_obj$running_variance + h_est_obj$exp_weight * 
        (x - h_est_obj$running_mean)^2
      x <- (x - h_est_obj$running_mean) / sqrt(h_est_obj$running_variance)
    }
  }
  if (is.na(x)){
    return(h_est_obj)
  }
  h_k <-
    as.vector(hermite_function(h_est_obj$N_param, x))
  if (is.na(h_est_obj$exp_weight)) {
    h_est_obj$coeff_vec <-
      (h_est_obj$coeff_vec * (h_est_obj$num_obs - 1) + h_k) / 
      h_est_obj$num_obs
  } else {
    h_est_obj$coeff_vec <-
      h_est_obj$coeff_vec * (1 - h_est_obj$exp_weight) + h_k * 
      h_est_obj$exp_weight
  }
  return(h_est_obj)
}

#' Updates the Hermite series based estimator sequentially
#'
#' This method can be applied in sequential estimation settings.
#' 
#'
#' @param h_est_obj A hermite_estimator_univar object.
#' @param x A numeric vector. A vector of observations to be incorporated into 
#' the estimator.
#' @return An object of class hermite_estimator_univar.
update_sequential.hermite_estimator_univar <- function(h_est_obj, x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x) < 1) {
    stop("x must contain at least one observation.")
  }
  if (any(is.na(x) | !is.finite(x))){
    stop("The sequential method is only 
         applicable to finite, non NaN, non NA values.")
  }
  for (idx in seq_along(x)) {
    h_est_obj <- update_sequential_univar_helper(h_est_obj,
                                     x[idx])
  }
  return(h_est_obj)
}

#' Initializes the Hermite series based estimator with a batch of data
#'
#' @param h_est_obj A hermite_estimator_univar object.
#' @param x A numeric vector. A vector of observations to be incorporated
#' into the estimator.
#' @return An object of class hermite_estimator_univar.
initialize_batch_univar <- function(h_est_obj, x) {
  h_est_obj$num_obs <- length(x)
  if (h_est_obj$standardize_obs == TRUE) {
    h_est_obj$running_mean <- mean(x)
    h_est_obj$running_variance <- stats::var(x) * (length(x) - 1)
    x <-
      (x - h_est_obj$running_mean) / sqrt(h_est_obj$running_variance / 
                                            (h_est_obj$num_obs - 1))
  }
  h_est_obj$coeff_vec <- hermite_function_sum_N(h_est_obj$N_param, x) /
    h_est_obj$num_obs
  return(h_est_obj)
}

calculate_running_std.hermite_estimator_univar<- function(h_est_obj){
  if (is.na(h_est_obj$exp_weight)) {
    running_std <- sqrt(h_est_obj$running_variance / (h_est_obj$num_obs - 1))
  } else {
    running_std <- sqrt(h_est_obj$running_variance)
  }
  return(running_std)
}

#' Estimates the probability density for a vector of x values
#'
#' This method calculates the probability density values at a vector of
#' x values using the hermite_estimator_univar object (h_est_obj).
#'
#' The object must be updated with observations prior to the use of the method.
#'
#' @param h_est_obj A hermite_estimator_univar object.
#' @param x A numeric vector. Values at which to estimate the probability
#' density.
#' @param clipped A boolean value. This value determines whether
#' probability densities are clipped to be bigger than zero.
#' @param accelerate_series A boolean value. This value determines whether
#' Hermite series acceleration is applied.
#' @return A numeric vector of probability density values.
dens.hermite_estimator_univar <- function(h_est_obj, x, clipped = FALSE, 
                                          accelerate_series = TRUE) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x) < 1) {
    stop("x must contain at least one value.")
  }
  if (h_est_obj$num_obs < 2) {
    return(rep(NA, length(x)))
  }
  factor <- 1
  if (h_est_obj$standardize_obs == TRUE) {
    running_std <- calculate_running_std(h_est_obj)
    x <- (x - h_est_obj$running_mean) / running_std
    factor <- 1 / running_std
  }
  h_k <-
    hermite_function(h_est_obj$N_param, x)
  pdf_val <- series_calculate(h_k, h_est_obj$coeff_vec, accelerate_series = 
                                accelerate_series) * factor
  if (clipped == TRUE) {
    pdf_val <- pmax(pdf_val, 1e-08)
  }
  return(as.vector(pdf_val))
}

#' Creates an object summarizing the PDF with associated generic methods 
#' print and plot.
#'
#' The hermite_estimator_univar, x must be updated with observations prior to 
#' the use of the method.
#'
#' @param x A hermite_estimator_univar object.
#' @param x_lower A numeric value. This value determines the lower limit of 
#' x values at which to evaluate the density.
#' @param x_upper A numeric value. This value determines the upper limit of 
#' x values at which to evaluate the density.
#' @param ... Additional arguments for the dens function.
#' @return A hdensity_univar object whose underlying structure is a list 
#' containing the following components.
#' 
#' x: The points at which the density is calculated.
#' density_vals: The density values at the points x.
#' num_obs: The number of observations used to form the Hermite density 
#' estimates.
#' N: The number of terms N in the Hermite series estimator.
#' @export
density.hermite_estimator_univar <- function(x, x_lower = NA, x_upper = NA, 
                                             ...) {
  if (x$standardize_obs == FALSE){
    if (is.na(x_lower) | is.na(x_upper)){
      stop("For non-standardized hermite_estimator objects, a lower and an
           upper x limit for the PDF summary must be provided i.e. x_lower and 
           x_upper arguments must be provided.")
    }
  } 
  if (!is.na(x_lower) & !is.na(x_upper)){
    x_step <- (x_upper - x_lower)/99
    x_vals <- seq(x_lower, x_upper, by = x_step)
  } else {
    x_vals <- quant(x, p = seq(0.01,0.99,0.01))
  }
  density_vals <- dens(x, x_vals, ...)
  result <- list(x = x_vals, density_vals = density_vals, 
                 num_obs = x$num_obs, N = x$N_param)
  class(result) <- c("hdensity_univar", "list")
  return(result)
}

#' Estimates the cumulative probability for a vector of x values
#'
#' This method calculates the cumulative probability values at a vector of
#' x values using the hermite_estimator_univar object (h_est_obj).
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_univar object.
#' @param x A numeric vector. Values at which to estimate the cumulative
#' probability
#' @param clipped A boolean value. This value determines whether cumulative
#' probabilities are clipped to lie within the range [0,1].
#' @param accelerate_series A boolean value. This value determines whether
#' Hermite series acceleration is applied.
#' @return A numeric vector of cumulative probability values.
cum_prob.hermite_estimator_univar <- function(h_est_obj, x, clipped = FALSE, 
                                              accelerate_series = TRUE) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x) < 1) {
    stop("x must contain at least one value.")
  }
  if (h_est_obj$num_obs < 2) {
    return(rep(NA, length(x)))
  }
  if (h_est_obj$standardize_obs == TRUE) {
    if (h_est_obj$running_variance == 0) {
      return(ifelse(h_est_obj$running_mean <= x,1,0))
    }
    running_std <- calculate_running_std(h_est_obj)
    x <- (x - h_est_obj$running_mean) / running_std
  }
  h_k <-
    hermite_function(h_est_obj$N_param, x)
  integrals_hermite <- hermite_int_lower(h_est_obj$N_param, x, 
                                         hermite_function_matrix = h_k)
  cdf_val <- series_calculate(integrals_hermite, h_est_obj$coeff_vec, 
                              accelerate_series = accelerate_series)
  if (clipped == TRUE) {
    cdf_val <- pmin(pmax(cdf_val, 1e-08), 1)
  }
  return(as.vector(cdf_val))
}

#' Creates an object summarizing the CDF with associated generic methods print,
#' plot and summary.
#'
#' The hermite_estimator_univar object, h_est_obj must be updated with 
#' observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_univar object.
#' @param clipped A boolean value. This value determines whether cumulative
#' probabilities are clipped to lie within the range [0,1].
#' @param accelerate_series A boolean value. This value determines whether
#' Hermite series acceleration is applied.
#' @param x_lower A numeric value. This value determines the lower limit of 
#' x values at which to evaluate the CDF.
#' @param x_upper A numeric value. This value determines the upper limit of 
#' x values at which to evaluate the CDF.
#' @return A hcdf_univar object whose underlying structure is a list 
#' containing the following components.
#' 
#' x: The points at which the cumulative probability is calculated.
#' cum_prob_vals: The cumulative probability values at the points x.
#' num_obs: The number of observations used to form the Hermite cumulative 
#' probability estimates.
#' N: The number of terms N in the Hermite series estimator.
#' @export
hcdf.hermite_estimator_univar <- function(h_est_obj, clipped = FALSE, 
                                          accelerate_series = TRUE, x_lower=NA,
                                          x_upper=NA) {
  if (h_est_obj$standardize_obs == FALSE){
    if (is.na(x_lower) | is.na(x_upper)){
      stop("For non-standardized hermite_estimator objects, a lower and an
           upper x limit for the CDF summary must be provided i.e. x_lower and 
           x_upper arguments must be provided.")
    }
  } 
  if (!is.na(x_lower) & !is.na(x_upper)){
    x_step <- (x_upper - x_lower)/99
    x_vals <- seq(x_lower, x_upper, by = x_step)
  } else {
    x_vals <- quant(h_est_obj, p = seq(0.01,0.99,0.01))
  }
  cum_prob_vals <- cum_prob(h_est_obj, x_vals, 
                            clipped, accelerate_series)
  result <- list(x = x_vals, cum_prob_vals = cum_prob_vals, 
                 num_obs = h_est_obj$num_obs, N = h_est_obj$N_param)
  class(result) <- c("hcdf_univar", "list")
  return(result)
}

# Estimates the quantile at a vector of probability values using a vectorized
# implementation of the bisection search root finding algorithm.
#
# This helper method is intended for internal use by the 
# hermite_estimator_univar class.
quantile_helper_bisection <- function(h_est_obj, p_vec, accelerate_series) {
  f_est <- function(x,p) {
    lower_idx <- which(x < x_lower)
    upper_idx <- which(x > x_upper)
    ambig_idx <- which(x >= x_lower & x <= x_upper)
    res <- rep(NA,length(x))
    if (length(lower_idx)>0){
      res[lower_idx] <- series_calculate(hermite_int_lower(h_est_obj$N_param,
         x[lower_idx]), 
         h_est_obj$coeff_vec, accelerate_series) - p[lower_idx]
    }
    if (length(upper_idx)>0){
      res[upper_idx] <- 1 - 
        series_calculate(hermite_int_upper(h_est_obj$N_param, x[upper_idx]), 
         h_est_obj$coeff_vec, accelerate_series) - p[upper_idx]
    }
    if (length(ambig_idx)>0){
      if (p_upper < p_lower){
        res[ambig_idx] <- (p_upper + (x[ambig_idx]-x_lower) * 
            as.numeric((p_lower-p_upper)/(x_upper-x_lower))) - p[ambig_idx]
      }
      else if (p_upper > p_lower){
        res[ambig_idx] <- (p_lower + (x[ambig_idx]-x_lower) * 
                        (p_upper-p_lower)/(x_upper-x_lower)) - p[ambig_idx]
      } else if (p_upper == p_lower){
        res[ambig_idx] <- p_upper - p[ambig_idx]
      }
    }
    return(res)
  }
  h_int_lower_zero_serialized_mat <- 
    h_int_lower_zero_serialized[1:(h_est_obj$N_param+1)]
  dim(h_int_lower_zero_serialized_mat) <- c(h_est_obj$N_param+1,1)
  h_int_upper_zero_serialized_mat <- 
    h_int_upper_zero_serialized[1:(h_est_obj$N_param+1)]
  dim(h_int_upper_zero_serialized_mat) <- c(h_est_obj$N_param+1,1)
  p_lower <- as.numeric(series_calculate(h_int_lower_zero_serialized_mat, 
                                      h_est_obj$coeff_vec, accelerate_series))
  p_upper <- 1-as.numeric(series_calculate(h_int_upper_zero_serialized_mat, 
                                      h_est_obj$coeff_vec, accelerate_series))
  if (is.na(p_lower) | is.na(p_upper)){
    return(rep(NA, length(p_vec)))
  }
  if (p_upper < p_lower){
    x_lower <- tryCatch({
      stats::uniroot(
        f = function(x) {
          series_calculate(hermite_int_lower(h_est_obj$N_param,x), 
               h_est_obj$coeff_vec, accelerate_series) - p_upper
        },
        interval = c(-100, 100)
      )$root
    },
    error = function(e) {NA})
    x_upper <- tryCatch({
      stats::uniroot(
        f = function(x) {
          1- series_calculate(hermite_int_upper(h_est_obj$N_param,x), 
             h_est_obj$coeff_vec, accelerate_series) - p_lower
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
    return(rep(NA, length(p_vec)))
  }
  max_steps <- 25
  eps_quant <- 2e-4
  x_0 <- rep(-50,length(p_vec))
  x_1 <- rep(50,length(p_vec))
  f_0 <- rep(0,length(p_vec)) - p_vec
  f_1 <- rep(1,length(p_vec)) - p_vec
  for (step in seq_len(max_steps)) {
    est  <- (x_0 + x_1)/2
    f_mid <- f_est(est,p_vec)
    mask_0 <- sign(f_mid) == sign(f_0)
    mask_1 <-  sign(f_mid) == sign(f_1)
    x_0 <- ifelse( mask_0, est, x_0)
    x_1 <- ifelse( mask_1, est, x_1)
    f_0 <- ifelse( mask_0, f_mid, f_0 )
    f_1 <- ifelse( mask_1, f_mid, f_1 )
    error_max <- max(abs(x_1 - x_0))
    if (error_max <= eps_quant) {break}
  }
  if (is.na(h_est_obj$exp_weight)) {
    est <-
    est * sqrt(h_est_obj$running_variance / (h_est_obj$num_obs - 1)) + 
      h_est_obj$running_mean
  } else {
    est <- est * sqrt(h_est_obj$running_variance) + h_est_obj$running_mean
  }
  return(est)
}

# Estimates the quantile at a vector of probability values using an 
# interpolation approximation.
#
# This helper method is intended for internal use by the 
# hermite_estimator_univar class.
quantile_helper_interpolate <- function(h_est_obj, p_vec, 
                                        accelerate_series = TRUE){
  result <- rep(NA,length(p_vec))
  p_all_vals <- cummax(shift_pvec + scale_pvec * 
    series_calculate(h_est_obj$h_int_full_domain_serialized,
                     h_est_obj$coeff_vec, accelerate_series))
  result <- stats::approx(p_all_vals,x_full_domain_serialized,
                   xout = p_vec, method="linear",ties = "ordered",
                   yleft = -10, yright = 10)$y
  if (is.na(h_est_obj$exp_weight)) {
    result <-
      result * sqrt(h_est_obj$running_variance / (h_est_obj$num_obs - 1)) +
      h_est_obj$running_mean
  } else {
    result <- result * sqrt(h_est_obj$running_variance) + h_est_obj$running_mean
  }
  return(result)
}

#' Estimates the quantiles at a vector of probability values
#' 
#' This method utilizes the estimator (13) in paper Stephanou, Michael, 
#' Varughese, Melvin and Iain Macdonald. "Sequential quantiles via Hermite 
#' series density estimation." Electronic Journal of Statistics 11.1 (2017): 
#' 570-607 <doi:10.1214/17-EJS1245>, with some modifications to improve the 
#' stability of numerical root finding. 
#'
#' @param h_est_obj A hermite_estimator_univar object.
#' @param p A numeric vector. A vector of probability values.
#' @param algorithm A string. Two possible values 'interpolate' which is faster
#' but may be less accurate or 'bisection' which is slower but potentially more
#' accurate.
#' @param accelerate_series A boolean value. If set to TRUE, the series 
#' acceleration methods described in:
#'
#' Boyd, John P., and Dennis W. Moore. "Summability methods for 
#' Hermite functions." Dynamics of atmospheres and oceans 10.1 (1986): 51-62. 
#'  
#' are applied. If set to FALSE, then standard summation is applied.
#' @return A numeric vector. The vector of quantile values associated with the
#' probabilities p.
quant.hermite_estimator_univar <- function(h_est_obj, p, 
                                   algorithm="interpolate", 
                                   accelerate_series = TRUE) {
  if (!is.numeric(p)) {
    stop("p must be numeric.")
  }
  if (length(p) < 1) {
    stop("p must contain at least one value.")
  }
  if (any(p>1) | any(p<0)) {
    stop("p must contain probabilities i.e. p>=0 and p<=1.")
  }
  if (h_est_obj$standardize_obs != TRUE) {
    stop("Quantile estimation requires standardization to be true.")
  }
  if (!(algorithm=="interpolate" | algorithm=="bisection") ) {
    stop("Algorithm must be either 'interpolate' or 'bisection'.")
  }
  result <- rep(NA, length(p))
  if (h_est_obj$num_obs < 2) {
    return(result)
  }
  if (h_est_obj$running_variance == 0) {
    return(rep(h_est_obj$running_mean, length(p)))
  }
  if (algorithm=="interpolate"){
    result <- quantile_helper_interpolate(h_est_obj,p,accelerate_series)
  }
  if (algorithm=="bisection"){
    result <- quantile_helper_bisection(h_est_obj,p,accelerate_series)
  }
  return(result)
}

#' Estimates the quantiles at a vector of probability values
#' 
#' This generic method is a convenience wrapper around the quant method
#'
#' @param x A hermite_estimator_univar object.
#' @param probs A numeric vector. A vector of probability values.
#' @param ... Optional additional arguments to the quant function namely 
#' algorithm and accelerate_series.
#' @return A numeric vector. The vector of quantile values associated with the
#' probabilities probs.
#' @export
quantile.hermite_estimator_univar <- function(x, probs = seq(0, 1, 0.25), ...){
  quant(x, p = probs, ...)
}

#' Estimates the median
#' 
#' This generic method is a convenience wrapper around the quant method to 
#' calculate the median.
#'
#' @param x A hermite_estimator_univar object.
#' @param ... Optional additional arguments to the quant function namely 
#' algorithm and accelerate_series.
#' @return A numeric value.
#' @export
median.hermite_estimator_univar <- function(x, ...){
  quant(x, p = 0.5, ...)
}

#' Estimates the Interquartile range (IQR)
#' 
#' This generic method is a convenience wrapper around the quant method to 
#' calculate the interquartile range.
#'
#' @param x A hermite_estimator_univar object.
#' @param ... Optional additional arguments to the quant function namely 
#' algorithm and accelerate_series.
#' @return A numeric value.
#' @export
IQR.hermite_estimator_univar <- function(x, ...){
  quartiles <- quant(x, p = c(0.25,0.75), ...)
  return(quartiles[2] - quartiles[1])
}

#' Prints univariate hermite_estimator object.
#' 
#'
#' @param x A hermite_estimator_univar object.
#' @param ... Unused
print.hermite_estimator_univar <- function(x, ...) {
  describe_estimator(x,"univariate")
}

summary_quantiles_helper <- function(cumulative_probs, quantiles_at_p, digits){
  cat("Estimated Quantiles:\n")
  cum_probs_size <- length(cumulative_probs)
  quantile_values <- matrix(round(quantiles_at_p,digits), nrow = 1, 
                            ncol = cum_probs_size, byrow = TRUE)
  colnames(quantile_values) <- paste0(cumulative_probs ,"%")
  rownames(quantile_values) <- ""
  print(quantile_values, quote = FALSE)
}

#' Summarizes univariate hermite_estimator object.
#' 
#' Outputs key parameters of a univariate hermite_estimator object along with
#' estimates of the mean, standard deviation and deciles of the data that
#' the object has been updated with.
#'
#' @param object A hermite_estimator_univar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Other arguments passed on to methods used in summary.
summary.hermite_estimator_univar <- function(object, 
                              digits = max(3, getOption("digits") - 3),...) {
  describe_estimator(object,"univariate")
  if (object$num_obs > 2){
    cat("\n")
    cat(paste0("Mean = ",round(object$running_mean,digits), "\n"))
    cat(paste0("Standard Deviation = ",
               round(calculate_running_std(object),digits), "\n"))
    if (object$standardize_obs == TRUE){
      cumulative_probs <- seq(10,90,10)
      quantiles_at_p <- quant(object,p=cumulative_probs/100)
      summary_quantiles_helper(cumulative_probs, quantiles_at_p, digits)
    }
  }
}

spearmans.hermite_estimator_univar <- function(h_est_obj, clipped) {
  stop("Spearman's Rho correlation estimation is not defined for the univariate 
       Hermite estimator")
}

kendall.hermite_estimator_univar <- function(h_est_obj, clipped) {
  stop("Kendall Tau correlation estimation is not defined for the univariate 
       Hermite estimator")
}
