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
#' @param N An integer between 0 and 75. The Hermite series based estimator
#' is truncated at N+1 terms.
#' @param standardize A boolean value. Determines whether the observations are
#' standardized, a transformation which often improves performance.
#' @param exp_weight_lambda A numerical value between 0 and 1. This parameter
#' controls the exponential weighting of the Hermite series based estimator.
#' If this parameter is NA, no exponential weighting is applied.
#' @return An S3 object of class hermite_estimator_univar, with methods for
#' density function, distribution function and quantile function estimation.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_univar(N = 10, standardize = TRUE)
hermite_estimator_univar <-
  function(N = 10,
           standardize = FALSE,
           exp_weight_lambda = NA) {
    if (!is.numeric(N)) {
      stop("N must be numeric.")
    }
    if (N < 0 | N > 75) {
      stop("N must be >= 0 and N <= 75.")
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
        normalization_hermite_vec = c(),
        h_int_lower_serialized = c(),
        h_int_upper_serialized = c()
      )
    this$normalization_hermite_vec <- h_norm_serialized[1:(this$N_param+1)]
    this$h_int_lower_serialized <- h_int_lower_serialized[1:(this$N_param+1),]
    this$h_int_upper_serialized <- h_int_upper_serialized[1:(this$N_param+1),]
    class(this) <- c("hermite_estimator_univar", "list")
    return(this)
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
#' @param this A hermite_estimator_univar object. The first Hermite series based
#' estimator.
#' @param hermite_estimator_other A hermite_estimator_univar object. The 
#' second Hermite series based estimator.
#' @return An object of class hermite_estimator_univar.
#' @export
#' @examples
#' hermite_est_1 <- hermite_estimator_univar(N = 10, standardize = FALSE)
#' hermite_est_1 <- update_batch(hermite_est_1, rnorm(30))
#' hermite_est_2 <- hermite_estimator_univar(N = 10, standardize = FALSE)
#' hermite_est_2 <- update_batch(hermite_est_2, rnorm(30))
#' hermite_merged <- merge_pair(hermite_est_1, hermite_est_2)
merge_pair.hermite_estimator_univar <-
  function(this, hermite_estimator_other) {
    if (!is(hermite_estimator_other, "hermite_estimator_univar")) {
      stop("merge_pair.hermite_estimator_univar can only be applied to 
           hermite_estimator_univar objects.")
    }
    if (this$N_param != hermite_estimator_other$N_param) {
      stop("N must be equal to merge estimators.")
    }
    if (this$standardize_obs != hermite_estimator_other$standardize_obs) {
      stop("Standardization setting must be the same to merge estimators.")
    }
    if (!is.na(this$exp_weight) |
        !is.na(hermite_estimator_other$exp_weight)) {
      stop("Cannot merge exponentially weighted estimators.")
    }
    if (this$standardize_obs == FALSE) {
      hermite_estimator_merged <- merge_moments_and_count_univar(this, 
                                                      hermite_estimator_other)
      hermite_estimator_merged$coeff_vec <-
        (
          this$coeff_vec * this$num_obs + hermite_estimator_other$coeff_vec 
          * hermite_estimator_other$num_obs
        ) / (this$num_obs + hermite_estimator_other$num_obs)
    } else {
      hermite_estimator_merged <- 
        merge_standardized_helper_univar(list(this,hermite_estimator_other))
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
#' @export
#' @examples
#' hermite_est_1 <- hermite_estimator_univar(N = 10, standardize = FALSE)
#' hermite_est_1 <- update_batch(hermite_est_1, rnorm(30))
#' hermite_est_2 <- hermite_estimator_univar(N = 10, standardize = FALSE)
#' hermite_est_2 <- update_batch(hermite_est_2, rnorm(30))
#' hermite_merged <- merge_hermite(list(hermite_est_1, hermite_est_2))
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

#' Updates the Hermite series based estimator sequentially
#'
#' This method can be applied in sequential estimation settings.
#' 
#'
#' @param this A hermite_estimator_univar object.
#' @param x A numeric value. An observation to be incorporated into the
#' estimator.
#' @return An object of class hermite_estimator_univar.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_univar(N = 10, standardize = TRUE)
#' hermite_est <- update_sequential(hermite_est, x = 2)
update_sequential.hermite_estimator_univar <- function(this, x) {
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
      prev_running_mean <- this$running_mean
      this$running_mean <-  (this$running_mean * (this$num_obs - 1) + x) /
        this$num_obs
      if (this$num_obs < 2) {
        return(this)
      }
      this$running_variance <- this$running_variance + (x - prev_running_mean) *
        (x - this$running_mean)
      x <- (x - this$running_mean) /  sqrt(this$running_variance /
                                            (this$num_obs - 1))
    } else {
      if (this$num_obs < 2){
        this$running_mean <- x
        this$running_variance <- 1
        return(this)
      }
      this$running_mean <-  (1 - this$exp_weight) * this$running_mean +
        this$exp_weight * x
      this$running_variance <- (1 - this$exp_weight) * this$running_variance +
        this$exp_weight * (x - this$running_mean)^2
      x <- (x - this$running_mean) / sqrt(this$running_variance)
    }
  }
  h_k <-
    as.vector(hermite_function_N(this$N_param, x, 
                                 this$normalization_hermite_vec))
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
#' @param this A hermite_estimator_univar object.
#' @param x A numeric vector. A vector of observations to be incorporated
#' into the estimator.
#' @return An object of class hermite_estimator_univar.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_univar(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, x = c(1, 2))
update_batch.hermite_estimator_univar <- function(this, x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (!is.na(this$exp_weight)) {
    stop("The Hermite estimator cannot be exponentially weighted.")
  }
  if (length(x) < 1) {
    stop("x must contain at least one value.")
  }
  this$num_obs <- length(x)
  if (this$standardize_obs == TRUE) {
    this$running_mean <- mean(x)
    this$running_variance <- stats::var(x) * (length(x) - 1)
    x <-
      (x - this$running_mean) / sqrt(this$running_variance / (this$num_obs - 1))
  }
  this$coeff_vec <- hermite_function_sum_N(this$N_param, x,
                                 this$normalization_hermite_vec) / this$num_obs
  return(this)
}

calculate_running_std.hermite_estimator_univar<- function(this){
  if (is.na(this$exp_weight)) {
    running_std <- sqrt(this$running_variance / (this$num_obs - 1))
  } else {
    running_std <- sqrt(this$running_variance)
  }
  return(running_std)
}

#' Estimates the probability density for a vector of x values
#'
#' This method calculates the probability density values at a vector of
#' x values using the hermite_estimator_univar object (this).
#'
#' The object must be updated with observations prior to the use of the method.
#'
#' @param this A hermite_estimator_univar object.
#' @param x A numeric vector. Values at which to estimate the probability
#' density.
#' @param clipped A boolean value. This value determines whether
#' probability densities are clipped to be bigger than zero.
#' @return A numeric vector of probability density values.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_univar(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, rnorm(30))
#' pdf_est <- dens(hermite_est, c(0, 0.5, 1))
dens.hermite_estimator_univar <- function(this, x, clipped = FALSE) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x) < 1) {
    stop("x must contain at least one value.")
  }
  if (this$num_obs < 2) {
    return(rep(NA, length(x)))
  }
  factor <- 1
  if (this$standardize_obs == TRUE) {
    running_std <- calculate_running_std(this)
    x <- (x - this$running_mean) / running_std
    factor <- 1 / running_std
  }
  h_k <-
    hermite_function_N(this$N_param, x, this$normalization_hermite_vec)
  pdf_val <- crossprod(h_k, this$coeff_vec) * factor
  if (clipped == TRUE) {
    pdf_val <- pmax(pdf_val, 1e-08)
  }
  return(as.vector(pdf_val))
}

#' Estimates the cumulative probability for a vector of x values
#'
#' This method calculates the cumulative probability values at a vector of
#' x values using the hermite_estimator_univar object (this).
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param this A hermite_estimator_univar object.
#' @param x A numeric vector. Values at which to estimate the cumulative
#' probability
#' @param clipped A boolean value. This value determines whether cumulative
#' probabilities are clipped to lie within the range [0,1].
#' @return A numeric vector of cumulative probability values.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_univar(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, rnorm(30))
#' cdf_est <- cum_prob(hermite_est, c(0, 0.5, 1))
cum_prob.hermite_estimator_univar <- function(this, x, clipped = FALSE) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x) < 1) {
    stop("x must contain at least one value.")
  }
  if (this$num_obs < 2) {
    return(rep(NA, length(x)))
  }
  if (this$standardize_obs == TRUE) {
    running_std <- calculate_running_std(this)
    x <- (x - this$running_mean) / running_std
  }
  h_k <-
    hermite_function_N(this$N_param, x, this$normalization_hermite_vec)
  integrals_hermite <- hermite_int_lower(this$N_param, x, 
                                         hermite_function_matrix = h_k)
  cdf_val <- crossprod(integrals_hermite, this$coeff_vec)
  if (clipped == TRUE) {
    cdf_val <- pmin(pmax(cdf_val, 1e-08), 1)
  }
  return(as.vector(cdf_val))
}

# Estimates the quantile at a single probability value
#
# This helper method is intended for internal use by the 
# hermite_estimator_univar class.
quantile_helper <- function(this, p_vec, p_lower, p_upper, x_lower,
                            x_upper) {
  f_est <- function(x,p) {
    lower_idx <- which(x < x_lower)
    upper_idx <- which(x > x_upper)
    ambig_idx <- which(x >= x_lower & x <= x_upper)
    res <- rep(NA,length(x))
    if (length(lower_idx)>0){
      res[lower_idx] <- crossprod(hermite_int_lower(this$N_param, 
           x[lower_idx], normalization_hermite = 
             this$normalization_hermite_vec), this$coeff_vec) - p[lower_idx]
    }
    if (length(upper_idx)>0){
      res[upper_idx] <- 1 - 
        crossprod(hermite_int_upper(this$N_param, 
           x[upper_idx], normalization_hermite = 
             this$normalization_hermite_vec), this$coeff_vec) - p[upper_idx]
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
  # Vectorized bisection search:
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
#' @param this A hermite_estimator_univar object.
#' @param p A numeric vector. A vector of probability values.
#' @return A numeric vector. The vector of quantile values associated with the
#' probabilities p.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_univar(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, rnorm(30))
#' quant_est <- quant(hermite_est, c(0.25, 0.5, 0.75))
quant.hermite_estimator_univar <- function(this, p) {
  if (!is.numeric(p)) {
    stop("p must be numeric.")
  }
  if (length(p) < 1) {
    stop("p must contain at least one value.")
  }
  if (any(p>1) | any(p<0)) {
    stop("p must contain probabilities i.e. p>=0 and p<=1.")
  }
  if (this$standardize_obs != TRUE) {
    stop("Quantile estimation requires standardization to be true.")
  }
  if (this$num_obs < 2) {
    return(rep(NA, length(p)))
  }
  # p_lower <- as.numeric(crossprod(this$coeff_vec, 
  #                             h_int_lower_zero_serialized[1:(this$N_param+1)]))
  # p_upper <- 1-as.numeric(crossprod(this$coeff_vec, 
  #                             h_int_upper_zero_serialized[1:(this$N_param+1)]))
  # if (is.na(p_lower) | is.na(p_upper)){
  #   return(rep(NA, length(p)))
  # }
  # if (p_upper < p_lower){
  #   x_lower <- tryCatch({
  #     stats::uniroot(
  #       f = function(x) {
  #         crossprod(this$coeff_vec,hermite_int_lower(this$N_param,x,
  #         normalization_hermite = this$normalization_hermite_vec)) - p_upper
  #       },
  #       interval = c(-100, 100)
  #     )$root
  #   },
  #   error = function(e) {NA})
  #   x_upper <- tryCatch({
  #     stats::uniroot(
  #       f = function(x) {
  #         1-as.numeric(crossprod(this$coeff_vec, 
  #                        hermite_int_upper(this$N_param,x,
  #         normalization_hermite = this$normalization_hermite_vec))) - p_lower
  #       },
  #       interval = c(-100, 100)
  #     )$root
  #   },
  #   error = function(e) {NA})
  # } else if (p_upper > p_lower) {
  #   x_lower <- -1e-6
  #   x_upper <- 1e-6
  # } else if (p_upper == p_lower){
  #   x_lower <- 0
  #   x_upper <- 0
  # }
  # if (is.na(x_lower) | is.na(x_upper)){
  #   return(rep(NA, length(p)))
  # }
  # result <- quantile_helper(this, p, p_lower, p_upper, x_lower,
  #                           x_upper)
  result <- rep(NA,length(p))
  coeffs <- as.numeric(this$coeff_vec)
  p_lower_vals <- crossprod(this$h_int_lower_serialized, coeffs)
  p_upper_vals <- 1-crossprod(this$h_int_upper_serialized, coeffs)
  # p_lower_vals <- arma_mm(h_int_lower_t,coeffs)
  # p_upper_vals <- 1-arma_mm(h_int_upper_t,coeffs)
  p_all_vals <- cummax(c(p_lower_vals,p_upper_vals))
  res <- findInterval(p,p_all_vals)
  result <- x_full_domain_serialized[res]
  if (is.na(this$exp_weight)) {
    result <-
      result * sqrt(this$running_variance / (this$num_obs - 1)) +
      this$running_mean
  } else {
    result <- result * sqrt(this$running_variance) + this$running_mean
  }
  return(result)
}

spearmans.hermite_estimator_univar <- function(this, clipped) {
  stop("Spearman's Rho correlation estimation is not defined for the univariate 
       Hermite estimator")
}

kendall.hermite_estimator_univar <- function(this, clipped) {
  stop("Kendall Tau correlation estimation is not defined for the univariate 
       Hermite estimator")
}
