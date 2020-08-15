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
      stop("N must be greater than or equal to 0 and smaller than or equal to 100.")
    }
    if (!(standardize == TRUE | standardize == FALSE)) {
      stop("standardize can only take on values TRUE or FALSE.")
    }
    if (!is.na(exp_weight_lambda)) {
      if (!is.numeric(exp_weight_lambda)) {
        stop("exp_weight_lambda must be numeric.")
      }
      if (exp_weight_lambda <= 0 | exp_weight_lambda > 1) {
        stop("exp_weight_lambda must be a real number greater than 0 and less than or equal to 1.")
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
#' equivalence will be approximate.
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
    hermite_estimator_combined <-
      hermite_estimator(N = this$N_param,
                        standardize = this$standardize_obs)
    hermite_estimator_combined$coeff_vec <-
      (
        this$coeff_vec * this$num_obs + hermite_estimator_other$coeff_vec * hermite_estimator_other$num_obs
      ) / (this$num_obs + hermite_estimator_other$num_obs)
    hermite_estimator_combined$num_obs <-
      this$num_obs + hermite_estimator_other$num_obs
    hermite_estimator_combined$running_mean <-
      (
        this$running_mean * this$num_obs + hermite_estimator_other$running_mean * hermite_estimator_other$num_obs
      ) / (this$num_obs + hermite_estimator_other$num_obs)
    hermite_estimator_combined$running_variance <-
      this$running_variance + hermite_estimator_other$running_variance
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
#' equivalence will be approximate.
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
  hermite_estimator_combined <- hermite_estimators[[1]]
  for (idx in c(2:length(hermite_estimators))) {
    hermite_estimator_combined <-
      combine_pair(hermite_estimator_combined, hermite_estimators[[idx]])
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
#' hermite_estimator <- update_sequential(hermite_estimator(N = 10, standardize = TRUE), x = 2)
update_sequential <- function(this, x) {
  UseMethod("update_sequential", this)
}

#' @export
update_sequential.hermite_estimator <- function(this, x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x) != 1) {
    stop("The sequential method is only applicable to one observation at a time.")
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
#' hermite_estimator <- update_batch(hermite_estimator(N = 10, standardize = TRUE), x = c(1, 2))
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
      (this$running_mean * (this$num_obs - length(x)) + length(x) * mean(x)) / this$num_obs
    this$running_variance <-
      (
        this$running_variance * (this$num_obs - length(x)) + length(x) * stats::var(x) * (length(x) - 1)
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

#' Estimates the cumulative probability for quantile estimation
#'
#' This helper method uses a modified distribution function estimator which
#' differs from the cum_prob.hermite_estimator.
#'
#' The modified distribution function estimator appears more accurate for
#' quantile estimation as validated empirically in:
#'
#' \url{https://projecteuclid.org/euclid.ejs/1488531636}
#'
#' This method is intended for internal use by the hermite_estimator class.
#'
#' @param this A hermite_estimator object.
#' @param x A numeric vector.
#' @return A numeric vector of cumulative probability values.
#' @export
cum_prob_quantile_helper <- function(this, x) {
  UseMethod("cum_prob_quantile_helper", this)
}

#' @export
cum_prob_quantile_helper.hermite_estimator <- function(this, x) {
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  h_k <-
    hermite_function(this$N_param, x, this$normalization_hermite_vec)
  integrals_hermite <-
    hermite_integral_val_quantile_adap(this$N_param, x, h_k)
  cdf_val <-
    1 - as.vector(as.vector(this$coeff_vec) %*% integrals_hermite)
  return(cdf_val)
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
  est <- tryCatch({
    stats::uniroot(
      f = function(x) {
        cum_prob_quantile_helper(this, x) - p
      },
      interval = c(-100, 100)
    )$root
  },
  warning = function(w) {
    
  },
  error = function(e) {
    NA
  },
  finally = {
    
  })
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
#' @param this A hermite_estimator object.
#' @param p A numeric vector. A vector of probability values.
#' @return A numeric vector. The vector of quantile values associated with the
#' probabilites p.
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
  for (idx in c(1:length(p))) {
    result[idx] <- quantile_helper(this, p[idx])
  }
  return(result)
}
