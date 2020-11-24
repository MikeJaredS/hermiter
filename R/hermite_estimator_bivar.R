#' A class to sequentially estimate bivariate pdfs, cdfs and nonparametric 
#' correlations
#'
#' This method constructs an S3 object with methods for nonparametric estimation
#' of bivariate pdfs and cdfs along with nonparametric correlations.
#'
#' The hermite_estimator_bivar class allows the sequential or one-pass batch
#' estimation of the full bivariate probability density function and cumulative 
#' distribution function along with the Spearman's rank correlation coefficient. 
#' It is well suited to streaming data (both stationary and non-stationary) and 
#' to efficient estimation in the context of massive or distributed data sets. 
#' Indeed,estimators constructed on different subsets of a distributed data set 
#' can be consistently combined.
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
#' @return An S3 object of class hermite_estimator_bivar, with methods for
#' density function and distribution function estimation along with Spearman's
#' rank correlation estimation.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_bivar(N = 10, standardize = TRUE)
hermite_estimator_bivar <- function(N = 10, standardize=FALSE, 
                                    exp_weight_lambda=NA)
{
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
  this <- list(
    N_param = N,
    standardize_obs=standardize,
    running_mean_x = 0,
    running_mean_y = 0,
    running_variance_x = 0,
    running_variance_y = 0,
    exp_weight = exp_weight_lambda,
    num_obs=0,
    coeff_mat_bivar=matrix(rep(0,(N+1)*(N+1)),nrow = (N+1), ncol=(N+1)),
    coeff_vec_x = rep(0,N+1),
    coeff_vec_y = rep(0,N+1),
    normalization_hermite_vec=c(),
    W=c(),
    z=c()
  )
  # this$normalization_hermite_vec <-
  #   hermite_normalization(this$N_param)
  this$normalization_hermite_vec <- h_norm_serialized[1:(this$N_param+1)]
  # this$W <- W_matrix(this$N_param,this$N_param,this$normalization_hermite_vec)
  # this$z <- z_vector(this$N_param, this$normalization_hermite_vec)
  this$W <- W_serialized[1:(this$N_param+1),1:(this$N_param+1)]
  this$z <- z_serialized[1:(this$N_param+1)]
  # class(this) <- append(class(this),"hermite_estimator_bivar")
  class(this) <- "hermite_estimator_bivar"
  return(this)
}

# Matrix \eqn{W_{kl} = \int_{-\infty}^{\infty} h_k(u) \int_{-\infty}^{u} h_l(v)dvdu}
# W_matrix <- function(num_r,num_s, hermite_norm){
#   integrand <- function(x,r,s, hermite_norm){
#     hermite_vec <- hermite_function(max(r,s),x,hermite_norm)
#     hermite_integral_vec <- hermite_integral_val(s,x,hermite_vec)
#     result <- hermite_integral_vec[s+1,]*hermite_vec[r+1,]
#     return(result)
#   }
#   result <- matrix(rep(0,(num_r+1)*(num_s+1)), nrow=(num_r+1),ncol=(num_s+1),
#                    byrow = TRUE)
#   for (r in c(0:num_r)) {
#     for (s in c(0:num_s)) {
#       result[r+1,s+1] <- stats::integrate(function(t){integrand(t,r,s, 
#                                                                 hermite_norm)}, lower=-Inf,upper=Inf)$value
#     }
#   }
#   return(result)
# }

# Vector \eqn{z_{k} = \int_{-\infty}^{\infty} h_k(u)du}
# z_vector <- function(num_r, hermite_norm){
#   result <- hermite_int_full_domain(num_r)
#   return(result)
# }

#' Internal method to consistently combine the number of observations, means and 
#' variances of two bivariate Hermite estimators
#' 
#' The algorithm to combine the variances consistently comes from
#' Schubert, Erich, and Michael Gertz. "Numerically stable parallel computation 
#' of (co-) variance." Proceedings of the 30th International Conference on 
#' Scientific and Statistical Database Management. 2018.
#'
#' @param hermite_estimator1 A hermite_estimator_bivar object.
#' @param hermite_estimator2 A hermite_estimator_bivar object.
#' @return An object of class hermite_estimator_bivar
combine_moments_and_count_bivar <- function(hermite_estimator1 
                                                              ,hermite_estimator2){
  num_obs_1 <- hermite_estimator1$num_obs
  num_obs_2 <- hermite_estimator2$num_obs
  hermite_combined <- hermite_estimator_bivar(hermite_estimator1$N_param, 
                                              hermite_estimator1$standardize_obs)
  hermite_combined$num_obs <- num_obs_1 + num_obs_2
  hermite_combined$running_mean_x <- (num_obs_1*hermite_estimator1$running_mean_x +
                                        num_obs_2*hermite_estimator2$running_mean_x)/(num_obs_1+num_obs_2)
  hermite_combined$running_variance_x <- (hermite_estimator1$running_variance_x +
                                            hermite_estimator2$running_variance_x) + ((num_obs_1 * num_obs_2) / 
                                                                                        (num_obs_1 + num_obs_2)) *(hermite_estimator1$running_mean_x 
                                                                                                                   - hermite_estimator2$running_mean_x)^2
  
  hermite_combined$running_mean_y <- (num_obs_1*hermite_estimator1$running_mean_y +
                                        num_obs_2*hermite_estimator2$running_mean_y)/(num_obs_1+num_obs_2)
  hermite_combined$running_variance_y <- (hermite_estimator1$running_variance_y +
                                            hermite_estimator2$running_variance_y) + ((num_obs_1 * num_obs_2) / 
                                                                                        (num_obs_1 + num_obs_2)) *(hermite_estimator1$running_mean_y 
                                                                                                                   - hermite_estimator2$running_mean_y)^2
  return(hermite_combined)
}

#' Internal method to combine a list of standardized bivariate Hermite estimators with 
#' improved accuracy
#'
#'
#' @param hermite_estimators A list of hermite_estimator_bivar objects.
#' @return An object of class hermite_estimator_bivar.
combine_standardized_helper_bivar <- function(hermite_estimators) {
  N <- hermite_estimators[[1]]$N_param
  hermite_estimator_combined <- base::Reduce(f=combine_moments_and_count_bivar, 
                                             x = hermite_estimators)
hermite_estimator_combined$coeff_vec_x <-
  sapply(1:(N+1),FUN=function(k){sum(sapply(hermite_estimators,
                                            FUN=function(x){(x$num_obs / hermite_estimator_combined$num_obs) *
                                                gauss_hermite_quad_100(function(t){integrand_coeff_univar(t, x,
                                                                                                          hermite_estimator_combined, k, 1)})}))})

hermite_estimator_combined$coeff_vec_y <-
  sapply(1:(N+1),FUN=function(k){sum(sapply(hermite_estimators,
                                            FUN=function(x){(x$num_obs / hermite_estimator_combined$num_obs) *
                                                gauss_hermite_quad_100(function(t){integrand_coeff_univar(t, x,
                                                                                                          hermite_estimator_combined, k, 2)})}))})
herm_mod <- function(t,N){hermite_polynomial(N, t) *
    hermite_normalization(N)}
g_mat <- function(N,t_1,t_2, const_mult){sqrt(2)*const_mult*tcrossprod(hermite_function_N(N,t_1),herm_mod(t_2,N))}
g_mat_scale_shift <- function(N,scale,shift){
  res <- matrix(rep(0,(N+1)*(N+1)), nrow=N+1,ncol=N+1)
  root_x <- root_x_serialized
  weight_w <- weight_w_serialized
  for (idx in seq_along(root_x)) {
    res <- res + g_mat(N,scale*sqrt(2)*root_x[idx]+shift,sqrt(2)*root_x[idx],weight_w[idx])
  }
  return(res)
}
coeff_bivar_combine <- function(hermite_est_current,hermite_estimator_combined){
  N_in <- hermite_est_current$N_param
  prev_mean_x <- hermite_est_current$running_mean_x
  prev_mean_y <- hermite_est_current$running_mean_y
  prev_sd_x <- sqrt(hermite_est_current$running_variance_x/(hermite_est_current$num_obs-1))
  prev_sd_y <- sqrt(hermite_est_current$running_variance_y/(hermite_est_current$num_obs-1))
  new_mean_x <- hermite_estimator_combined$running_mean_x
  new_mean_y <- hermite_estimator_combined$running_mean_y
  new_sd_x <- sqrt(hermite_estimator_combined$running_variance_x/(hermite_estimator_combined$num_obs-1))
  new_sd_y <- sqrt(hermite_estimator_combined$running_variance_y/(hermite_estimator_combined$num_obs-1))
  res_1 <- g_mat_scale_shift(N_in,prev_sd_x/new_sd_x,(prev_mean_x-new_mean_x)/new_sd_x)
  res_2 <- g_mat_scale_shift(N_in,prev_sd_y/new_sd_y,(prev_mean_y-new_mean_y)/new_sd_y)
  return(res_1%*%hermite_est_current$coeff_mat_bivar%*%t(res_2))
}
coeff_mat_lst <- lapply(hermite_estimators, FUN = function(x){x$num_obs/hermite_estimator_combined$num_obs * coeff_bivar_combine(x,hermite_estimator_combined)} )
hermite_estimator_combined$coeff_mat_bivar <- Reduce(f = "+",coeff_mat_lst)
  return(hermite_estimator_combined)
}

#' Combines two bivariate Hermite estimators
#'
#' This method allows a pair of Hermite based estimators of class
#' hermite_estimator_bivar to be consistently combined.
#'
#' Note that the N and standardize arguments must be the same for the two
#' estimators in order to combine them. In addition, note that exponentially
#' weighted estimators cannot be combined. If the Hermite estimators are not
#' standardized, the combined estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator_bivar inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param this A hermite_estimator_bivar object. The first Hermite series based
#' estimator.
#' @param hermite_estimator_other A hermite_estimator_bivar object. The second 
#' Hermite series based estimator.
#' @return An object of class hermite_estimator_bivar.
#' @export
combine_pair.hermite_estimator_bivar <-
  function(this, hermite_estimator_other) {
    if (!is(hermite_estimator_other, "hermite_estimator_bivar")) {
      stop("combine_pair can only be applied to hermite_estimator_bivar objects.")
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
      hermite_estimator_combined <- combine_moments_and_count_bivar(this, 
                                                              hermite_estimator_other)
      hermite_estimator_combined$coeff_vec_x <-
        (
          this$coeff_vec_x * this$num_obs + hermite_estimator_other$coeff_vec_x 
          * hermite_estimator_other$num_obs
        ) / (this$num_obs + hermite_estimator_other$num_obs)
      hermite_estimator_combined$coeff_vec_y <-
        (
          this$coeff_vec_y * this$num_obs + hermite_estimator_other$coeff_vec_y 
          * hermite_estimator_other$num_obs
        ) / (this$num_obs + hermite_estimator_other$num_obs)
      hermite_estimator_combined$coeff_mat_bivar <-
        (
          this$coeff_mat_bivar * this$num_obs + hermite_estimator_other$coeff_mat_bivar 
          * hermite_estimator_other$num_obs
        ) / (this$num_obs + hermite_estimator_other$num_obs)
    } else {
      hermite_estimator_combined <- combine_standardized_helper_bivar(list(this,
                                                                     hermite_estimator_other))
    }
    return(hermite_estimator_combined)
  }

#' @export
combine_hermite_bivar <- function(hermite_estimators) {
  UseMethod("combine_hermite_bivar", hermite_estimators)
}

#' Combines a list of bivariate Hermite estimators
#'
#' This method allows a list of Hermite based estimators of class
#' hermite_estimator_bivar to be consistently combined.
#'
#' Note that the N and standardize arguments must be the same for all estimators
#' in order to combine them. In addition, note that exponentially weighted
#' estimators cannot be combined. If the Hermite estimators are not
#' standardized, the combined estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator_bivar inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param hermite_estimators A list of hermite_estimator_bivar objects.
#' @return An object of class hermite_estimator_bivar.
#' @export
combine_hermite_bivar.list <- function(hermite_estimators) {
  if (length(hermite_estimators) == 1) {
    return(hermite_estimators[[1]])
  }
  if (hermite_estimators[[1]]$standardize_obs==FALSE){
    hermite_estimator_combined <- base::Reduce(combine_pair,hermite_estimators)
  } else {
    hermite_estimator_combined <- combine_standardized_helper_bivar(hermite_estimators)
  }
  return(hermite_estimator_combined)
}

#' Updates the Hermite series based estimator sequentially
#'
#' This method can be applied in sequential estimation settings.
#' 
#'
#' @param this A hermite_estimator_bivar object.
#' @param x A numeric vector of length 2. A bivariate observation to be 
#' incorporated into the estimator.
#' @return An object of class hermite_estimator_bivar.
#' @export
#' @examples
#' hermite_estimator <- hermite_estimator_bivar(N = 10, standardize = TRUE)
#' hermite_estimator <- update_sequential(hermite_estimator, x = c(1,2))
update_sequential.hermite_estimator_bivar <- function(this, x)
{
  y <- x[2]
  x <- x[1]
  this$num_obs <- this$num_obs + 1
  if (this$standardize_obs == TRUE) {
    if (is.na(this$exp_weight)) {
      processed_vec <-
        standardizeInputs(x,
                          this$num_obs,
                          this$running_mean_x,
                          this$running_variance_x)
      this$running_mean_x <- processed_vec[1]
      this$running_variance_x <- processed_vec[2]
      x <- processed_vec[3]
      processed_vec <-
        standardizeInputs(y,
                          this$num_obs,
                          this$running_mean_y,
                          this$running_variance_y)
      this$running_mean_y <- processed_vec[1]
      this$running_variance_y <- processed_vec[2]
      y <- processed_vec[3]
      if (this$num_obs < 2) {
        return(this)
      }
    } else {
      processed_vec <-
        standardizeInputsEW(x,
                            this$num_obs,
                            this$exp_weight,
                            this$running_mean_x,
                            this$running_variance_x)
      this$running_mean_x <- processed_vec[1]
      this$running_variance_x <- processed_vec[2]
      x <- processed_vec[3]
      processed_vec <-
        standardizeInputsEW(y,
                            this$num_obs,
                            this$exp_weight,
                            this$running_mean_y,
                            this$running_variance_y)
      this$running_mean_y <- processed_vec[1]
      this$running_variance_y <- processed_vec[2]
      y <- processed_vec[3]
      if (this$num_obs < 2) {
        return(this)
      }
    }
  }
  h_x <-
    as.vector(hermite_function(this$N_param, x, this$normalization_hermite_vec))
  h_y <-
    as.vector(hermite_function(this$N_param, y, this$normalization_hermite_vec))
  if (is.na(this$exp_weight)) {
    this$coeff_vec_x <-
      (this$coeff_vec_x * (this$num_obs - 1) + h_x) / this$num_obs
    this$coeff_vec_y <-
      (this$coeff_vec_y * (this$num_obs - 1) + h_y) / this$num_obs
  } else {
    this$coeff_vec_x <-
      this$coeff_vec_x * (1 - this$exp_weight) + h_x * this$exp_weight
    this$coeff_vec_y <-
      this$coeff_vec_y * (1 - this$exp_weight) + h_y * this$exp_weight
  }
  if (is.na(this$exp_weight)){
    this$coeff_mat_bivar <- (this$coeff_mat_bivar*(this$num_obs-1) +
                               tcrossprod(h_x,h_y))/(this$num_obs)
  } else {
    this$coeff_mat_bivar<- this$coeff_mat_bivar*(1-this$exp_weight) +
      tcrossprod(h_x,h_y)*this$exp_weight
  }
  return(this)
}

#' Updates the Hermite series based estimator with a batch of data
#'
#' This method can be applied in one-pass batch estimation settings. This
#' method cannot be used with an exponentially weighted estimator.
#'
#' @param this A hermite_estimator_bivar object.
#' @param x A numeric matrix. A matrix of bivariate observations to be 
#' incorporated into the estimator. Each row corresponds to a single bivariate
#' observation.
#' @return An object of class hermite_estimator_bivar.
#' @export
#' @examples
#' hermite_estimator <- hermite_estimator_bivar(N = 10, standardize = TRUE)
#' hermite_estimator <- update_batch(hermite_estimator, x = matrix(c(1, 2, 3, 4,
#'  5, 6),nrow=3, ncol=2, byrow = TRUE))
update_batch.hermite_estimator_bivar <- function(this, x) {
  if (is.null(nrow(x))){
    dim(x) <- c(1,2)
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (!is.na(this$exp_weight)) {
    stop("The Hermite estimator cannot be exponentially weighted.")
  }
  this$num_obs <- nrow(x)
  if (this$standardize_obs == TRUE) {
    this$running_mean_x <-mean(x[,1])
    this$running_variance_x <-stats::var(x[,1]) * (nrow(x) - 1)
    x[,1] <-
      (x[,1] - this$running_mean_x) /
      sqrt(this$running_variance_x / (this$num_obs - 1))
    this$running_mean_y <- mean(x[,2])
    this$running_variance_y <- stats::var(x[,2]) * (nrow(x) - 1)
    x[,2] <-
      (x[,2] - this$running_mean_y) / sqrt(this$running_variance_y /
                                             (this$num_obs - 1))
  }
  h_x <-
    hermite_function(this$N_param, x[,1], this$normalization_hermite_vec)
  h_y <-
    hermite_function(this$N_param, x[,2], this$normalization_hermite_vec)
  this$coeff_vec_x <- rowSums(h_x) / this$num_obs
  this$coeff_vec_y <- rowSums(h_y) / this$num_obs
  this$coeff_mat_bivar <- tcrossprod(h_x,h_y) / this$num_obs
  return(this)
}

# Internal helper method to calculate the probability density at a single 2-d
# x value
dens_helper <- function(this, x, clipped)
{
  UseMethod("dens_helper",this)
}

dens_helper.hermite_estimator_bivar <- function(this,x, clipped = FALSE){
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (this$num_obs < 2) {
    return(NA)
  }
  y <- x[2]
  x <- x[1]
  factor <- 1
  if (this$standardize_obs == TRUE) {
    if (is.na(this$exp_weight)) {
      running_std_x <- sqrt(this$running_variance_x / (this$num_obs - 1))
      running_std_y <- sqrt(this$running_variance_y / (this$num_obs - 1))
    } else {
      running_std_x <- sqrt(this$running_variance_x)
      running_std_y <- sqrt(this$running_variance_y)
    }
    x <- (x - this$running_mean_x) / running_std_x
    y <- (y - this$running_mean_y) / running_std_y
    factor <- 1 / (running_std_x * running_std_y)
  }
  return(factor * t(hermite_function(this$N_param, x, this$normalization_hermite_vec)) %*%
           this$coeff_mat_bivar  %*% hermite_function(this$N_param, y,
                                                      this$normalization_hermite_vec))
}

#' Estimates the probability densities for a matrix of 2-d x values
#'
#' This method calculates the probability density values for a matrix of 
#' 2-d x vector values using the hermite_estimator_bivar object (this).
#'
#' The object must be updated with observations prior to the use of the method.
#'
#' @param this A hermite_estimator_bivar object.
#' @param x A numeric matrix. Each row corresponds to a 2-d coordinate.
#' @param clipped A boolean value. This value determines whether
#' probability densities are clipped to be bigger than zero.
#' @return A numeric vector of probability density values.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_bivar(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, matrix(rnorm(30*2), nrow=30, ncol=2, 
#' byrow = TRUE))
#' cdf_est <- dens(hermite_est, matrix(c(0, 0, 1, 1, 2, 2), nrow=3, ncol=2, 
#' byrow = TRUE))
dens.hermite_estimator_bivar <- function(this,x, clipped = FALSE){
  if (is.null(nrow(x))){
    dim(x) <- c(1,2)
  }
  result <- rep(0,nrow(x))
  for (idx in 1:nrow(x)) {
    result[idx] <- dens_helper(this,x[idx,])
  }
  if (clipped == TRUE) {
    result <- pmax(result, 1e-08)
  }
  return(result)
}

# Internal helper method to calculate the cumulative probability at a single 2-d
# x value
cum_prob_helper <- function(this, x, clipped = FALSE)
{
  UseMethod("cum_prob_helper",this)
}

cum_prob_helper.hermite_estimator_bivar <- function(this,x, clipped = FALSE){
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (this$num_obs < 2) {
    return(NA)
  }
  y <- x[2]
  x <- x[1]
  if (this$standardize_obs == TRUE) {
    if (is.na(this$exp_weight)) {
      running_std_x <- sqrt(this$running_variance_x / (this$num_obs - 1))
      running_std_y <- sqrt(this$running_variance_y / (this$num_obs - 1))
    } else {
      running_std_x <- sqrt(this$running_variance_x)
      running_std_y <- sqrt(this$running_variance_y)
    }
    x <- (x - this$running_mean_x) / running_std_x
    y <- (y - this$running_mean_y) / running_std_y
  }
  return(t(hermite_int_lower(N = this$N_param,x = x)) %*%
           this$coeff_mat_bivar  %*% hermite_int_lower(N = this$N_param,x = y))
}

#' Estimates the cumulative probabilities for a matrix of 2-d x values
#'
#' This method calculates the cumulative probability values for a matrix of 
#' 2-d x vector values using the hermite_estimator_bivar object (this).
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param this A hermite_estimator_bivar object.
#' @param x A numeric matrix. Each row corresponds to a 2-d coordinate.
#' @param clipped A boolean value. This value determines whether cumulative
#' probabilities are clipped to lie within the range [0,1].
#' @return A numeric vector of cumulative probability values.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_bivar(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, matrix(rnorm(30*2), nrow=30, ncol=2, 
#' byrow = TRUE))
#' cdf_est <- cum_prob(hermite_est, matrix(c(0, 0, 1, 1, 2, 2), nrow=3, ncol=2, 
#' byrow = TRUE))
cum_prob.hermite_estimator_bivar <- function(this,x, clipped = FALSE){
  if (is.null(nrow(x))){
    dim(x) <- c(1,2)
  }
  result <- rep(0,nrow(x))
  for (idx in 1:nrow(x)) {
    result[idx] <- cum_prob_helper(this,x[idx,])
  }
  if (clipped == TRUE) {
    result <- pmin(pmax(result, 1e-08), 1)
  }
  return(result)
}

#' Estimates the Spearman's rank correlation coefficient
#'
#' This method calculates the Spearman's rank correlation coefficient value
#' using the hermite_estimator_bivar object (this).
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param this A hermite_estimator_bivar object.
#' @return A numeric value.
#' @export
#' @examples
#' hermite_est <- hermite_estimator_bivar(N = 10, standardize = TRUE)
#' hermite_est <- update_batch(hermite_est, matrix(rnorm(30*2), nrow=30, ncol=2, 
#' byrow = TRUE))
#' spearmans_est <- spearmans(hermite_est)
spearmans <- function(this, clipped = FALSE)
{
  UseMethod("spearmans",this)
}

#' @export
spearmans.hermite_estimator_bivar <- function(this, clipped = FALSE)
{
  if (this$num_obs < 2) {
    return(NA)
  }
  W_transpose <- t(this$W)
  result <- 12*t(this$coeff_vec_x) %*% W_transpose %*%this$coeff_mat_bivar %*% 
    this$W%*%this$coeff_vec_y +
    -6 * t(this$coeff_vec_x)%*%W_transpose%*%this$coeff_mat_bivar%*%this$z +
    -6 *t(this$z)%*% this$coeff_mat_bivar%*%this$W%*%this$coeff_vec_y +
    3 * t(this$z)%*% this$coeff_mat_bivar%*%this$z
  if (clipped == TRUE) {
    result <- pmin(pmax(result, -1), 1)
  }
  return(as.numeric(result))
}
