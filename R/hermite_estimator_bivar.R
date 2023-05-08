#' A class to sequentially estimate bivariate pdfs, cdfs and nonparametric 
#' correlations
#'
#' This method constructs an S3 object with methods for nonparametric estimation
#' of bivariate pdfs and cdfs along with nonparametric correlations.
#'
#' The hermite_estimator_bivar class allows the sequential or one-pass batch
#' estimation of the full bivariate probability density function and cumulative 
#' distribution function along with the Spearman's rank correlation 
#' coefficient. It is well suited to streaming data (both stationary and 
#' non-stationary) and to efficient estimation in the context of massive or 
#' distributed data sets. Indeed,estimators constructed on different subsets 
#' of a distributed data set can be consistently merged.
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
#' @param observations A numeric matrix. A matrix of bivariate observations to 
#' be incorporated into the estimator. Each row corresponds to a single 
#' bivariate observation.
#' @return An S3 object of class hermite_estimator_bivar, with methods for
#' density function and distribution function estimation along with Spearman's
#' rank correlation estimation.
hermite_estimator_bivar <- function(N = 30, standardize=TRUE, 
                                    exp_weight_lambda=NA, observations = c())
{
  h_est_obj <- list(
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
  h_est_obj$normalization_hermite_vec <- 
    h_norm_serialized[1:(h_est_obj$N_param+1)]
  class(h_est_obj) <- c("hermite_estimator_bivar", "list")
  if (length(observations) > 0) {
    if (!is.numeric(observations)) {
      stop("x must be numeric.")
    }
    if (is.null(nrow(observations))){
      if (length(observations) != 2){
        stop("x must be comprised of bivariate observations.")
      }
      dim(observations) <- c(1,2)
    }
    if (ncol(observations)!=2 | nrow(observations) < 1){
      stop("x must be comprised at least one bivariate observation.")
    }
    if (any(is.na(observations)) | any(!is.finite(observations))){
      stop("The batch initialization method is only 
         applicable to finite, non NaN, non NA values.")
    }
    h_est_obj <- initialize_batch_bivar(h_est_obj, observations)
  }
  return(h_est_obj)
}

#' Internal method to consistently merge the number of observations, means and 
#' variances of two bivariate Hermite estimators
#' 
#' The algorithm to merge the variances consistently comes from
#' Schubert, Erich, and Michael Gertz. "Numerically stable parallel computation 
#' of (co-) variance." Proceedings of the 30th International Conference on 
#' Scientific and Statistical Database Management. 2018.
#'
#' @param hermite_estimator1 A hermite_estimator_bivar object.
#' @param hermite_estimator2 A hermite_estimator_bivar object.
#' @return An object of class hermite_estimator_bivar
merge_moments_and_count_bivar <- function(hermite_estimator1 
                                            ,hermite_estimator2){
  num_obs_1 <- hermite_estimator1$num_obs
  num_obs_2 <- hermite_estimator2$num_obs
  hermite_merged <- hermite_estimator_bivar(hermite_estimator1$N_param, 
                                            hermite_estimator1$standardize_obs)
  hermite_merged$num_obs <- num_obs_1 + num_obs_2
  hermite_merged$running_mean_x <- (num_obs_1 * 
                                        hermite_estimator1$running_mean_x +
                                        num_obs_2 * 
                                      hermite_estimator2$running_mean_x) / 
    (num_obs_1+num_obs_2)
  hermite_merged$running_variance_x <- (hermite_estimator1$running_variance_x + 
                                        hermite_estimator2$running_variance_x)+
    ((num_obs_1 * num_obs_2) /(num_obs_1 + num_obs_2)) * 
    (hermite_estimator1$running_mean_x - hermite_estimator2$running_mean_x)^2
  
  hermite_merged$running_mean_y <- (num_obs_1 * 
                                        hermite_estimator1$running_mean_y + 
                                      num_obs_2 * 
                                        hermite_estimator2$running_mean_y) / 
    (num_obs_1+num_obs_2)
  hermite_merged$running_variance_y <-(hermite_estimator1$running_variance_y +
                                      hermite_estimator2$running_variance_y) + 
    ((num_obs_1 * num_obs_2) /(num_obs_1 + num_obs_2)) * 
    (hermite_estimator1$running_mean_y - hermite_estimator2$running_mean_y)^2
  return(hermite_merged)
}

#' Internal method to merge a list of standardized bivariate Hermite 
#' estimators
#'
#'
#' @param hermite_estimators A list of hermite_estimator_bivar objects.
#' @return An object of class hermite_estimator_bivar.
merge_standardized_helper_bivar <- function(hermite_estimators) {
  all_N <- lapply(hermite_estimators, FUN =
                    function(x){return(x$N_param)})
  if (length(unique(all_N)) >1) {
    stop("List must contain Hermite estimators with a consistent N")
  }
  N <- hermite_estimators[[1]]$N_param
  hermite_estimator_merged <- base::Reduce(f=merge_moments_and_count_bivar, 
                                             x = hermite_estimators)
  hermite_estimator_merged$coeff_vec_x <-
    vapply(1:(N+1),FUN=function(k){sum(vapply(hermite_estimators,
       FUN=function(x){(x$num_obs / hermite_estimator_merged$num_obs) *
      gauss_hermite_quad_100(function(t){integrand_coeff_univar(t, x, 
                                       hermite_estimator_merged, k, 1)})}, 
      FUN.VALUE=numeric(1)))}, FUN.VALUE=numeric(1))
  hermite_estimator_merged$coeff_vec_y <-
    vapply(1:(N+1),FUN=function(k){sum(vapply(hermite_estimators,
       FUN=function(x){(x$num_obs / hermite_estimator_merged$num_obs) *
      gauss_hermite_quad_100(function(t){integrand_coeff_univar(t, x, 
                                      hermite_estimator_merged, k, 2)})}, 
      FUN.VALUE=numeric(1)))}, FUN.VALUE=numeric(1))
  herm_mod <- function(t,N){hermite_polynomial(N, t) *
      hermite_normalization(N)}
  g_mat <- function(N,t_1,t_2, const_mult){sqrt(2) * const_mult * 
      tcrossprod(hermite_function(N,t_1),herm_mod(t_2,N))}
  g_mat_scale_shift <- function(N,scale,shift){
    res <- matrix(rep(0,(N+1)*(N+1)), nrow=N+1,ncol=N+1)
    root_x <- root_x_serialized
    weight_w <- weight_w_serialized
    for (idx in seq_along(root_x)) {
      res <- res + g_mat(N,scale*sqrt(2) * root_x[idx] + shift,sqrt(2) * 
                           root_x[idx],weight_w[idx])
    }
    return(res)
  }
  coeff_bivar_merge <- function(hermite_est_current,hermite_estimator_merged){
    N_in <- hermite_est_current$N_param
    prev_mean_x <- hermite_est_current$running_mean_x
    prev_mean_y <- hermite_est_current$running_mean_y
    prev_sd_x <- sqrt(hermite_est_current$running_variance_x / 
                        (hermite_est_current$num_obs-1))
    prev_sd_y <- sqrt(hermite_est_current$running_variance_y / 
                        (hermite_est_current$num_obs-1))
    new_mean_x <- hermite_estimator_merged$running_mean_x
    new_mean_y <- hermite_estimator_merged$running_mean_y
    new_sd_x <- sqrt(hermite_estimator_merged$running_variance_x / 
                       (hermite_estimator_merged$num_obs-1))
    new_sd_y <- sqrt(hermite_estimator_merged$running_variance_y / 
                       (hermite_estimator_merged$num_obs-1))
    res_1 <- g_mat_scale_shift(N_in,prev_sd_x/new_sd_x,(prev_mean_x-new_mean_x)/
                                 new_sd_x)
    res_2 <- g_mat_scale_shift(N_in,prev_sd_y/new_sd_y,(prev_mean_y-new_mean_y)/
                                 new_sd_y)
    return(res_1%*%(hermite_est_current$coeff_mat_bivar%*%t(res_2)))
  }
  coeff_mat_lst <- lapply(hermite_estimators, FUN = function(x){x$num_obs / 
      hermite_estimator_merged$num_obs * 
      coeff_bivar_merge(x,hermite_estimator_merged)} )
  hermite_estimator_merged$coeff_mat_bivar <- Reduce(f = "+",coeff_mat_lst)
  return(hermite_estimator_merged)
}

#' Merges two bivariate Hermite estimators
#'
#' This method allows a pair of Hermite based estimators of class
#' hermite_estimator_bivar to be consistently merged.
#'
#' Note that the N and standardize arguments must be the same for the two
#' estimators in order to merge them. In addition, note that exponentially
#' weighted estimators cannot be merged. If the Hermite estimators are not
#' standardized, the merged estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator_bivar inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param h_est_obj A hermite_estimator_bivar object. The first Hermite series 
#' based estimator.
#' @param hermite_estimator_other A hermite_estimator_bivar object. The second 
#' Hermite series based estimator.
#' @return An object of class hermite_estimator_bivar.
merge_pair.hermite_estimator_bivar <-
  function(h_est_obj, hermite_estimator_other) {
    if (!is(hermite_estimator_other, "hermite_estimator_bivar")) {
      stop("merge_pair can only be applied to hermite_estimator_bivar 
           objects.")
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
      hermite_estimator_merged <- merge_moments_and_count_bivar(h_est_obj, 
                                                       hermite_estimator_other)
      hermite_estimator_merged$coeff_vec_x <-
        (
          h_est_obj$coeff_vec_x * h_est_obj$num_obs + 
           hermite_estimator_other$coeff_vec_x * hermite_estimator_other$num_obs
        ) / (h_est_obj$num_obs + hermite_estimator_other$num_obs)
      hermite_estimator_merged$coeff_vec_y <-
        (
          h_est_obj$coeff_vec_y * h_est_obj$num_obs + 
           hermite_estimator_other$coeff_vec_y * hermite_estimator_other$num_obs
        ) / (h_est_obj$num_obs + hermite_estimator_other$num_obs)
      hermite_estimator_merged$coeff_mat_bivar <-
        (
          h_est_obj$coeff_mat_bivar * h_est_obj$num_obs + 
      hermite_estimator_other$coeff_mat_bivar * hermite_estimator_other$num_obs
        ) / (h_est_obj$num_obs + hermite_estimator_other$num_obs)
    } else {
     hermite_estimator_merged <- merge_standardized_helper_bivar(list(h_est_obj,
                                                    hermite_estimator_other))
    }
    return(hermite_estimator_merged)
  }

#' Merges a list of bivariate Hermite estimators
#'
#' This method allows a list of Hermite based estimators of class
#' hermite_estimator_bivar to be consistently merged.
#'
#' Note that the N and standardize arguments must be the same for all estimators
#' in order to merge them. In addition, note that exponentially weighted
#' estimators cannot be merged. If the Hermite estimators are not
#' standardized, the merged estimator will be exactly equivalent to
#' constructing a single estimator on the data set formed by combining the
#' data sets used to update the respective hermite_estimator_bivar inputs.
#' If the input Hermite estimators are standardized however, then the
#' equivalence will be approximate but still accurate in most cases.
#'
#' @param hermite_estimators A list of hermite_estimator_bivar objects.
#' @return An object of class hermite_estimator_bivar.
merge_hermite_bivar <- function(hermite_estimators) {
  if (length(hermite_estimators) == 1) {
    return(hermite_estimators[[1]])
  }
  if (hermite_estimators[[1]]$standardize_obs==FALSE){
    hermite_estimator_merged <- base::Reduce(merge_pair,hermite_estimators)
  } else {
    hermite_estimator_merged <- 
      merge_standardized_helper_bivar(hermite_estimators)
  }
  return(hermite_estimator_merged)
}

update_sequential_bivar_helper <- function(h_est_obj,x){
  h_est_obj$num_obs <- h_est_obj$num_obs + 1
  if (h_est_obj$standardize_obs == TRUE) {
    if (is.na(h_est_obj$exp_weight)) {
      prev_mean <- c(h_est_obj$running_mean_x,h_est_obj$running_mean_y)
      upd_mean <- (prev_mean*(h_est_obj$num_obs-1) + x)/h_est_obj$num_obs
      h_est_obj$running_mean_x <- upd_mean[1]
      h_est_obj$running_mean_y <- upd_mean[2]
      if (h_est_obj$num_obs < 2) {
        return(h_est_obj)
      }
      upd_var <-c(h_est_obj$running_variance_x,h_est_obj$running_variance_y) +
        (x - prev_mean) * (x - upd_mean)
      h_est_obj$running_variance_x <- upd_var[1]
      h_est_obj$running_variance_y <- upd_var[2]
      x <- (x - upd_mean)/sqrt(upd_var/(h_est_obj$num_obs-1))
    } else {
      prev_mean <- c(h_est_obj$running_mean_x,h_est_obj$running_mean_y)
      upd_mean <- prev_mean*(1-h_est_obj$exp_weight) + x*h_est_obj$exp_weight
      h_est_obj$running_mean_x <- upd_mean[1]
      h_est_obj$running_mean_y <- upd_mean[2]
      if (h_est_obj$num_obs < 2) {
        h_est_obj$running_mean_x <- x[1]
        h_est_obj$running_mean_y <- x[2]
        h_est_obj$running_variance_x <- 1
        h_est_obj$running_variance_y <- 1
        return(h_est_obj)
      }
      upd_var <- (1 - h_est_obj$exp_weight) * 
        c(h_est_obj$running_variance_x,h_est_obj$running_variance_y) + 
        h_est_obj$exp_weight * (x - upd_mean)^2
      h_est_obj$running_variance_x <- upd_var[1]
      h_est_obj$running_variance_y <- upd_var[2]
      x <- (x - upd_mean)/sqrt(upd_var)
    }
  }
  if (any(is.na(x))){
    return(h_est_obj)
  }
  h_x <-
    as.vector(hermite_function_N(h_est_obj$N_param, x[1]))
  h_y <-
    as.vector(hermite_function_N(h_est_obj$N_param, x[2]))
  if (is.na(h_est_obj$exp_weight)) {
    h_est_obj$coeff_vec_x <-
      (h_est_obj$coeff_vec_x * (h_est_obj$num_obs - 1) + h_x) / 
      h_est_obj$num_obs
    h_est_obj$coeff_vec_y <-
      (h_est_obj$coeff_vec_y * (h_est_obj$num_obs - 1) + h_y) / 
      h_est_obj$num_obs
  } else {
    h_est_obj$coeff_vec_x <-
      h_est_obj$coeff_vec_x * (1 - h_est_obj$exp_weight) + h_x * 
      h_est_obj$exp_weight
    h_est_obj$coeff_vec_y <-
      h_est_obj$coeff_vec_y * (1 - h_est_obj$exp_weight) + h_y * 
      h_est_obj$exp_weight
  }
  if (is.na(h_est_obj$exp_weight)){
    h_est_obj$coeff_mat_bivar <- (h_est_obj$coeff_mat_bivar * 
              (h_est_obj$num_obs-1) + tcrossprod(h_x,h_y))/(h_est_obj$num_obs)
  } else {
    h_est_obj$coeff_mat_bivar<- h_est_obj$coeff_mat_bivar * 
      (1-h_est_obj$exp_weight) + tcrossprod(h_x,h_y)*h_est_obj$exp_weight
  }
  return(h_est_obj)
}

#' Updates the Hermite series based estimator sequentially
#'
#' This method can be applied in sequential estimation settings.
#' 
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param x A numeric vector of length 2 or a n x 2 matrix with n bivariate 
#' observations to be incorporated into the estimator.
#' @return An object of class hermite_estimator_bivar.
update_sequential.hermite_estimator_bivar <- function(h_est_obj, x)
{
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (is.null(nrow(x))){
    if (length(x) != 2){
      stop("x must be a vector of length 2 (bivariate observation).")
    } else {
      dim(x) <- c(1,2)
    }
  } else {
    if (ncol(x) != 2 | nrow(x) < 1){
      stop("x must be a matrix with 2 columns (bivariate observations) and
           at least one row.")
    }
  }
  if (any(is.na(x)) | any(!is.finite(x))){
    stop("The sequential update method is only 
         applicable to finite, non NaN, non NA values.")
  }
  for (idx in seq_len(nrow(x))) {
    h_est_obj <- update_sequential_bivar_helper(h_est_obj,
                                     x[idx,])
  }
  return(h_est_obj)
}

#' Initializes the Hermite series based estimator with a batch of data
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param x A numeric matrix. A matrix of bivariate observations to be 
#' incorporated into the estimator. Each row corresponds to a single bivariate
#' observation.
#' @return An object of class hermite_estimator_bivar.
initialize_batch_bivar <- function(h_est_obj, x) {
  h_est_obj$num_obs <- nrow(x)
  if (h_est_obj$standardize_obs == TRUE) {
    h_est_obj$running_mean_x <-mean(x[,1])
    h_est_obj$running_variance_x <-stats::var(x[,1]) * (nrow(x) - 1)
    x[,1] <-
      (x[,1] - h_est_obj$running_mean_x) /
      sqrt(h_est_obj$running_variance_x / (h_est_obj$num_obs - 1))
    h_est_obj$running_mean_y <- mean(x[,2])
    h_est_obj$running_variance_y <- stats::var(x[,2]) * (nrow(x) - 1)
    x[,2] <-
      (x[,2] - h_est_obj$running_mean_y) / sqrt(h_est_obj$running_variance_y /
                                             (h_est_obj$num_obs - 1))
  }
  h_x <-
    hermite_function(h_est_obj$N_param, x[,1])
  h_y <-
    hermite_function(h_est_obj$N_param, x[,2])
  h_est_obj$coeff_vec_x <- rowSums(h_x) / h_est_obj$num_obs
  h_est_obj$coeff_vec_y <- rowSums(h_y) / h_est_obj$num_obs
  h_est_obj$coeff_mat_bivar <- tcrossprod(h_x,h_y) / h_est_obj$num_obs
  return(h_est_obj)
}

# Internal helper method to calculate running standard deviations
calculate_running_std.hermite_estimator_bivar <- function(h_est_obj){
  if (is.na(h_est_obj$exp_weight)) {
    running_std_x <- sqrt(h_est_obj$running_variance_x / 
                            (h_est_obj$num_obs - 1))
    running_std_y <- sqrt(h_est_obj$running_variance_y / 
                            (h_est_obj$num_obs - 1))
  } else {
    running_std_x <- sqrt(h_est_obj$running_variance_x)
    running_std_y <- sqrt(h_est_obj$running_variance_y)
  }
  return(c(running_std_x,running_std_y))
}

# Internal helper method to calculate the probability density at a single 2-d
# x value
dens_helper <- function(h_est_obj, x, clipped)
{
  UseMethod("dens_helper",h_est_obj)
}

dens_helper.hermite_estimator_bivar <- function(h_est_obj,x, clipped = FALSE){
  if (h_est_obj$num_obs < 2) {
    return(NA)
  }
  factor <- 1
  if (h_est_obj$standardize_obs == TRUE) {
    running_std_vec <- calculate_running_std(h_est_obj)
    x <- (x - c(h_est_obj$running_mean_x,h_est_obj$running_mean_y)) / 
      running_std_vec
    factor <- 1 / (prod(running_std_vec))
  }
  return(factor * t(hermite_function(h_est_obj$N_param, x[1])) %*%
           h_est_obj$coeff_mat_bivar %*%
           hermite_function(h_est_obj$N_param, x[2]))
}

#' Estimates the probability densities for a matrix of 2-d x values
#'
#' This method calculates the probability density values for a matrix of 
#' 2-d x vector values using the hermite_estimator_bivar object (h_est_obj).
#'
#' The object must be updated with observations prior to the use of the method.
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param x A numeric matrix. Each row corresponds to a 2-d coordinate.
#' @param clipped A boolean value. This value determines whether
#' probability densities are clipped to be bigger than zero.
#' @param accelerate_series A boolean value. Series acceleration has not yet 
#' been implemented for bivariate estimators.
#' @return A numeric vector of probability density values.
dens.hermite_estimator_bivar <- function(h_est_obj,x, clipped = FALSE, 
                                         accelerate_series = FALSE){
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (is.null(nrow(x))){
    if (length(x) != 2){
      stop("vector input for x must be of length 2.")
    }
    dim(x) <- c(1,2)
  }
  if (ncol(x)!=2 | nrow(x) < 1){
    stop("matrix input for x must have 2 columns and at least 1 row.")
  }
  result <- rep(0,nrow(x))
  for (idx in seq_len(nrow(x))) {
    result[idx] <- dens_helper(h_est_obj,x[idx,])
  }
  if (clipped == TRUE) {
    result <- pmax(result, 1e-08)
  }
  return(result)
}


#' Creates an object summarizing the bivariate PDF with associated generic 
#' methods print and plot.
#'
#' The hermite_estimator_bivar object, x must be updated with observations 
#' prior to the use of the method.
#'
#' @param x A hermite_estimator_bivar object.
#' @param x_lower A numeric vector. This vector determines the lower limit of 
#' x values at which to evaluate the density.
#' @param x_upper A numeric vector. This vector determines the upper limit of 
#' x values at which to evaluate the density.
#' @param ... Additional arguments for the dens function.
#' @return A hdensity_bivar object whose underlying structure is a list 
#' containing the following components.
#' 
#' x: The points at which the density is calculated.
#' x_vals_1: Marginal quantiles of first random variable, used for plotting.
#' x_vals_2: Marginal quantiles of second random variable, used for plotting.
#' density_vals: The density values at the points x.
#' num_obs: The number of observations used to form the Hermite density 
#' estimates.
#' N: The number of terms N in the Hermite series estimator.
#' @export
density.hermite_estimator_bivar <- function(x, x_lower=NA, x_upper=NA, ...) {
  if (x$standardize_obs == FALSE){
    if (any(is.na(x_lower)) | any(is.na(x_upper))){
      stop("For non-standardized hermite_estimator objects, a lower and 
           upper x limits for the PDF summary must be provided i.e. x_lower and 
           x_upper arguments must be provided.")
    }
  } 
  if (all(!(is.na(x_lower))) & all(!is.na(x_upper))){
    x_step_1 <- (x_upper[1] - x_lower[1])/19
    x_vals_1 <- seq(x_lower[1], x_upper[1], by = x_step_1)
    x_step_2 <- (x_upper[2] - x_lower[2])/19
    x_vals_2 <- seq(x_lower[2], x_upper[2], by = x_step_2)
  } else {
    h_est_marginal_x <- hermite_estimator(N=x$N_param, 
                                          exp_weight_lambda = x$exp_weight)
    h_est_marginal_x$coeff_vec <- x$coeff_vec_x
    h_est_marginal_x$num_obs <- x$num_obs
    h_est_marginal_x$running_mean <- x$running_mean_x
    h_est_marginal_x$running_variance <- x$running_variance_x
    h_est_marginal_y <- hermite_estimator(N=x$N_param, 
                                          exp_weight_lambda = x$exp_weight)
    h_est_marginal_y$coeff_vec <- x$coeff_vec_y
    h_est_marginal_y$num_obs <- x$num_obs
    h_est_marginal_y$running_mean <- x$running_mean_y
    h_est_marginal_y$running_variance <- x$running_variance_y
    x_vals_1 <- quant(h_est_marginal_x, p = seq(0.05,0.95,0.05))
    x_vals_2 <- quant(h_est_marginal_y, p = seq(0.05,0.95,0.05))
  }
  x_vals <- cbind(rep(x_vals_1, each = length(x_vals_2)), 
                  rep(x_vals_2, length(x_vals_1)))
  density_vals <- dens(h_est_obj = x, x = x_vals,...)
  result <- list(x = x_vals, x_vals_1 = x_vals_1,
                 x_vals_2 = x_vals_2, 
                 density_vals = density_vals, num_obs = x$num_obs,
                 N = x$N_param)
  class(result) <- c("hdensity_bivar", "list")
  return(result)
}

# Internal helper method to calculate the cumulative probability at a single 
# 2-d x value
cum_prob_helper <- function(h_est_obj, x, clipped = FALSE)
{
  UseMethod("cum_prob_helper",h_est_obj)
}

cum_prob_helper.hermite_estimator_bivar <- function(h_est_obj,x, 
                                                    clipped = FALSE){
  if (h_est_obj$num_obs < 2) {
    return(NA)
  }
  if (h_est_obj$standardize_obs == TRUE) {
    running_std_vec <- calculate_running_std(h_est_obj)
    x <- (x - c(h_est_obj$running_mean_x,h_est_obj$running_mean_y)) / 
      running_std_vec
  }
  return(t(hermite_int_lower(N = h_est_obj$N_param,x = x[1])) 
         %*% h_est_obj$coeff_mat_bivar %*% 
           hermite_int_lower(N = h_est_obj$N_param,x = x[2]))
}

#' Estimates the cumulative probabilities for a matrix of 2-d x values
#'
#' This method calculates the cumulative probability values for a matrix of 
#' 2-d x vector values using the hermite_estimator_bivar object (h_est_obj).
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param x A numeric matrix. Each row corresponds to a 2-d coordinate.
#' @param clipped A boolean value. This value determines whether cumulative
#' probabilities are clipped to lie within the range [0,1].
#' @param accelerate_series A boolean value. Series acceleration has not yet 
#' been implemented for bivariate estimators.
#' @return A numeric vector of cumulative probability values.
cum_prob.hermite_estimator_bivar <- function(h_est_obj,x, clipped = FALSE,
                                             accelerate_series = FALSE){
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (is.null(nrow(x))){
    if (length(x) != 2){
      stop("vector input for x must be of length 2.")
    }
    dim(x) <- c(1,2)
  }
  if (ncol(x)!=2 | nrow(x) < 1){
    stop("matrix input for x must have 2 columns and at least 1 row.")
  }
  result <- rep(0,nrow(x))
  for (idx in seq_len(nrow(x))) {
    result[idx] <- cum_prob_helper(h_est_obj,x[idx,])
  }
  if (clipped == TRUE) {
    result <- pmin(pmax(result, 1e-08), 1)
  }
  return(result)
}

#' Creates an object summarizing the bivariate CDF with associated generic 
#' methods print, plot and summary.
#'
#' The hermite_estimator_bivar object h_est_obj must be updated with 
#' observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param clipped A boolean value. This value determines whether cumulative
#' probabilities are clipped to lie within the range [0,1].
#' @param accelerate_series A boolean value. This value determines whether
#' Hermite series acceleration is applied.
#' @param x_lower A numeric vector. This vector determines the lower limit of 
#' x values at which to evaluate the CDF.
#' @param x_upper A numeric value. This vector determines the upper limit of 
#' x values at which to evaluate the CDF.
#' @return A hcdf_bivar object whose underlying structure is a list 
#' containing the following components.
#' 
#' x: The points at which the cumulative probability is calculated.
#' x_vals_1: Marginal quantiles of first random variable, used for plotting.
#' x_vals_2: Marginal quantiles of second random variable, used for plotting.
#' cum_prob_vals: The cumulative probability values at the points x.
#' num_obs: The number of observations used to form the Hermite cumulative 
#' probability estimates.
#' N: The number of terms N in the Hermite series estimator.
#' @export
hcdf.hermite_estimator_bivar <- function(h_est_obj, clipped = FALSE, 
                                          accelerate_series = TRUE,
                                         x_lower=NA,
                                         x_upper=NA) {
  
  
  if (h_est_obj$standardize_obs == FALSE){
    if (any(is.na(x_lower)) | any(is.na(x_upper))){
      stop("For non-standardized hermite_estimator objects, a lower and 
           upper x limits for the PDF summary must be provided i.e. x_lower and 
           x_upper arguments must be provided.")
    }
  } 
  if (all(!is.na(x_lower)) & all(!is.na(x_upper))){
    x_step_1 <- (x_upper[1] - x_lower[1])/19
    x_vals_1 <- seq(x_lower[1], x_upper[1], by = x_step_1)
    x_step_2 <- (x_upper[2] - x_lower[2])/19
    x_vals_2 <- seq(x_lower[2], x_upper[2], by = x_step_2)
  } else {
    h_est_marginal_x <- hermite_estimator(N=h_est_obj$N_param, 
                                      exp_weight_lambda = h_est_obj$exp_weight)
    h_est_marginal_x$coeff_vec <- h_est_obj$coeff_vec_x
    h_est_marginal_x$num_obs <- h_est_obj$num_obs
    h_est_marginal_x$running_mean <- h_est_obj$running_mean_x
    h_est_marginal_x$running_variance <- h_est_obj$running_variance_x
    h_est_marginal_y <- hermite_estimator(N=h_est_obj$N_param, 
                                      exp_weight_lambda = h_est_obj$exp_weight)
    h_est_marginal_y$coeff_vec <- h_est_obj$coeff_vec_y
    h_est_marginal_y$num_obs <- h_est_obj$num_obs
    h_est_marginal_y$running_mean <- h_est_obj$running_mean_y
    h_est_marginal_y$running_variance <- h_est_obj$running_variance_y
    x_vals_1 <- quant(h_est_marginal_x, p = seq(0.05,0.95,0.05))
    x_vals_2 <- quant(h_est_marginal_y, p = seq(0.05,0.95,0.05))
  }
  x_vals <- cbind(rep(x_vals_1, each = length(x_vals_2)), rep(x_vals_2, 
                                                            length(x_vals_1)))
  cum_prob_vals <- cum_prob(h_est_obj, x_vals, clipped, 
                       accelerate_series)
  
  result <- list(x = x_vals, x_vals_1 = x_vals_1,
                 x_vals_2 = x_vals_2, 
                 cum_prob_vals = cum_prob_vals, num_obs = h_est_obj$num_obs,
                 N = h_est_obj$N_param)
  class(result) <- c("hcdf_bivar", "list")
  return(result)
}

#' Estimates the Spearman's rank correlation coefficient
#'
#' This method calculates the Spearman's rank correlation coefficient value
#' using the hermite_estimator_bivar object (h_est_obj).
#' 
#' The method utilizes the estimator defined in the paper Stephanou, Michael 
#' and Varughese, Melvin. "Sequential Estimation of Nonparametric Correlation 
#' using Hermite Series Estimators." arXiv Preprint (2020), 
#' https://arxiv.org/abs/2012.06287
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param clipped A boolean value. Indicates whether to clip Spearman's rank 
#' correlation estimates to lie between -1 and 1.
#' @return A numeric value.
spearmans.hermite_estimator_bivar <- function(h_est_obj, clipped = FALSE)
{
  if (h_est_obj$num_obs < 2) {
    return(NA)
  }
  W <- W_serialized[1:(h_est_obj$N_param+1),1:(h_est_obj$N_param+1)]
  z <- z_serialized[1:(h_est_obj$N_param+1)]
  W_transpose <- t(W)
  result <- 12*(t(h_est_obj$coeff_vec_x) %*% W_transpose) %*% 
    h_est_obj$coeff_mat_bivar %*% (W %*% h_est_obj$coeff_vec_y) +
    -6 * (t(h_est_obj$coeff_vec_x) %*% W_transpose) %*% 
    (h_est_obj$coeff_mat_bivar %*% z) +
    -6 * (t(z)%*% h_est_obj$coeff_mat_bivar) %*% (W %*% h_est_obj$coeff_vec_y) +
    3 * t(z) %*% (h_est_obj$coeff_mat_bivar%*%z)
  if (clipped == TRUE) {
    result <- pmin(pmax(result, -1), 1)
  }
  return(as.numeric(result))
}

#' Estimates the Kendall rank correlation coefficient
#'
#' This method calculates the Kendall rank correlation coefficient value
#' using the hermite_estimator_bivar object (h_est_obj).
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param clipped A boolean value. Indicates whether to clip the Kendall rank 
#' correlation estimates to lie between -1 and 1.
#' @return A numeric value.
#' @export
kendall.hermite_estimator_bivar <- function(h_est_obj, clipped = FALSE)
{
  if (h_est_obj$num_obs < 2) {
    return(NA)
  }
  W <- W_serialized[1:(h_est_obj$N_param+1),1:(h_est_obj$N_param+1)]
  W_transpose <- t(W)
  result <- 4*sum(diag(W_transpose%*%t(h_est_obj$coeff_mat_bivar) %*% 
                         W%*%h_est_obj$coeff_mat_bivar)) - 1
  if (clipped == TRUE) {
    result <- pmin(pmax(result, -1), 1)
  }
  return(as.numeric(result))
}

#' Prints bivariate hermite_estimator object.
#' 
#'
#' @param x A hermite_estimator_bivar object.
#' @param ... Other arguments passed on to methods used in printing.
print.hermite_estimator_bivar <- function(x, ...) {
  describe_estimator(x,"bivariate")
}

#' Summarizes bivariate hermite_estimator object.
#' 
#' Outputs key parameters of a bivariate hermite_estimator object along with
#' estimates of the mean and standard deviation of the first and second 
#' dimensions of the bivariate data that the object has been updated with.
#' Also outputs the Spearman's Rho and Kendall Tau of the bivariate data that 
#' the object has been updated with.
#'
#' @param object A hermite_estimator_bivar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Other arguments passed on to methods used in summary.
summary.hermite_estimator_bivar <- function(object, 
                              digits = max(3, getOption("digits") - 3), ...) {
  describe_estimator(object,"bivariate")
  if (object$num_obs > 2){
    running_std <- calculate_running_std(object)
    cat("\n")
    cat(paste0("Mean x = ",round(object$running_mean_x,digits), "\n"))
    cat(paste0("Mean y = ",round(object$running_mean_y,digits), "\n"))
    cat(paste0("Standard Deviation x = ", 
               round(running_std[1],digits), "\n"))
    cat(paste0("Standard Deviation y = ", 
               round(running_std[2],digits), "\n"))
    cat(paste0("Spearman's Rho = ",round(spearmans(object),digits), "\n"))
    cat(paste0("Kendall Tau = ",round(kendall(object),digits), "\n"))
  }
}

quant.hermite_estimator_bivar <- 
  function(h_est_obj, p, algorithm, accelerate_series) {
  stop("Quantile estimation is not defined for the bivariate Hermite 
       estimator")
}
