#' Convenience function to output Hermite normalization factors 
#' 
#' The method returns numeric normalization factors that, when multiplied by 
#' the physicist Hermite polynomials times a Gaussian factor i.e.
#' \eqn{\exp{x^2/2}H_k(x)}, yields orthonormal Hermite functions \eqn{h_k(x)} 
#' for \eqn{k=0,\dots,N}.
#'
#' @author Michael Stephanou <michael.stephanou@gmail.com>
#'
#' @param N An integer number.
#' @return A numeric vector of length N+1
#' @export
hermite_normalization_N <- function(N){
if (N < 0) {
  stop("N must be >= 0.")
}
return(hermite_normalization(N))
}

#' Convenience function to output physicist Hermite polynomials
#'  
#'  
#' The method calculates the physicist version of Hermite polynomials, 
#' \eqn{H_k(x)} from \eqn{k=0,\dots,N} for the vector of values, x.
#' 
#' @param N An integer number.
#' @param x A numeric vector.
#' @return A numeric matrix with N+1 rows and length(x) columns.
#' @export
hermite_polynomial_N <- function(N,x){
  if (N < 0) {
    stop("N must be >= 0.")
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x)<1) {
    stop("x must contain at least one value.")
  }
  return(hermite_polynomial(N,x))
}

#' Convenience function to output orthonormal Hermite functions
#'  
#'  
#' The method calculates the orthonormal Hermite functions, \eqn{h_k(x)} 
#' from \eqn{k=0,\dots,N} for the vector of values, x.
#' 
#' @param N An integer number.
#' @param x A numeric vector.
#' @return A numeric matrix with N+1 rows and length(x) columns.
#' @export
hermite_function_N <- function(N,x){
  if (N < 0) {
    stop("N must be >= 0.")
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x)<1) {
    stop("x must contain at least one value.")
  }
  return(hermite_function(N,x))
}

#' Convenience function to output the sum of orthonormal Hermite functions
#'  
#'  
#' The method calculates the sum of orthonormal Hermite functions, 
#' \eqn{\sum_{i} h_k(x_{i})} from \eqn{k=0,\dots,N} for the vector of values,
#' x.
#' 
#'
#' @param N An integer number.
#' @param x A numeric vector.
#' @return A numeric vector of length N+1.
#' @export
hermite_function_sum_N <- function(N,x){
  if (N < 0) {
    stop("N must be >= 0.")
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x)<1) {
    stop("x must contain at least one value.")
  }
  if (N < 2){
    return(hermite_function_sum_serial(N,x))
  } else {
    if (getOption("hermiter.parallel", TRUE) == TRUE){
      return(hermite_function_sum_parallel(N,x))
    } else {
      return(hermite_function_sum_serial(N,x))
    }
  }
}

#' Convenience function to output a definite integral of the orthonormal 
#' Hermite functions
#' 
#' The method calculates \eqn{\int_{-\infty}^{x} h_k(t) dt} 
#' for \eqn{k=0,\dots,N} and the vector of values x.
#' 
#' @param N An integer number.
#' @param x A numeric vector.
#' @param hermite_function_matrix A numeric matrix. A matrix of Hermite 
#' function values. 
#' @return A numeric matrix with N+1 rows and length(x) columns.
#' @export
hermite_int_lower <- function(N,x,hermite_function_matrix=NULL){
  if (N < 0) {
    stop("N must be >= 0.")
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x)<1) {
    stop("x must contain at least one value.")
  }
  if (is.null(hermite_function_matrix)){
    hermite_function_matrix <- hermite_function_N(N,x)
  } else {
    if (nrow(hermite_function_matrix) != (N+1)){
      stop("Hermite function matrix must have N+1 rows.")
    }
  }
  return(hermite_integral_val(N,x,hermite_function_matrix))
}

#' Convenience function to output a definite integral of the orthonormal 
#' Hermite functions
#' 
#' The method calculates \eqn{\int_{x}^{\infty} h_k(t) dt} 
#' for \eqn{k=0,\dots,N} and the vector of values x.
#' 
#' @param N An integer number.
#' @param x A numeric vector.
#' @param hermite_function_matrix A numeric matrix. A matrix of Hermite 
#' function values.
#' @return A numeric matrix with N+1 rows and length(x) columns.
#' @export
hermite_int_upper <- function(N,x, hermite_function_matrix=NULL){
  if (N < 0) {
    stop("N must be >= 0.")
  }
  if (!is.numeric(x)) {
    stop("x must be numeric.")
  }
  if (length(x)<1) {
    stop("x must contain at least one value.")
  }
  if (is.null(hermite_function_matrix)){
    hermite_function_matrix <- hermite_function_N(N,x)
  } else {
    if (nrow(hermite_function_matrix) != (N+1)){
      stop("Hermite function matrix must have N+1 rows.")
    }
  }
  return(hermite_integral_val_upper(N,x,hermite_function_matrix))
}

#' Convenience function to output the integral of the orthonormal Hermite 
#' functions on the full domain
#' 
#' The method calculates \eqn{\int_{-\infty}^{\infty} h_k(t) dt} 
#' for \eqn{k=0,\dots,N}.
#' 
#' @param N An integer number.
#' @return A numeric matrix with N+1 rows and 1 columns.
#' @export
hermite_int_full <- function(N){
  if (N < 0) {
    stop("N must be >= 0.")
  }
  return(hermite_int_full_domain(N))
}

#' Calculates \eqn{\int_{-\infty}^{\infty} f(x) e^{-x^2} dx} using 
#' Gauss-Hermite quadrature with 100 terms.
#'  
#' @param f A function.
#' @return A numeric value.
#' @export
gauss_hermite_quad_100 <- function(f){
  result <- sum(f(root_x_serialized)*weight_w_serialized)
  return(result)
}

# Helper method for merging univariate hermite_estimator objects.
integrand_coeff_univar <- function(t,hermite_est_current, 
                                   hermite_estimator_merged, current_k, 
                                   dimension = NA){
  normalization_hermite <- hermite_est_current$normalization_hermite_vec
  t <- sqrt(2) * t
  herm_mod <- hermite_polynomial_N(hermite_est_current$N_param, t) *
    normalization_hermite
  if (is.na(dimension)){
    original_sd <- sqrt(hermite_est_current$running_variance / 
                          (hermite_est_current$num_obs-1))
    original_mean <- hermite_est_current$running_mean
    new_sd <- sqrt(hermite_estimator_merged$running_variance / 
                     (hermite_estimator_merged$num_obs-1))
    new_mean <- hermite_estimator_merged$running_mean
    original_coeff_vec <- hermite_est_current$coeff_vec
  } else {
    if (dimension==1){
      original_sd <- sqrt(hermite_est_current$running_variance_x / 
                            (hermite_est_current$num_obs-1))
      original_mean <- hermite_est_current$running_mean_x
      new_sd <- sqrt(hermite_estimator_merged$running_variance_x / 
                       (hermite_estimator_merged$num_obs-1))
      new_mean <- hermite_estimator_merged$running_mean_x
      original_coeff_vec <- hermite_est_current$coeff_vec_x
    } else if (dimension==2) {
      original_sd <- sqrt(hermite_est_current$running_variance_y / 
                            (hermite_est_current$num_obs-1))
      original_mean <- hermite_est_current$running_mean_y
      new_sd <- sqrt(hermite_estimator_merged$running_variance_y / 
                       (hermite_estimator_merged$num_obs-1))
      new_mean <- hermite_estimator_merged$running_mean_y
      original_coeff_vec <- hermite_est_current$coeff_vec_y
    } 
  }
  return(sqrt(2) * hermite_function_N(current_k-1,((t*original_sd + 
                                            original_mean) -  new_mean)/new_sd) 
                             [current_k, ] * 
           as.vector(crossprod(herm_mod, original_coeff_vec)))
}

# Helper method for estimating Hermite series sums, with or without series 
# acceleration. These series acceleration techniques are drawn from the below 
# reference:
# 
# Boyd, John P., and Dennis W. Moore. "Summability methods for 
# Hermite functions." Dynamics of atmospheres and oceans 10.1 (1986): 51-62. 
# 
# Note that h_input is a numeric matrix of N+1 rows and length(x) columns,
# coeffs is a numeric vector of length N+1
#
# This helper method is intended for internal use by the 
# hermite_estimator_univar class.
series_calculate <- function(h_input, coeffs, accelerate_series = TRUE){
  N <- nrow(h_input) - 1
  if (length(coeffs) < 3 | accelerate_series == FALSE){
    return(as.numeric(crossprod(h_input,coeffs)))
  }
  if (length(coeffs) >=3 & length(coeffs) < 6){
    result <- crossprod(h_input[1:N,,drop=FALSE], coeffs[1:N])
    result <- result + 1/2 * h_input[N+1,] * coeffs[N+1]
    return(as.numeric(result))
  }
  if (length(coeffs) >=6 & length(coeffs) < 12){
    result <- crossprod(h_input[1:(N-3),,drop=FALSE], coeffs[1:(N-3)])
    result <- result + 15/16*h_input[N-2,]*coeffs[N-2]+
      11/16*h_input[N-1,]*coeffs[N-1]+
      5/16*h_input[N,]*coeffs[N]+
      1/16*h_input[N+1,]*coeffs[N+1]
    return(as.numeric(result))
  }
  if (length(coeffs) >= 12){
    result <- crossprod(h_input[1:(N-7),,drop=FALSE], coeffs[1:(N-7)])
    result <- result +
    255/256*h_input[N-6,]*coeffs[N-6]+
    247/256*h_input[N-5,]*coeffs[N-5]+
    219/256*h_input[N-4,]*coeffs[N-4]+
    163/256*h_input[N-3,]*coeffs[N-3]+
    93/256*h_input[N-2,]*coeffs[N-2]+
    37/256*h_input[N-1,]*coeffs[N-1]+
    9/256*h_input[N,]*coeffs[N]+
    1/256*h_input[N+1,]*coeffs[N+1]
    return(as.numeric(result))
  }
}
