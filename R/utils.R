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
  normalization_hermite <- hermite_normalization(N)
  return(hermite_function(N,x,normalization_hermite))
}

#' Convenience function to output a definite integral of the orthonormal 
#' Hermite functions
#' 
#' The method calculates \eqn{\int_{-\infty}^{x} h_k(t) dt} 
#' for \eqn{k=0,\dots,N} and the vector of values x.
#' 
#' @param N An integer number.
#' @param x A numeric vector.
#' @return A numeric matrix with N+1 rows and length(x) columns.
#' @export
hermite_int_lower <- function(N,x){
  normalization_hermite <- hermite_normalization(N)
  hermite_function_matrix <- hermite_function(N,x,normalization_hermite)
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
#' @return A numeric matrix with N+1 rows and length(x) columns.
#' @export
hermite_int_upper <- function(N,x){
  normalization_hermite <- hermite_normalization(N)
  hermite_function_matrix <- hermite_function(N,x,normalization_hermite)
  return(hermite_integral_val_upper(N,x,hermite_function_matrix))
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

integrand_coeff_univar <- function(t,hermite_est_current, 
                                   hermite_estimator_merged, current_k, 
                                   dimension = NA){
  normalization_hermite <- hermite_est_current$normalization_hermite_vec
  t <- sqrt(2) * t
  herm_mod <- hermite_polynomial(hermite_est_current$N_param, t) *
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
  return(sqrt(2) * hermite_function(current_k,((t*original_sd +original_mean) - 
                                                 new_mean)/new_sd, 
                             normalization_hermite[1:current_k])[current_k, ] * 
           as.vector(crossprod(herm_mod, original_coeff_vec)))
}
