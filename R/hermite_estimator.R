#' A class to sequentially estimate univariate and bivariate pdfs and cdfs along
#' with quantile functions in the univariate setting and nonparametric 
#' correlations in the bivariate setting.
#'
#' The hermite_estimator class provides a unified interface to the univariate 
#' and bivariate Hermite series based estimators, leveraging generic methods and
#' S3 dispatch. Methods are included for the sequential or one-pass batch 
#' estimation of the full probability density function and cumulative 
#' distribution function in the univariate and bivariate settings. Sequential 
#' or one-pass batch estimation methods are also provided for the full quantile 
#' function in the univariate setting along with the Spearman and Kendall 
#' correlation coefficients in the bivariate setting. Note that RcppParallel is
#' utilized to speed up batch updating in the univariate case. If one wishes to 
#' switch to serial batch updating (typically slower), utilize 
#' options(hermiter.parallel = FALSE).
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
#' @param est_type A string value. Options are "univariate" or "bivariate".
#' @param observations A numeric vector or a numeric matrix. Note that for 
#' univariate estimators, x is a numeric vector of observations to be 
#' incorporated. For bivariate estimators, x is a numeric matrix with n rows for
#' n observations and 2 columns.
#' @return An S3 object of class hermite_estimator_univar or 
#' hermite_estimator_bivar. 
#' @export
#' @examples
#' \dontrun{
#' hermite_est <- hermite_estimator(N = 50, standardize = TRUE,
#' est_type="univariate")
#' hermite_est <- hermite_estimator(N = 50, standardize = TRUE,
#' est_type="univariate", observations = c(1,2,3))
#' hermite_est <- hermite_estimator(N = 30, standardize = TRUE, 
#' est_type="bivariate", observations = matrix(c(1,1,2,2,3,3), 
#' nrow=3, ncol=2,byrow=TRUE))
#' }
hermite_estimator <-
  function(N = NA,
           standardize = TRUE,
           exp_weight_lambda = NA, est_type = "univariate", observations = c()){
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
    if (is.na(N)){
      if (est_type=="univariate"){
        if (is.na(exp_weight_lambda)){
          return(hermite_estimator_univar(N = 50,standardize,exp_weight_lambda,
                                        observations))
        } else {
          return(hermite_estimator_univar(N = 20,standardize,exp_weight_lambda,
                                          observations))
        }
      }
      else if (est_type=="bivariate"){
        if (is.na(exp_weight_lambda)){
          return(hermite_estimator_bivar(N = 30,standardize,exp_weight_lambda,
                                       observations))
        } else {
          return(hermite_estimator_bivar(N = 20,standardize,exp_weight_lambda,
                                          observations))
        }
      } else {
        stop("Unknown estimator type.")
      }
    } else {
      if (!is.numeric(N) | is.nan(N)) {
        stop("N must be numeric.")
      }
      if (N < 0 | N > 75) {
        stop("N must be >= 0 and N <= 75.")
      }
      if (est_type=="univariate"){
        return(hermite_estimator_univar(N, standardize,exp_weight_lambda,
                                        observations))
      }
      else if (est_type=="bivariate"){
        return(hermite_estimator_bivar(N, standardize,exp_weight_lambda,
                                       observations))
      } else {
        stop("Unknown estimator type.")
      }
    }
  }

#' Merges two Hermite estimators
#'
#' Note that the estimators must be of the same type to be merged i.e. both 
#' estimators must have a consistent est_type, either "univariate" or 
#' "bivariate". In addition, the N and standardize arguments must be the same 
#' for both estimators in order to merge them. Finally, note that exponentially 
#' weighted estimators cannot be merged. If the Hermite estimators are not 
#' standardized, the merged estimator will be exactly equivalent to constructing
#' a single estimator on the data set formed by combining the data sets used to 
#' update the respective hermite_estimator inputs. If the input Hermite 
#' estimators are standardized however, then the equivalence will be approximate
#' but still reasonably accurate in most cases.
#'
#' @param h_est_obj A hermite_estimator_univar or hermite_estimator_bivar 
#' object. 
#' The first Hermite series based estimator.
#' @param hermite_estimator_other A hermite_estimator_univar or
#' hermite_estimator_bivar object. The second Hermite series based estimator.
#' @return An object of class hermite_estimator_univar or 
#' hermite_estimator_bivar.
#' @export
#' @examples
#' \dontrun{
#' hermite_est_1 <- hermite_estimator(N = 10, standardize = FALSE, 
#' observations = rnorm(30))
#' hermite_est_2 <- hermite_estimator(N = 10, standardize = FALSE, 
#' observations = rnorm(30))
#' hermite_merged <- merge_pair(hermite_est_1, hermite_est_2)
#' }
merge_pair <- function(h_est_obj, hermite_estimator_other) {
  UseMethod("merge_pair", h_est_obj)
}

#' Merges a list of Hermite estimators
#'
#' Note that the estimators must be of the same type to be merged i.e. all 
#' estimators must have a consistent est_type, either "univariate" or 
#' "bivariate". In addition, the N and standardize arguments must be the same 
#' for all estimators in order to merge them. Finally, note that exponentially 
#' weighted estimators cannot be merged. If the Hermite estimators are not 
#' standardized, the merged estimator will be exactly equivalent to constructing
#' a single estimator on the data set formed by combining the data sets used to 
#' update the respective hermite_estimator inputs. If the input Hermite 
#' estimators are standardized however, then the equivalence will be approximate
#' but still reasonably accurate in most cases.
#'
#' @param hermite_estimators A list of hermite_estimator_univar or 
#' hermite_estimator_bivar objects. 
#' @return An object of class hermite_estimator_univar or 
#' hermite_estimator_bivar.
#' @export
#' @examples
#' \dontrun{
#' hermite_est_1 <- hermite_estimator(N = 10, standardize = FALSE, 
#' observations = rnorm(30))
#' hermite_est_2 <- hermite_estimator(N = 10, standardize = FALSE, 
#' observations = rnorm(30))
#' hermite_merged <- merge_hermite(list(hermite_est_1, hermite_est_2))
#' }
merge_hermite <- function(hermite_estimators) {
  UseMethod("merge_hermite", hermite_estimators)
}

merge_hermite.list <- function(hermite_estimators){
  if (length(hermite_estimators) == 0) {
    stop("List must contain at least one Hermite estimator.")
  }
  if (length(hermite_estimators) == 1) {
    return(hermite_estimators[[1]])
  }
  all_classes <- lapply(hermite_estimators, FUN =
                          function(x){return(class(x)[[1]])})
  if (length(unique(all_classes)) >1) {
    stop("List must contain Hermite estimators of a consistent type")
  }
  if (class(hermite_estimators[[1]])[1]=="hermite_estimator_univar"){
    return(merge_hermite_univar(hermite_estimators))
  } else {
    return(merge_hermite_bivar(hermite_estimators))
  }
}

#' Updates the Hermite series based estimator sequentially
#'
#' This method can be applied in sequential estimation settings.
#' 
#'
#' @param h_est_obj A hermite_estimator_univar or hermite_estimator_bivar 
#' object.
#' @param x A numeric vector or matrix. Observations to be incorporated into 
#' the estimator. Note that for univariate estimators, x is a numeric vector 
#' whereas for bivariate estimators, x is a numeric vector of length 2 or a 
#' n x 2 matrix with n bivariate observations to be incorporated into the 
#' estimator.
#' @return An object of class hermite_estimator_univar or 
#' hermite_estimator_bivar.
#' @export
#' @examples
#' \dontrun{
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
#' est_type="univariate")
#' hermite_est <- update_sequential(hermite_est, x = 2)
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
#' est_type="bivariate")
#' hermite_est <- update_sequential(hermite_est, x = c(1,2))
#' }
update_sequential <- function(h_est_obj, x) {
  UseMethod("update_sequential", h_est_obj)
}

# An internal method to calculate running standard deviation from the scaled
# running variance.
calculate_running_std <- function(h_est_obj)
{
  UseMethod("calculate_running_std",h_est_obj)
}

#' Estimates the probability density at one or more x values
#'
#' This method calculates the probability density values at a vector of
#' x values in the univariate case. In the bivariate case, the method calculates
#' the probability density values for a matrix of x values, each row of which 
#' represents a 2-d point. 
#'
#' The object must be updated with observations prior to the use of the method.
#'
#' @param h_est_obj A hermite_estimator_univar or hermite_estimator_bivar 
#' object.
#' @param x A numeric vector (univariate) or a numeric matrix (bivariate) of
#' values at which to calculate the probability density.
#' @param clipped A boolean value. This value determines whether
#' probability densities are clipped to be bigger than zero.
#' @param accelerate_series A boolean value. This value determines whether
#' Hermite series acceleration is applied.
#' @return A numeric vector of probability density values.
#' @export
#' @examples
#' \dontrun{
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
#' est_type="univariate", observations = rnorm(30))
#' pdf_est <- dens(hermite_est, c(0, 0.5, 1))
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
#' est_type="bivariate", observations = matrix(rnorm(60), 
#' nrow=30, ncol=2,byrow=TRUE))
#' pdf_est <- dens(hermite_est, matrix(c(0,0,0.5,0.5,1,1),nrow=3,
#' ncol=2,byrow=TRUE))
#' }
dens <- function(h_est_obj, x, clipped, accelerate_series = TRUE) {
  UseMethod("dens", h_est_obj)
}

#' Estimates the cumulative probability at one or more x values
#'
#' This method calculates the cumulative probability at a vector of
#' x values in the univariate case. In the bivariate case, the method calculates
#' the probability density values for a matrix of x values, each row of which 
#' represents a 2-d point. 
#'
#' The object must be updated with observations prior to the use of the method.
#'
#' @param h_est_obj A hermite_estimator_univar or hermite_estimator_bivar 
#' object.
#' @param x A numeric vector (univariate) or a numeric matrix (bivariate).
#' Values at which to calculate the cumulative probability.
#' @param clipped A boolean value. This value determines whether
#' cumulative probabilities are clipped to lie between 0 and 1.
#' @param accelerate_series A boolean value. This value determines whether
#' Hermite series acceleration is applied.
#' @return A numeric vector of cumulative probability values.
#' @export
#' @examples
#' \dontrun{
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
#' est_type="univariate", observations = rnorm(30))
#' cdf_est <- cum_prob(hermite_est, c(0, 0.5, 1))
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
#' est_type="bivariate", observations = matrix(rnorm(60), 
#' nrow=30, ncol=2,byrow=TRUE))
#' cdf_est <- cum_prob(hermite_est, matrix(c(0,0,0.5,0.5,1,1),nrow=3,
#' ncol=2,byrow=TRUE))
#' }
cum_prob <- function(h_est_obj, x, clipped, accelerate_series = TRUE) {
  UseMethod("cum_prob", h_est_obj)
}

#' Creates an object summarizing the CDF with associated generic methods print,
#' plot and summary.
#'
#' The h_est_obj object must be updated with observations prior to the use of 
#' the method.
#'
#' @param h_est_obj A hermite_estimator_univar or hermite_estimator_bivar 
#' object.
#' @param clipped A boolean value. This value determines whether
#' cumulative probabilities are clipped to lie between 0 and 1.
#' @param accelerate_series A boolean value. This value determines whether
#' Hermite series acceleration is applied.
#' @param x_lower A numeric value (univariate) or a numeric vector (bivariate).
#' This value determines the lower limit of x values at which to evaluate 
#' the CDF.
#' @param x_upper A numeric value (univariate) or a numeric vector (bivariate).
#' This value determines the upper limit of x values at which to evaluate 
#' the CDF.
#' @return A hcdf_univar or hcdf_bivar object.
hcdf <- function(h_est_obj, clipped = FALSE, accelerate_series = TRUE, 
                 x_lower = NA, x_upper = NA) {
  UseMethod("hcdf", h_est_obj)
}

#' Estimates the quantiles at a vector of probability values
#' 
#' This method utilizes the estimator (13) in paper Stephanou, Michael, 
#' Varughese, Melvin and Iain Macdonald. "Sequential quantiles via Hermite 
#' series density estimation." Electronic Journal of Statistics 11.1 (2017): 
#' 570-607 <doi:10.1214/17-EJS1245>, with some modifications to improve the 
#' stability of numerical root finding when using the bisection algorithm. 
#' Note that this method is only applicable to the univariate Hermite 
#' estimator i.e. est_type = "univariate".
#' 
#' The object must be updated with observations prior to the use of this method.
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
#' @export
#' @examples
#' \dontrun{
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE, 
#' est_type="univariate", observations = rnorm(30))
#' quant_est <- quant(hermite_est, c(0.25, 0.5, 0.75))
#' }
quant <- function(h_est_obj, p, algorithm="interpolate", 
                  accelerate_series = TRUE) {
  UseMethod("quant", h_est_obj)
}

#' Estimates the Spearman's rank correlation coefficient
#' 
#' This method utilizes the estimator (8) in the paper Stephanou, Michael and 
#' Varughese, Melvin. "Sequential estimation of Spearman rank correlation using 
#' Hermite series estimators." Journal of Multivariate Analysis (2021) 
#' <doi:10.1016/j.jmva.2021.104783> to calculate the Spearman rank correlation 
#' coefficient. It is only applicable to the bivariate Hermite estimator 
#' i.e. est_type = "bivariate".
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param clipped A boolean value. Indicates whether to clip Spearman's rank 
#' correlation estimates to lie between -1 and 1.
#' @return A numeric value.
#' @export
#' @examples
#' \dontrun{
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE,
#' est_type="bivariate", observations = matrix(rnorm(30*2), nrow=30, 
#' ncol=2, byrow = TRUE))
#' spearmans_est <- spearmans(hermite_est)
#' }
spearmans <- function(h_est_obj, clipped = FALSE)
{
  UseMethod("spearmans",h_est_obj)
}

#' Estimates the Kendall rank correlation coefficient
#'
#' This method calculates the Kendall rank correlation coefficient value. It 
#' is only applicable to the bivariate Hermite estimator i.e. est_type = 
#' "bivariate".
#'
#' The object must be updated with observations prior to the use of this method.
#'
#' @param h_est_obj A hermite_estimator_bivar object.
#' @param clipped A boolean value. Indicates whether to clip the Kendall rank 
#' correlation estimates to lie between -1 and 1.
#' @return A numeric value.
#' @export
#' @examples
#' \dontrun{
#' hermite_est <- hermite_estimator(N = 10, standardize = TRUE,
#' est_type="bivariate", observations = matrix(rnorm(30*2), nrow=30, 
#' ncol=2, byrow = TRUE))
#' kendall_est <- kendall(hermite_est)
#' }
kendall <- function(h_est_obj, clipped = FALSE)
{
  UseMethod("kendall",h_est_obj)
}

# Internal helper function for print methods.
describe_estimator <- function(h_est_obj, est_type){
  if (est_type == "univariate") {
    cat("Univariate Hermite Estimator:\n")
  } else {
    cat("Bivariate Hermite Estimator:\n")
  }
  cat(paste0("N = ",h_est_obj$N_param,"\n"))
  cat(paste0("Standardize observations = ",h_est_obj$standardize_obs,"\n"))
  if (!is.na(h_est_obj$exp_weight)){
    cat(paste0("Exponential weighting for coefficents = TRUE, lambda = ", 
               h_est_obj$exp_weight,"\n"))
  } else {
    cat(paste0("Exponential weighting for coefficents = FALSE","\n"))
  }
  cat(paste0("Number of observations = ",h_est_obj$num_obs,"\n"))
}
