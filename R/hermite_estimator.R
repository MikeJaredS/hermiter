hermite_estimator <-
  function(N = 10,
           standardize = FALSE,
           exp_weight_lambda = NA, est_type = "univariate") {
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
    if (est_type=="univariate"){
      return(hermite_estimator_univar(N,standardize,exp_weight_lambda))
    }
    else if (est_type=="bivariate"){
      return(hermite_estimator_bivar(N,standardize,exp_weight_lambda))
    } else {
      stop("Unknown estimator type.")
    }
  }

#' @export
combine_hermite <- function(hermite_estimators) {
  UseMethod("combine_hermite", hermite_estimators)
}

combine_hermite.list <- function(hermite_estimators){
  if (length(hermite_estimators) == 0) {
    stop("List must contain at least one Hermite estimator.")
  }
  if (length(hermite_estimators) == 1) {
    return(hermite_estimators[[1]])
  }
  # all_classes <- lapply(hermite_estimators, FUN= function(x){return(class(x))})
  if (class(hermite_estimators[[1]])=="hermite_estimator_univar"){
    return(combine_hermite_univar(hermite_estimators))
  } else {
    return(combine_hermite_bivar(hermite_estimators))
  }
}

#' @export
combine_pair <- function(this, hermite_estimator_other) {
  UseMethod("combine_pair", this)
}

#' @export
update_sequential <- function(this, x) {
  UseMethod("update_sequential", this)
}

#' @export
update_batch <- function(this, x) {
  UseMethod("update_batch", this)
}

calculate_running_std <- function(this)
{
  UseMethod("calculate_running_std",this)
}

#' @export
dens <- function(this, x, clipped) {
  UseMethod("dens", this)
}

#' @export
cum_prob <- function(this, x, clipped) {
  UseMethod("cum_prob", this)
}

#' @export
quant <- function(this, p) {
  UseMethod("quant", this)
}

#' @export
spearmans <- function(this, clipped = FALSE)
{
  UseMethod("spearmans",this)
}