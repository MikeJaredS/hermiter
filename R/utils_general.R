#' Prints the hcdf_univar object as output by the hcdf function when evaluated
#' on a hermite_estimator_univar object.
#' 
#' Mirrors the print method of the stats::ecdf function 
#'
#' @param x A hcdf_univar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Unused
#' @export
print.hcdf_univar <- function (x, 
                               digits = getOption("digits") - 2L, ...) 
{
  numform <- function(x) paste(formatC(x, digits = digits), 
                               collapse = ", ")
  n <- length(x$x)
  i1 <- 1L:min(3L, n)
  i2 <- if (n >= 4L) 
    max(4L, n - 1L):n
  else integer()
  cat(" x[1:", n, "] = ", numform(x$x[i1]), if (n > 3L) 
    ", ", if (n > 5L) 
      " ..., ", numform(x$x[i2]), "\n", sep = "")
  cat(" p[1:", n, "] = ", numform(x$cum_prob_vals[i1]), if (n > 3L) 
    ", ", if (n > 5L) 
      " ..., ", numform(x$cum_prob_vals[i2]), "\n", sep = "")
  invisible(x)
}


#' Plots the hcdf_univar object as output by the hcdf function when evaluated
#' on a hermite_estimator_univar object.
#'
#' @param x A hcdf_univar object.
#' @param main A string, title for plot.
#' @param xlab A string, x label for plot.
#' @param ylab A string, y label for plot.
#' @param ... Additional parameters for plotting.
#' @export
plot.hcdf_univar <- function (x, main="Hermite CDF", xlab = "x", 
                              ylab = "F(x)", ...) 
{
  plot.default(x=x$x,y=x$cum_prob_vals, main=main, xlab=xlab, 
               ylab = ylab,type="l",...)
  abline(h = c(0, 1), lty = 2)
}

#' Summarizes the hcdf_univar object as output by the hcdf function when 
#' evaluated on a hermite_estimator_univar object.
#'
#' @param object A hcdf_univar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Unused.
#' @export
summary.hcdf_univar <- function (object, 
                                 digits = getOption("digits") - 2L, ...) 
{
  cat("Hermite CDF estimates:\n")
  cat(object$num_obs,"observations\n")
  n <- length(object$x)
  cat(n,"evaluation points\n")
  print.hcdf_univar(object, 
                    digits = getOption("digits") - 2L, ...)
  invisible(NULL)
}

#' Prints the hcdf_bivar object as output by the hcdf function when evaluated
#' on a hermite_estimator_bivar object.
#' 
#'
#' @param x A hcdf_bivar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Additional parameters for printing.
#' @export
print.hcdf_bivar <- function (x, 
                              digits = getOption("digits") - 2L, ...) 
{
  n <- length(x$cum_prob_vals)
  if (n==4){
    peak_head <- cbind(x$x[1:4,,drop=FALSE],x$cum_prob_vals[1:4])
    rownames(peak_head) <- 1:4
    colnames(peak_head) <- c("x1","x2","p")
    print.default(peak_head, digits=digits,...)
    return(invisible(x))
  } else {
    i1 <- 1L:min(3L, n)
    i2 <- if (n >= 4L) 
      max(4L, n - 1L):n
    else integer()
  }
  peak_head <- cbind(x$x[i1,,drop=FALSE],x$cum_prob_vals[i1])
  rownames(peak_head) <- i1
  colnames(peak_head) <- c("x1","x2","p")
  print.default(peak_head, digits=digits,...)
  if (length(i2) > 0){
    peak_tail <- cbind(x$x[i2,,drop=FALSE],x$cum_prob_vals[i2])
    rownames(peak_tail) <- i2
    colnames(peak_tail) <- rep("",3)
    if (n > 5L) {cat("...\n")}
    print.default(peak_tail, digits=digits,...)  
  }
  invisible(x)
}

#' Plots the hcdf_bivar object as output by the hcdf function when evaluated
#' on a hermite_estimator_bivar object.
#'
#' @param x A hcdf_bivar object.
#' @param main A string, title for plot.
#' @param xlab A string, x label for plot.
#' @param ylab A string, y label for plot.
#' @param ... Unused.
#' @export
plot.hcdf_bivar <- function (x, main="Hermite CDF",xlab = "X", ylab="Y",...) 
{
  z <- matrix(x$cum_prob_vals,byrow = FALSE,nrow=length(x$x_vals_1))
  filled.contour(x = x$x_vals_1,y = x$x_vals_2, z= z,
                 color.palette = function(n) hcl.colors(n, "Oslo", rev = TRUE),
                 main=main,xlab=xlab,ylab=ylab)
}

#' Summarizes the hcdf_bivar object as output by the hcdf function when 
#' evaluated on a hermite_estimator_bivar object.
#'
#' @param object A hcdf_bivar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Additional parameters for printing.
#' @export
summary.hcdf_bivar <- function (object, 
                                digits = getOption("digits") - 2L, ...) 
{
  cat("Hermite CDF estimates:\n")
  cat(object$num_obs,"observations\n")
  n <- nrow(object$x)
  cat(n,"evaluation points\n")
  print.hcdf_bivar(object,digits = getOption("digits") - 2L, ...)
  invisible(NULL)
}

#' Prints the hdensity_univar object as output by the density function when
#' evaluated on a hermite_estimator_univar object.
#' 
#'
#' @param x A hdensity_univar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Unused
#' @export
print.hdensity_univar <- function (x, 
                                   digits = getOption("digits") - 2L, ...) 
{
  numform <- function(x) paste(formatC(x, digits = digits), 
                               collapse = ", ")
  n <- length(x$x)
  i1 <- 1L:min(3L, n)
  i2 <- if (n >= 4L) 
    max(4L, n - 1L):n
  else integer()
  cat(" x[1:", n, "] = ", numform(x$x[i1]), if (n > 3L) 
    ", ", if (n > 5L) 
      " ..., ", numform(x$x[i2]), "\n", sep = "")
  cat(" d[1:", n, "] = ", numform(x$density_vals[i1]), if (n > 3L) 
    ", ", if (n > 5L) 
      " ..., ", numform(x$density_vals[i2]), "\n", sep = "")
  invisible(x)
}

#' Plots the hdensity_univar object as output by the density function when 
#' evaluated on a hermite_estimator_univar object.
#'
#' @param x A hdensity_univar object.
#' @param main A string, title for plot.
#' @param xlab A string, x label for plot.
#' @param ylab A string, y label for plot.
#' @param ... Additional parameters for plotting.
#' @export
plot.hdensity_univar <- function (x, main="Hermite PDF", xlab = "x", 
                                  ylab = "Density",...) 
{
  plot.default(x=x$x,y=x$density_vals, main=main, xlab=xlab, 
               ylab = ylab,type="l",...)
  abline(h = 0, lty = 2)
}

#' Prints the hdensity_bivar object as output by the density function when
#' evaluated on a hermite_estimator_bivar object.
#' 
#'
#' @param x A hdensity_bivar object.
#' @param digits A numeric value. Number of digits to round to.
#' @param ... Additional parameters for printing.
#' @export
print.hdensity_bivar <- function (x, 
                                  digits = getOption("digits") - 2L, ...) 
{
  n <- length(x$density_vals)
  if (n==4){
    peak_head <- cbind(x$x[1:4,,drop=FALSE],x$density_vals[1:4])
    rownames(peak_head) <- 1:4
    colnames(peak_head) <- c("x1","x2","d")
    print.default(peak_head, digits=digits,...)
    return(invisible(x))
  } else {
    i1 <- 1L:min(3L, n)
    i2 <- if (n >= 4L) 
      max(4L, n - 1L):n
    else integer()
  }
  peak_head <- cbind(x$x[i1,,drop=FALSE],x$density_vals[i1])
  rownames(peak_head) <- i1
  colnames(peak_head) <- c("x1","x2","d")
  print.default(peak_head, digits=digits,...)
  if (length(i2) > 0){
    peak_tail <- cbind(x$x[i2,,drop=FALSE],x$density_vals[i2])
    rownames(peak_tail) <- i2
    colnames(peak_tail) <- rep("",3)
    if (n > 5L) {cat("...\n")}
    print.default(peak_tail, digits=digits,...)  
  }
  invisible(x)
}

#' Plots the hdensity_bivar object as output by the density function when 
#' evaluated on a hermite_estimator_bivar object.
#'
#' @param x A hdensity_bivar object.
#' @param main A string, title for plot.
#' @param xlab A string, x label for plot.
#' @param ylab A string, y label for plot.
#' @param ... Unused.
#' @export
plot.hdensity_bivar <- function(x, main="Hermite PDF",xlab = "X", ylab = "Y", 
                                ...) 
{
  z <- matrix(x$density_vals,byrow = FALSE,nrow=length(x$x_vals_1))
  cols <- hcl.colors(10, "Oslo")
  filled.contour(x = x$x_vals_1,y = x$x_vals_2, z= z,
                 color.palette = function(n) hcl.colors(n, "Oslo", rev = TRUE),
                 main=main, xlab = xlab, ylab = ylab)
}

#' A wrapper around the stats::cor function adding two additional methods, 
#' namely method = "hermite.spearman" and method = "hermite.kendall" (can be
#' abbreviated). The input parameters and output value semantics closely match 
#' the stats::cor method for easy interchange. If neither the 
#' "hermite.spearman" nor the "hermite.kendall" method is selected, then this 
#' function will call stats::cor with the arguments provided.
#'
#' @param x a numeric vector, matrix or data frame. 
#' @param y NULL (default) or a vector, matrix or data frame with 
#' compatible dimensions to x. The default is equivalent to y = x 
#' (but more efficient).
#' @param use not used by hermite.spearman and hermite.kendall methods. For 
#' stats::cor this is an optional character string giving a method for 
#' computing covariances in the presence of missing values. 
#' This must be (an abbreviation of) one of the strings "everything", 
#' "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method a character string indicating which correlation coefficient 
#'  is to be computed. One of "pearson" (default), "kendall", "spearman", 
#'  "hermite.spearman" or "hermite.kendall": can be abbreviated.
#' @param ... Additional arguments for the constructor of the hermite_estimator
#' object if method is "hermite.spearman" or "hermite.kendall".
#' @export
cor <- function(x, y= NULL, use = "everything", method="pearson", ...){
  if (!startsWith(method,"hermite.k") & !startsWith(method,"hermite.s")){
    return(stats::cor(x,y,use,method))
  }
  # The error trapping below draws from the stats::cor function
  if (is.data.frame(y)) 
    y <- as.matrix(y)
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y)) 
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y))) 
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
    if (length(x) == 0L || length(y) == 0L) 
      stop("both 'x' and 'y' must be non-empty")
    matrix_result <- is.matrix(x) || is.matrix(y)
    if (!is.matrix(x)) 
      x <- matrix(x, ncol = 1L)
    if (!is.matrix(y)) 
      y <- matrix(y, ncol = 1L)
  }
  if (startsWith(method,"hermite.k")){
    result_mat <- matrix(rep(NA,ncol(x)^2),nrow=ncol(x))
    for (i in seq_len(ncol(x))) {
      for (j in seq_len(ncol(x))) {
        if (!is.null(y)){
          h_est <- hermite_estimator(est_type = "bivariate",
                                     observations = cbind(x[,i], y[,j], 
                                                        deparse.level = 0),...)
          result_mat[i,j] <- kendall(h_est)
        } else {
          if (i == j){
            result_mat[i,i] <- 1
          } else if (j > i) {
            h_est <- hermite_estimator(est_type = "bivariate",
                                       observations = x[,c(i,j)],...)
            correl_element <- kendall(h_est)
            result_mat[i,j] <- correl_element
            result_mat[j,i] <- correl_element
          }
        }
      }
    }
    if (!is.null(y)){
      rownames(result_mat) <- colnames(x)
      colnames(result_mat) <- colnames(y)
      if (matrix_result) 
        return(result_mat)
      else return(drop(result_mat))
    } else {
      rownames(result_mat) <- colnames(x)
      colnames(result_mat) <- colnames(x)
    }
    return(result_mat)
  } else if (startsWith(method,"hermite.s")){
    result_mat <- matrix(rep(NA,ncol(x)^2),nrow=ncol(x))
    for (i in seq_len(ncol(x))) {
      for (j in seq_len(ncol(x))) {
        if (!is.null(y)){
          h_est <- hermite_estimator(est_type = "bivariate",
                                     observations = cbind(x[,i], y[,j], 
                                                          deparse.level = 0))
          result_mat[i,j] <- spearmans(h_est)
        } else {
          if (i == j){
            result_mat[i,i] <- 1
          } else if (j > i) {
            h_est <- hermite_estimator(est_type = "bivariate",
                                       observations = x[,c(i,j)], ...)
            correl_element <- spearmans(h_est)
            result_mat[i,j] <- correl_element
            result_mat[j,i] <- correl_element
          }
        }
      }
    }
    if (!is.null(y)){
      rownames(result_mat) <- colnames(x)
      colnames(result_mat) <- colnames(y)
      if (matrix_result) 
        return(result_mat)
      else return(drop(result_mat))
    } else {
      rownames(result_mat) <- colnames(x)
      colnames(result_mat) <- colnames(x)
    }
    return(result_mat)
  }
}

#' Estimates the Interquartile range (IQR)
#' 
#' This creates a default generic method for the stats::IQR function.
#'
#' @param x A numeric vector.
#' @param ... Optional additional arguments to the stats::IQR function.
#' @return A numeric value.
#' @export
IQR.default <- function(x, ...) stats::IQR(x, ...)

#' Estimates the Interquartile range (IQR)
#' 
#' This generic method dispatches to the stats::IQR function or the 
#' IQR.hermite_estimator_univar function depending on the class of x.
#'
#' @param x A numeric vector or hermite_estimator_univar object.
#' @param ... Optional additional arguments.
#' @return A numeric value.
IQR <- function(x, ...) UseMethod("IQR")
