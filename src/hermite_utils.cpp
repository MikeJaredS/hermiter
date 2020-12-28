// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/erf.hpp>

using namespace Rcpp;

//' Outputs physicist version of Hermite Polynomials
//' 
//' 
//' The method calculates the physicist version of Hermite polynomials, 
//' \eqn{H_k(x)} from \eqn{k=0,\dots,N} for the vector of values, x.
//' 
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//'
//' @param N An integer number.
//' @param x A numeric vector.
//' @return A numeric matrix with N+1 rows and length(x) columns.
//' @export
// [[Rcpp::export]]
NumericMatrix hermite_polynomial(int N, NumericVector x) {
  int x_size = x.size();
  NumericMatrix hermite(N + 1, x_size);
  for(int i = 0; i < x_size; ++i) {
    hermite(0,i) = 1;
  }
  if (N==0){
    return hermite;
  }
  for(int i = 0; i < x_size; ++i) {
    hermite(1,i) = 2 * x[i];
  }
  if (N==1){
    return hermite;
  }
  for(int j = 0; j < x_size; ++j) {
    for(int i = 2; i <= N; ++i) {
      hermite(i,j) = 2 * x[j] * hermite(i - 1,j) - 2 *
       (((double)i) - 1) * hermite(i - 2,j);
    }
  }
  return hermite;
}

//' Outputs Hermite normalization factors 
//' 
//' The method returns numeric normalization factors that, when multiplied by 
//' the physicist Hermite polynomials \eqn{H_k(x)}, yield orthonormal 
//' Hermite functions \eqn{h_k(x)} for \eqn{k=0,\dots,N}.
//'
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//'
//' @param N An integer number.
//' @return A numeric vector of length N+1
//' @export
// [[Rcpp::export]]
NumericVector hermite_normalization(int N) {
  NumericVector out(N + 1);
  double sqrt_pi = sqrt(M_PI);
  for(int i = 0; i <= N; ++i) {
    out(i) = 1/sqrt(pow((double)2,(double)i) * 
      boost::math::factorial<double>((double)i) * sqrt_pi); 
  }
  return out;
}

//' Outputs orthonormal Hermite functions
//' 
//' 
//' The method calculates the orthonormal Hermite functions, \eqn{h_k(x)} 
//' from \eqn{k=0,\dots,N} for the vector of values, x.
//' 
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//'
//' @param N An integer number.
//' @param x A numeric vector.
//' @param normalization A numeric vector of normalization factors generated by
//' the hermite_normalization function.
//' @return A numeric matrix with N+1 rows and length(x) columns.
//' @export
// [[Rcpp::export]]
NumericMatrix hermite_function(int N, NumericVector x, 
                               NumericVector normalization) {
  int x_size = x.size();
  NumericMatrix hermite(N + 1, x_size);
  NumericMatrix out(N + 1, x_size);
  NumericVector expFac(x_size);
  for(int i = 0; i < x_size; ++i) {
    hermite(0,i) = 1;
    expFac(i) = exp(-1 * x[i] * x[i] / 2);
    out(0,i) = hermite(0,i) * normalization[0] * expFac(i);
  }
  if (N==0){
    return out;
  }
  for(int i = 0; i < x_size; ++i) {
    hermite(1,i) = 2 * x[i];
    out(1,i) = hermite(1,i) * normalization[1] * expFac(i);
  }
  if (N==1){
    return out;
  }
  for(int j = 0; j < x_size; ++j) {
    for(int i = 2; i <= N; ++i) {
      hermite(i,j) = 2 * x[j] * hermite(i - 1,j) - 2 *
       (((double)i) - 1) * hermite(i - 2,j);
      out(i,j) = hermite(i,j) * normalization[i] * expFac(j);
    }
  }
  return out;
}

//' Outputs lower integral of the orthonormal Hermite functions
//' 
//' The method calculates \eqn{\int_{-\infty}^{x} h_k(t) dt} 
//' for \eqn{k=0,\dots,N} and the vector of values x.
//' 
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//'
//' @param N An integer number.
//' @param x A numeric vector.
//' @param hermite_function_mat A numeric matrix of Hermite function values 
//' generated by the function hermite_function.
//' @return A numeric matrix with N+1 rows and length(x) columns.
//' @export
// [[Rcpp::export]]
NumericMatrix hermite_integral_val(int N, NumericVector x, 
                                   NumericMatrix hermite_function_mat) {
  int x_size = x.size();
  NumericMatrix out(N + 1, x_size);
  for(int i = 0; i < x_size; ++i) {
    out(0,i) = pow(M_PI,0.25)/sqrt((double) 2) *
     boost::math::erfc<double>((double) -1 * (1 / sqrt((double)2)) * x[i]); 
  }
  if (N == 0){
    return out;
  }
  for(int i = 0; i < x_size; ++i) {
    out(1,i) = -1 * sqrt((double) 2)/pow(M_PI,0.25) *exp(-1 * x[i] * x[i] / 2); 
  }
  if (N == 1){
    return out;
  }
  for(int i = 2; i <= N; ++i) {
    for(int j = 0; j < x_size; ++j) {
      out(i,j) = -1 * sqrt(2/((double)i)) * hermite_function_mat(i - 1,j) + 
       sqrt((((double)i) - 1)/((double)i)) * out(i - 2,j);
    }
  }
  return out;
}

//' Outputs upper integral of the orthonormal Hermite functions
//' 
//' The method calculates \eqn{\int_{x}^{\infty} h_k(t) dt} 
//' for \eqn{k=0,\dots,N} and the vector of values x.
//' 
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//'
//' @param N An integer number.
//' @param x A numeric vector.
//' @param hermite_function_mat A numeric matrix of Hermite function values 
//' generated by the function hermite_function.
//' @return A numeric matrix with N+1 rows and length(x) columns.
//' @export
// [[Rcpp::export]]
NumericMatrix hermite_integral_val_upper(int N,
                         NumericVector x, NumericMatrix hermite_function_mat) {
  int x_size = x.size();
  NumericMatrix out(N + 1, x_size);
  for(int i = 0; i < x_size; ++i) {
    out(0,i) = pow(M_PI,0.25)/sqrt((double) 2) *
     boost::math::erfc<double>((double) (1 / sqrt((double) 2)) * x[i]); 
  }
  if (N == 0){
    return out;
  }
  for(int i = 0; i < x_size; ++i) {
    out(1,i) = sqrt((double) 2)/pow(M_PI,0.25) *exp(-1 * x[i] * x[i] / 2); 
  }
  if (N == 1){
    return out;
  }
  for(int i = 2; i <= N; ++i) {
    for(int j = 0; j < x_size; ++j) {
      out(i,j) = sqrt(((double)2) / ((double)i)) * hermite_function_mat(i-1,j)+ 
       sqrt((((double)i)-1)/((double)i)) * out(i-2,j);
    }
  }
  return out;
}

//' Outputs integral of the orthonormal Hermite functions on the full domain
//' 
//' The method calculates \eqn{\int_{-\infty}^{\infty} h_k(t) dt} 
//' for \eqn{k=0,\dots,N}.
//' 
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//'
//' @param N An integer number.
//' @return A numeric matrix with N+1 rows and 1 columns.
//' @export
// [[Rcpp::export]]
NumericMatrix hermite_int_full_domain(int N) {
  NumericMatrix out(N + 1, 1);
  out(0,0) = pow(M_PI,0.25)*sqrt((double) 2); 
  if (N == 0){
    return out;
  }
  out(1,0) = 0; 
  if (N == 1){
    return out;
  }
  for(int i = 2; i <= N; ++i) {
      out(i,0) = sqrt((((double)i)-1)/((double)i)) * out(i-2,0);
  }
  return out;
}

//' Standardizes the observation x and updates the online moment inputs
//' 
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//' 
//' @param x A numeric value.
//' @param n_obs A numeric value. The number of observations.
//' @param current_mean A numeric value. 
//' @param current_var A numeric value. 
//' @return A numeric vector. The first element is the updated mean. The
//' second element is the updated variance times n_obs. The third element is the
//' updated, standardized value of x.
//' @export
// [[Rcpp::export]]
NumericVector standardizeInputs(double x, double n_obs, double current_mean,
                                double current_var) {
  NumericVector outputVec(3);
  double prev_running_mean = current_mean;
  outputVec[0] =  (current_mean * (n_obs - 1) + x) / n_obs;
  if (n_obs < 2){
    return outputVec;
  }
  outputVec[1] = current_var + (x - prev_running_mean) * (x - outputVec[0]);
  double running_std = sqrt(outputVec[1] / (n_obs - 1));
  outputVec[2] = (x - outputVec[0]) / running_std;
  return outputVec;
}

//' Standardizes the observation x and updates the online moment inputs
//' 
//' The online moments are updated via exponential weighting.
//' 
//' @author Michael Stephanou <michael.stephanou@gmail.com>
//'
//' @param x A numeric value.
//' @param n_obs A numeric value. The number of observations.
//' @param lambda A numeric value.
//' @param current_mean A numeric value.
//' @param current_var A numeric value.  
//' @return A numeric vector. The first element is the updated mean. The
//' second element is the updated variance times n_obs. The third element is the
//' updated, standardized value of x.
//' @export
// [[Rcpp::export]]
NumericVector standardizeInputsEW(double x, double n_obs,double lambda,
                                  double current_mean, double current_var) {
  NumericVector outputVec(3);
  if (n_obs < 2){
    current_mean = x;
    current_var = 1;
    outputVec[0] =  current_mean;
    outputVec[1] = current_var;
    return outputVec;
  }
  outputVec[0] =  (1 - lambda) * current_mean + lambda * x;
  current_var = (1 - lambda) * current_var + lambda * pow((x - outputVec[0]),2);
  outputVec[1] = current_var;
  double running_std = sqrt(current_var);
  outputVec[2] = (x - outputVec[0]) / running_std;
  return outputVec;
}
