% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{hermite_int_upper}
\alias{hermite_int_upper}
\title{Convenience function to output a definite integral of the orthonormal 
Hermite functions}
\usage{
hermite_int_upper(
  N,
  x,
  hermite_function_matrix = NULL,
  normalization_hermite = NULL
)
}
\arguments{
\item{N}{An integer number.}

\item{x}{A numeric vector.}

\item{hermite_function_matrix}{A numeric matrix. A matrix of Hermite 
function values.}

\item{normalization_hermite}{A numeric vector. A vector of normalization 
values necessary in the calculation of the Hermite functions.}
}
\value{
A numeric matrix with N+1 rows and length(x) columns.
}
\description{
The method calculates \eqn{\int_{x}^{\infty} h_k(t) dt} 
for \eqn{k=0,\dots,N} and the vector of values x.
}
