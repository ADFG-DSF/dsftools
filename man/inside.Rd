% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inside_subset.R
\name{inside}
\alias{inside}
\alias{\%inside\%}
\alias{\%inside()\%}
\alias{\%inside(]\%}
\alias{\%inside[)\%}
\title{Inside}
\usage{
x \%inside\% a

x \%inside()\% a

x \%inside(]\% a

x \%inside[)\% a
}
\arguments{
\item{x}{Vector input}

\item{a}{Vector expressing the range of values.  The user does not need to
manually supply the min and max of the range, as the min and max will be
automatically calculated from this vector.}
}
\value{
Logical vector of the same length as argument \code{x}
}
\description{
Similar to \code{\%in\%}, but for numeric input.  This function will
return a logical vector with length equal to the first operand, with values
of \code{TRUE} if \code{x} has a numeric value within the range specified by \code{a}.

The square and round brackets in the function names denote closed, open, or
clopen intervals; that is, whether the endpoints are included.  See examples
below.

Note: for consistency with \code{\%in\%}, NA values in \code{x} are currently given a values
of \code{FALSE}.
}
\examples{
xx <- 1:10
aa <- 1:3

## can manually supply interval min and max as a vector
xx \%inside\% c(1, 3)

## or can use a vector as input
xx \%inside\% aa

## differences between interval closure:
xx \%inside()\% aa
xx \%inside[)\% aa
xx \%inside(]\% aa

## handling of NA is similar to \%in\%
c(xx, NA) \%inside\% aa
}
\author{
Matt Tyers
}
