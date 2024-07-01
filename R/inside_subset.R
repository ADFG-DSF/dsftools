#' Inside
#' @name inside
#' @rdname inside
#' @description Similar to `%in%`, but for numeric input.  This function will
#' return a logical vector with length equal to the first operand, with values
#' of `TRUE` if `x` has a numeric value within the range specified by `a`.
#'
#' The square and round brackets in the function names denote closed, open, or
#' clopen intervals; that is, whether the endpoints are included.  See examples
#' below.
#'
#' Note: for consistency with `%in%`, NA values in `x` are currently given a values
#' of `FALSE`.
#' @param x Vector input
#' @param a Vector expressing the range of values.  The user does not need to
#' manually supply the min and max of the range, as the min and max will be
#' automatically calculated from this vector.
#' @return Logical vector of the same length as argument `x`
#' @author Matt Tyers
#' @examples
#' xx <- 1:10
#' aa <- 1:3
#'
#' ## can manually supply interval min and max as a vector
#' xx %inside% c(1, 3)
#'
#' ## or can use a vector as input
#' xx %inside% aa
#'
#' ## differences between interval closure:
#' xx %inside()% aa
#' xx %inside[)% aa
#' xx %inside(]% aa
#'
#' ## handling of NA is similar to %in%
#' c(xx, NA) %inside% aa
NULL

#' @rdname inside
#' @export
`%inside%` <- function(x, a) {
  out1 <- (x >= min(a, na.rm=TRUE)) & (x <= max(a, na.rm=TRUE))
  out1[is.na(out1)] <- FALSE
  return(out1)
}

#' @rdname inside
#' @export
`%inside()%` <- function(x, a) {
  out1 <- (x > min(a, na.rm=TRUE)) & (x < max(a, na.rm=TRUE))
  out1[is.na(out1)] <- FALSE
  return(out1)
}

#' @rdname inside
#' @export
`%inside(]%` <- function(x, a) {
  out1 <- (x > min(a, na.rm=TRUE)) & (x <= max(a, na.rm=TRUE))
  out1[is.na(out1)] <- FALSE
  return(out1)
}

#' @rdname inside
#' @export
`%inside[)%` <- function(x, a) {
  out1 <- (x >= min(a, na.rm=TRUE)) & (x < max(a, na.rm=TRUE))
  out1[is.na(out1)] <- FALSE
  return(out1)
}




#' Vector subsetting shorthand
#' @name subsets
#' @rdname subsets
#' @description Shortcut functions returning a subset of some vector.  These
#' may be handy when dealing with vectors or data.frame columns with long names.
#'
#' `%s_l%`, `%s_leq%`, `%s_g%`, and `%s_geq%` denote the subset less than, less
#' than or equal to, greater than, and greater than or equal to some quantity.
#'
#' `x %s_l% a` is evaluated simply as `x[x < a]`, preserving all vector-wise
#' behavior if it is desirable to exploit this.
#'
#' `%s_inside%` and friends denote the subset bounded by the range expressed by
#' some vector (see \link{inside}).
#' @param x Vector input
#' @param a Number or vector for logical comparison
#' @return Subset of vector `x`
#' @author Matt Tyers
#' @examples
#' ## subsets in comparison to a single number
#' 1:10 %s_l% 5
#' 1:10 %s_leq% 5
#' 1:10 %s_g% 5
#' 1:10 %s_geq% 5
#'
#' ## subsets in comparison to a range
#' 1:10 %s_inside% 3:7
#' 1:10 %s_inside()% 3:7
#' 1:10 %s_inside[)% 3:7
#' 1:10 %s_inside(]% 3:7
NULL


#' @rdname subsets
#' @export
`%s_l%` <- function(x, a) x[x < a]

#' @rdname subsets
#' @export
`%s_g%` <- function(x, a) x[x > a]

#' @rdname subsets
#' @export
`%s_leq%` <- function(x, a) x[x <= a]

#' @rdname subsets
#' @export
`%s_geq%` <- function(x, a) x[x >= a]

#' @rdname subsets
#' @export
`%s_inside%` <- function(x, a) x[x %inside% a]

#' @rdname subsets
#' @export
`%s_inside()%` <- function(x, a) x[x %inside()% a]

#' @rdname subsets
#' @export
`%s_inside(]%` <- function(x, a) x[x %inside(]% a]

#' @rdname subsets
#' @export
`%s_inside[)%` <- function(x, a) x[x %inside[)% a]



# s <- function(x, FUN, a, na.rm=FALSE) {
#   logi1 <- FUN(x, a)
#   if(na.rm) logi1[is.na(x)] <- FALSE
#   return(x[logi1])
# }

