#' Logit
#' @description Logit, defined as
#' \eqn{log\left(\frac{x}{1-x}\right)}
#' @param x Numeric vector
#' @return Numeric vector
#' @seealso \link{expit}
#' @author Matt Tyers
#' @examples
#' logit(0.5)
#'
#' x <- seq(from=0.01, to=0.99, by=0.01)
#' plot(x, logit(x))
#' @export
logit <- function(x) log(x/(1-x))


#' Expit, or inverse logit
#' @description Inverse logit, where logit is defined as \eqn{log\left(\frac{x}{1-x}\right)}.
#'
#' Expit (inverse logit) is defined as \eqn{\frac{exp(x)}{1+exp(x)}}
#' @param x Numeric vector
#' @return Numeric vector
#' @seealso \link{logit}
#' @author Matt Tyers
#' @examples
#' expit(0)
#'
#' x <- seq(from=-5, to=5, by=0.1)
#' plot(x, expit(x))
#' @export
expit <- function(x) exp(x)/(1+exp(x))



#' Standard Error
#' @description Standard error, defined as
#' \eqn{\frac{\sigma_x}{\sqrt{n}}}
#' @param x Numeric vector
#' @param na.rm Logical.  Should missing values be removed?  Defaults to `FALSE`
#' for consistency with \link{sd}.  Note: if `na.rm==TRUE`, the denominator will
#' be the number of non-NA entries.
#' @return Numeric of length 1
#' @author Matt Tyers
#' @examples
#' a <- c(8, 6, 7, 5, 3, 0, 9, NA)
#'
#' se(a)
#' se(a, na.rm=TRUE)
#' sd(a, na.rm=TRUE)/sqrt(7)
#' @export
se <- function(x, na.rm=FALSE) {
  sd(x, na.rm=na.rm)/sqrt(sum(!is.na(x)))
}
