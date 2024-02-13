#' Logit
#' @description Logit, defined as
#' \deqn{log\left(\frac{x}{1-x}\right)}
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
#' Expit (inverse logit) is defined as \deqn{\frac{exp(x)}{1+exp(x)}}
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

