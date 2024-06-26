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



#' Relative Precision from Desired Confidence or Accuracy
#' @description Calculates (relative) accuracy from desired confidence, or
#' confidence from desired (relative) accuracy, from a vector of simulated values
#' and the true value.
#'
#' This function is intended to fill in the typical language in an Operational
#' Plan, stating "...estimate (some parameter), such that the estimate is within
#' (`accuracy`) of the true value (`confidence`) percent of the time."
#'
#' This function can calculate (relative) accuracy from desired confidence, or
#' confidence from desired (relative) accuracy, depending on whether the `accuracy=`
#' argument is specified.  See examples below for both usage cases.
#'
#' Depending on the use case, it may be more appropriate to calculate (or
#' calculate from) relative accuracy, which is defined as
#' \eqn{|\frac{estimate - truth}{truth}|}
#' or absolute accuracy, which is defined as
#' \eqn{|estimate - truth|}.
#' @param sim_vec Numeric vector of simulations
#' @param true_val Assumed true value, which may be given as a single number or
#' corresponding vector if needed
#' @param confidence Vector of desired confidence levels.  Defaults to
#' `c(0.8, 0.9, 0.95, 0.99)`.
#' @param accuracy Vector of desired (relative) accuracy.  Defaults to `NA`, which
#' will calculate (relative) accuracy from precision.  If this argument is given a
#' numeric input, the function will calculate precision from (relative) accuracy.
#' @param relative Whether to calculate (or calculate from) RELATIVE accuracy or
#' ABSOLUTE accuracy.  Defaults to `TRUE`, which will specify relative accuracy.
#' @return Named vector of values of (relative) accuracy or confidence, depending
#' on usage.
#' @author Matt Tyers
#' @examples
#' ## first simulating a random vector
#' xx <- rnorm(10000, mean=10, sd=1)
#'
#' ## calculating RELATIVE accuracy from desired confidence
#' rp(xx, true_val=10)
#'
#' ## calculating confidence from desired RELATIVE accuracy
#' rp(xx, true_val=10, accuracy=seq(from=0.1, to=0.5, by=0.1))
#'
#' ## calculating ABSOLUTE accuracy from desired confidence
#' rp(xx, true_val=10, relative=FALSE)
#'
#' ## calculating confidence from desired ABSOLUTE accuracy
#' rp(xx, true_val=10, accuracy=seq(from=0.5, to=2.5, by=0.5), relative=FALSE)
#'
#' ## plotting both methods and showing agreement
#' plot(x = seq(.01, .5, by=.01),
#'      y = rp(xx, true_val=10, accuracy=seq(.01, .5, by=.01)),
#'      xlab = "relative accuracy", ylab="confidence")
#' points(x = rp(xx, true_val=10, confidence=c(0.8, 0.9, 0.95, 0.99)),
#'        y = c(0.8, 0.9, 0.95, 0.99), pch=16, col=2)
#' @importFrom stats quantile
#' @export
rp <- function(sim_vec, true_val,
               confidence = c(0.8, 0.9, 0.95, 0.99), accuracy = NA,
               relative = TRUE) {
  if(!is.numeric(sim_vec)) stop("Non-numeric input to sim_vec=")
  if(!is.numeric(true_val)) stop("Non-numeric input to true_val=")
  if(relative) {
    comparison_vec <- abs((sim_vec - true_val)/true_val)
  } else {
    comparison_vec <- abs(sim_vec - true_val)
  }
  if(all(is.na(accuracy))) {
    # if relative accuracy is calculated from confidence

    if(!is.numeric(confidence)) stop("Non-numeric input to confidence=")
    out <- quantile(comparison_vec, p=confidence)

  } else {
    # if confidence is calculated from relative accuracy

    if(!is.numeric(accuracy)) stop("Non-numeric input to accuracy=")
    out <- colMeans(outer(comparison_vec, accuracy, FUN="<="))
    names(out) <- accuracy
  }
  return(out)
}
