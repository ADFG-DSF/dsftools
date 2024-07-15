#' Relative Precision from Desired Confidence or Accuracy, from simulated values
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



#' Detection Probability From a Given Sample Size, Assuming Binomial
#' @description In the context of a telemetry study, this function calculates
#' the probability of detecting a SINGLE area used by some proportion of the marked
#' population, given the sample size of instrumented fish, assuming a Binomial
#' distribution with random sampling.
#'
#' This is actually a simple wrapper of the function call below, and re-written
#' for convenience:
#'
#' `pbinom(q = observe_at_least-1,`
#'
#' `       size = round(n_raw*assumed_survival),`
#'
#' `       prob = prop_usedby,`
#'
#' `       lower.tail = FALSE)`
#'
#' Because of this, vectors may be used in arguments, but care should be taken
#' to ensure the function evaluates as desired (see examples below).
#'
#' If the probability of simultaneous detection among multiple areas is desired,
#' \link{multinomial_detection} may be considered.
#' @param n_raw Trial value of sample size, before accounting for mortality.
#' @param prop_usedby Hypothetical proportion of the population using the area
#' considered.
#' @param assumed_survival Assumed survival (or 1 - data loss proportion).
#' Defaults to `1`.
#' @param observe_at_least Minimum number of marked individuals to consider as
#' detection of an area.  Defaults to `1`, but a larger number may be used as
#' necessary, depending on criterion used to define detection of an aggregation.
#' @return A single number reflecting the probability of detection.
#' @seealso [multinomial_detection]
#' @author Matt Tyers
#' @importFrom stats pbinom rmultinom
#' @examples
#' ## The probability of detecting an area used by 5% of the population, given
#' ## a sample size of n=100 and assuming 80% survival.
#' binomial_detection(n_raw=100, prop_usedby = 0.05, assumed_survival = .8)
#'
#' ## examples with vector-valued input:
#' binomial_detection(n_raw=c(80, 100), prop_usedby = 0.05, assumed_survival = .8)
#' binomial_detection(n_raw=100, prop_usedby = c(0.025, 0.05), assumed_survival = .8)
#'
#' ## if multiple vectors are used, evaluation will be strictly vector-wise
#' binomial_detection(n_raw=c(80, 100), prop_usedby = c(0.025, 0.05))
#'
#' ## but outer() may be used if a 2d matrix is desired?
#' n_raw_trial <- c(80, 100)
#' prop_usedby_trial <- c(0.025, 0.05)
#' names(n_raw_trial) <- n_raw_trial              ## optional but useful
#' names(prop_usedby_trial) <- prop_usedby_trial  ## optional but useful
#' outer(n_raw_trial, prop_usedby_trial, FUN=binomial_detection)
#' outer(n_raw_trial, prop_usedby_trial, FUN=binomial_detection, assumed_survival=0.8)
#'
#' ## using outer() with different arguments
#' n_raw_trial <- c(80, 100)
#' survival_trial <- c(1, 0.9, 0.8)
#' outer(X=n_raw_trial, Y=survival_trial,
#'       FUN=function(X,Y) binomial_detection(n_raw=X,
#'                                            prop_usedby=0.05,
#'                                            assumed_survival=Y))
#' @export
binomial_detection <- function(n_raw, prop_usedby, assumed_survival=1, observe_at_least=1) {
  pbinom(q = observe_at_least-1,
         size = round(n_raw*assumed_survival),
         prob = prop_usedby,
         lower.tail = FALSE)
}

# # binomial_detection(135, 0.05, .8)

# binomial_detection(n_raw=80, prop_usedby = 0.05, assumed_survival = .8)
# binomial_detection(n_raw=80, prop_usedby = 0.05, assumed_survival = .8, observe_at_least = 2)
# binomial_detection(n_raw=80, prop_usedby = 0.025, assumed_survival = .8)
# binomial_detection(n_raw=80, prop_usedby = 0.025, assumed_survival = .8, observe_at_least = 2)
# binomial_detection(n_raw=100, prop_usedby = 0.05, assumed_survival = .8)
# binomial_detection(n_raw=100, prop_usedby = 0.05, assumed_survival = .8, observe_at_least = 2)
# binomial_detection(n_raw=100, prop_usedby = 0.025, assumed_survival = .8)
# binomial_detection(n_raw=100, prop_usedby = 0.025, assumed_survival = .8, observe_at_least = 2)




#' Simultaneous Detection Probability From a Given Sample Size, Assuming Multinomial
#' @description In the context of a telemetry study, this function uses simulation
#' to estimate the probability of detecting a single area used by some proportion
#' of the marked population, given the sample size of instrumented fish, as well
#' as the probability of detecting MULTIPLE areas, assuming random sampling from
#' a Multinomial worst-case scenario.
#'
#' The multinomial worst-case scenario is defined as a distribution with equal
#' probabilities of the proportion specified by `prop_usedby=`.
#' For example, if the proportion supplied is 5%, a multinomial distribution with
#' twenty (=1/0.05) categories with equal probabilities will be simulated.
#'
#' If the probability of detection of just a single area is desired,
#' \link{binomial_detection} may be considered as an exact solution.
#'
#' @param n_raw Trial value of sample size, before accounting for mortality.
#' @param prop_usedby Hypothetical proportion of the population using the area
#' considered.
#' @param assumed_survival Assumed survival (or 1 - data loss proportion).
#' Defaults to `1`.
#' @param observe_at_least Minimum number of marked individuals to consider as
#' detection of an area.  Defaults to `1`, but a larger number may be used as
#' necessary, depending on criterion used to define detection of an aggregation.
#' @param prop_ofareas Proportion of areas to simultaneously detect.  It may be
#' desirable to structure a precision statement in terms of detection of some
#' percentage of areas, see examples below.  Defaults to `1`, which can be
#' interpreted as simultaneous detection of all such areas.
#' @return Named list of estimated probabilities:
#'
#' * Element `$avg_p_detected` gives the average estimated detection probability
#' over all simulated areas, which will be a good estimate of the detection
#' probability of a SINGLE given area.
#'
#' * Element `$p_all_detected` gives the estimated SIMULTANEOUS detection probability
#' of the proportion of areas specified.
#' @seealso [binomial_detection]
#' @author Matt Tyers
#' @examples
#' ## The probability of detecting all areas used by 5% of the population, given
#' ## a sample size of n=80 and assuming 80% survival, is at least 44%.
#' multinomial_detection(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8)
#'
#'
#' ## The probability of detecting 90% of areas used by 5% of the population, given
#' ## a sample size of n=80 and assuming 80% survival, is at least 97%.
#' multinomial_detection(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8,
#'                       prop_ofareas = 0.9)
#' @export
multinomial_detection <- function(n_raw, prop_usedby, assumed_survival=0.8, observe_at_least=1,
                                  prop_ofareas=1) {
  if(length(n_raw) > 1 |
     length(prop_usedby) > 1 |
     length(assumed_survival) > 1 |
     length(observe_at_least) > 1 |
     length(prop_ofareas) > 1) {
    stop("Evaluation of inputs with length > 1 is not currently implemented.")
  }

  nsim <- 100000
  p <- rep(1, round(1/prop_usedby))
  observed <- rmultinom(n=nsim, size=round(n_raw*assumed_survival), prob=p) >= observe_at_least
  p_detected <- mean(observed)
  p_all_detected <- mean(colMeans(observed) >= prop_ofareas)
  return(list(avg_p_detected=p_detected, p_all_detected=p_all_detected))
}
# multinomial_detection(n_raw=80, prop_usedby = 0.05, assumed_survival = .8)
# multinomial_detection(n_raw=80, prop_usedby = 0.05, assumed_survival = .8, observe_at_least = 2)
# multinomial_detection(n_raw=80, prop_usedby = 0.025, assumed_survival = .8)
# multinomial_detection(n_raw=80, prop_usedby = 0.025, assumed_survival = .8, observe_at_least = 2)
# multinomial_detection(n_raw=100, prop_usedby = 0.05, assumed_survival = .8)
# multinomial_detection(n_raw=100, prop_usedby = 0.05, assumed_survival = .8, observe_at_least = 2)
# multinomial_detection(n_raw=100, prop_usedby = 0.025, assumed_survival = .8)
# multinomial_detection(n_raw=100, prop_usedby = 0.025, assumed_survival = .8, observe_at_least = 2)
