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



#' Detection Probability From a Given Sample Size, Assuming Binomial (Deprecated)
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
#' @seealso [detection_probability]
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
  .Deprecated("detection_probability")
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




#' Simultaneous Detection Probability From a Given Sample Size, Assuming Multinomial (Deprecated)
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
#' @seealso [detection_probability]
#' @author Matt Tyers
#' @examples
#' ## The probability of detecting all areas used by 5% of the population, given
#' ## a sample size of n=80 and assuming 80% survival, is at least 44%.
#' multinomial_detection(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8)
#'
#'
#' ## The probability of detecting 90% of areas used by 5% of the population, given
#' ## a sample size of n=80 and assuming 80% survival, is at least 82.8%.
#' multinomial_detection(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8,
#'                       prop_ofareas = 0.9)
#' @export
multinomial_detection <- function(n_raw, prop_usedby, assumed_survival=1, observe_at_least=1,
                                  prop_ofareas=1) {
  .Deprecated("detection_probability")
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


#' Detection Probability From a Given Sample Size, Assuming Random Sampling
#' @description In the context of a telemetry study, this function estimates the
#' probability of detecting a single area used by some proportion
#' of the marked population, given the sample size of instrumented fish, and
#' optionally the probability of detecting MULTIPLE areas, assuming random sampling from
#' a worst-case scenario.
#'
#' The function can use either the Binomial or Poisson probability models, and
#' the probability of detecting a given (single) use area is calculated from
#'
#' `pbinom(q = observe_at_least-1,`
#'
#' `       size = round(n_raw*assumed_survival),`
#'
#' `       prob = prop_usedby,`
#'
#' `       lower.tail = FALSE)`
#'
#' or
#'
#' `ppois(q = observe_at_least-1,`
#'
#' `      lambda = n_raw*assumed_survival*prop_usedby,`
#'
#' `      lower.tail = FALSE)`
#'
#' for single inputs.
#'
#' The probability of simultaneously detecting all (or some proportion of) areas
#' used by a given proportion of the marked population is estimated by considering
#' a worst-case scenario in which there are many such use areas, used by the
#' proportion specified by `prop_usedby=`.  For example, if the proportion supplied
#' is 5%, a scenario with twenty (=1/0.05) areas with equal probabilities is
#' considered.  It can be shown that simultaneous detection
#' is itself Binomially distributed, and the probability of simultaneous detection
#' can be given by:
#'
#' \deqn{p(n \geq mq) = \sum_{n \in \{mq,...,m\}}p^n(1-p)^{m-n} {m \choose n}}
#'
#' in which \emph{p} gives the detection probability of a single area,
#' \emph{m} gives the number of areas in a worst-case scenario, and
#' \emph{q} gives the proportion of areas desired to detect.
#'
#' It should be noted that the probability model used for simultaneous detection
#' assumes mutual independence among areas (that is, independent Binomial or
#' Poisson trials, as opposed to Multinomial).  However, this may be a better
#' reflection of reality, as it better allows for additional use areas beyond
#' those considered (that is, those with less proportional usage than the
#' value considered).
#'
#' Multiple values of all arguments may be supplied, in which case the output
#' will be a table rather than a single value.  See examples for details.
#'
#' @param n_raw Trial value of sample size, before accounting for mortality.
#' @param prop_usedby Hypothetical proportion of the population using the area
#' considered.
#' @param assumed_survival Assumed survival (or 1 - data loss proportion).
#' Defaults to `1`.
#' @param observe_at_least Minimum number of marked individuals to consider as
#' detection of an area.  Defaults to `1`, but a larger number may be used as
#' necessary, depending on criterion used to define detection of an aggregation.
#' @param model Assumed underlying probability model.  Allowed values are
#' `"binomial"` and `"poisson"`.  Defaults to `"binomial"`.
#' @param simplify Whether to simplify the output table to only show inputs with
#' multiple values (see examples below).  Defaults to `TRUE`.
#' @param prop_ofareas If simultaneous detection is desired , this gives the
#' proportion of areas to simultaneously detect.  It may be
#' desirable to structure a precision statement in terms of detection of some
#' percentage of areas, see examples below.  If the default `NA` is accepted,
#' simultaneous detection probability will not be calculated.
#' @return Either a single value or table of inputs and calculated probabilities.
#'
#' * Column `$p_singlearea` gives the detection
#' probability of a SINGLE given area.
#'
#' * Element `$p_multipleareas` gives the estimated SIMULTANEOUS detection probability
#' of the proportion of areas specified.
#' @author Matt Tyers
#' @examples
#' ## The probability of detecting a given area used by 5% of the population,
#' ## given a sample size of n=80 and assuming 80% survival, is 96.2%.
#' detection_probability(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8)
#'
#'
#' ## The probability of detecting at least 95% of areas used by 5% of the
#' ## population, given a sample size of n=80 and assuming 80% survival, is at
#' ## least 97%.
#' detection_probability(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8,
#'                       prop_ofareas = 0.95)
#'
#'
#' ## Multiple inputs can be specified, in which case a table will be returned.
#' detection_probability(n_raw = c(80, 100),
#'                       model = c("binomial", "poisson"),
#'                       prop_usedby = 0.05,
#'                       assumed_survival = 0.8,
#'                       prop_ofareas = c(0.9, 1))
#'
#'
#' ## The output table may be expanded with simplify=FALSE.
#' detection_probability(n_raw = c(80, 100),
#'                       model = c("binomial", "poisson"),
#'                       prop_usedby = 0.05,
#'                       assumed_survival = 0.8,
#'                       prop_ofareas = c(0.9, 1),
#'                       simplify = FALSE)
#' @export
detection_probability <- function(n_raw, prop_usedby,
                                  assumed_survival=1, observe_at_least=1,
                                  model="binomial", simplify=TRUE,
                                  prop_ofareas=NA) {
  if(!any(model %in% c("binomial","Binomial","poisson","Poisson"))) {
    stop("model= argument must be \"binomial\" or \"poisson\"")
  }

  # expanding a grid of all possible combinations of input arguments except prop_ofareas
  argmat <- expand.grid(n_raw=n_raw,
                        prop_usedby=prop_usedby,
                        assumed_survival=assumed_survival,
                        observe_at_least=observe_at_least)

  # initializing output objects
  out_binom <- out_pois <- NULL

  # defining possible output objects
  if("binomial" %in% model | "Binomial" %in% model) {
    out_binom <- cbind(argmat,
                       model = "Binomial",
                       p_singlearea = with(argmat, pbinom(q = observe_at_least-1,
                                                          size = round(n_raw*assumed_survival),
                                                          prob = prop_usedby,
                                                          lower.tail = FALSE)))
  }
  if("poisson" %in% model | "Poisson" %in% model) {
    out_pois <- cbind(argmat,
                      model = "Poisson",
                      p_singlearea = with(argmat, ppois(q = observe_at_least-1,
                                                        lambda = n_raw*assumed_survival*prop_usedby,
                                                        lower.tail = FALSE)))
  }
  pmat <- rbind(out_binom, out_pois)   # combining

  # if we also want simultaneous detection probability of multiple areas
  if(all(!is.na(prop_ofareas))) {

    # expanding the prop_ofareas awkwardly
    # maybe it's bad that I'm modifying pmat in place
    pmat <- cbind(do.call(rbind, replicate(length(prop_ofareas), pmat, simplify=FALSE)),
                  prop_ofareas = rep(prop_ofareas, each=nrow(pmat)),
                  p_multipleareas = NA)

    # calculating multiple detection for each row of above!
    for(irow in 1:nrow(pmat)) {  # I'm sure there's a more elegant way to do this, but loops work too
      n_areas <- with(pmat[irow, ], round(1/prop_usedby))   # total number of resulting areas
      n_atleast <- (1:n_areas)[(1:n_areas)/n_areas >= pmat$prop_ofareas[irow]]  # results that could satisfy
      pp <- pmat$p_singlearea[irow]  # just a simplified vbl name

      # actual calculation
      pmat$p_multipleareas[irow] <- sum((pp^n_atleast)*((1-pp)^(n_areas-n_atleast))*choose(n_areas, n_atleast))
    }
  }

  # simplifying
  if(simplify) {
    if(nrow(pmat) > 1) {
      pmat <- pmat[, apply(pmat, 2, function(x) length(unique(x)) > 1)]
    } else {
      pmat <- pmat[, names(pmat) %in% c("p_singlearea", "p_multipleareas")]
    }
  }
  return(pmat)
}
