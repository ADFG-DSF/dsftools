% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_size.R
\name{detection_probability}
\alias{detection_probability}
\title{Detection Probability From a Given Sample Size, Assuming Random Sampling}
\usage{
detection_probability(
  n_raw,
  prop_usedby,
  assumed_survival = 1,
  observe_at_least = 1,
  model = "binomial",
  simplify = TRUE,
  prop_ofareas = NA
)
}
\arguments{
\item{n_raw}{Trial value of sample size, before accounting for mortality.}

\item{prop_usedby}{Hypothetical proportion of the population using the area
considered.}

\item{assumed_survival}{Assumed survival (or 1 - data loss proportion).
Defaults to \code{1}.}

\item{observe_at_least}{Minimum number of marked individuals to consider as
detection of an area.  Defaults to \code{1}, but a larger number may be used as
necessary, depending on criterion used to define detection of an aggregation.}

\item{model}{Assumed underlying probability model.  Allowed values are
\code{"binomial"} and \code{"poisson"}.  Defaults to \code{"binomial"}.}

\item{simplify}{Whether to simplify the output table to only show inputs with
multiple values (see examples below).  Defaults to \code{TRUE}.}

\item{prop_ofareas}{If simultaneous detection is desired , this gives the
proportion of areas to simultaneously detect.  It may be
desirable to structure a precision statement in terms of detection of some
percentage of areas, see examples below.  If the default \code{NA} is accepted,
simultaneous detection probability will not be calculated.}
}
\value{
Either a single value or table of inputs and calculated probabilities.
\itemize{
\item Column \verb{$p_singlearea} gives the detection
probability of a SINGLE given area.
\item Element \verb{$p_multipleareas} gives the estimated SIMULTANEOUS detection probability
of the proportion of areas specified.
}
}
\description{
In the context of a telemetry study, this function estimates the
probability of detecting a single area used by some proportion
of the marked population, given the sample size of instrumented fish, and
optionally the probability of detecting MULTIPLE areas, assuming random sampling from
a worst-case scenario.

The function can use either the Binomial or Poisson probability models, and
the probability of detecting a given (single) use area is calculated from

\verb{pbinom(q = observe_at_least-1,}

\verb{       size = round(n_raw*assumed_survival),}

\verb{       prob = prop_usedby,}

\verb{       lower.tail = FALSE)}

or

\verb{ppois(q = observe_at_least-1,}

\verb{      lambda = n_raw*assumed_survival*prop_usedby,}

\verb{      lower.tail = FALSE)}

for single inputs.

The probability of simultaneously detecting all (or some proportion of) areas
used by a given proportion of the marked population is estimated by considering
a worst-case scenario in which there are many such use areas, used by the
proportion specified by \verb{prop_usedby=}.  For example, if the proportion supplied
is 5\%, a scenario with twenty (=1/0.05) areas with equal probabilities is
considered.  It can be shown that simultaneous detection
is itself Binomially distributed, and the probability of simultaneous detection
can be given by:

\deqn{p(n \geq mq) = \sum_{n \in \{mq,...,m\}}p^n(1-p)^{m-n} {m \choose n}}

in which \emph{p} gives the detection probability of a single area,
\emph{m} gives the number of areas in a worst-case scenario, and
\emph{q} gives the proportion of areas desired to detect.

It should be noted that the probability model used for simultaneous detection
assumes mutual independence among areas (that is, independent Binomial or
Poisson trials, as opposed to Multinomial).  However, this may be a better
reflection of reality, as it better allows for additional use areas beyond
those considered (that is, those with less proportional usage than the
value considered).

Multiple values of all arguments may be supplied, in which case the output
will be a table rather than a single value.  See examples for details.
}
\examples{
## The probability of detecting a given area used by 5\% of the population,
## given a sample size of n=80 and assuming 80\% survival, is 96.2\%.
detection_probability(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8)


## The probability of detecting at least 95\% of areas used by 5\% of the
## population, given a sample size of n=80 and assuming 80\% survival, is at
## least 97\%.
detection_probability(n_raw = 80, prop_usedby = 0.05, assumed_survival = .8,
                      prop_ofareas = 0.95)


## Multiple inputs can be specified, in which case a table will be returned.
detection_probability(n_raw = c(80, 100),
                      model = c("binomial", "poisson"),
                      prop_usedby = 0.05,
                      assumed_survival = 0.8,
                      prop_ofareas = c(0.9, 1))


## The output table may be expanded with simplify=FALSE.
detection_probability(n_raw = c(80, 100),
                      model = c("binomial", "poisson"),
                      prop_usedby = 0.05,
                      assumed_survival = 0.8,
                      prop_ofareas = c(0.9, 1),
                      simplify = FALSE)
}
\author{
Matt Tyers
}
