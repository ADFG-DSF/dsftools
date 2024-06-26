% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asl_sim.R
\name{rp_ASL_table}
\alias{rp_ASL_table}
\title{Relative precision for the Methods Used in \code{ASL_table()}}
\usage{
rp_ASL_table(
  conf_target = 0.9,
  sim = NULL,
  case = NULL,
  nstrata,
  nage,
  nt,
  Nt,
  stratum_weights = NULL,
  se_Nt = NULL,
  mn_length = NULL,
  sd_length = NULL,
  ptz = NULL,
  nNA = 0,
  nsim = 1000,
  plot_pop = TRUE,
  verbose = TRUE,
  print_table = FALSE,
  FPC = c("ifknown", "always", "never")
)
}
\arguments{
\item{conf_target}{The target level of confidence desired.  Defaults to \code{0.9}.}

\item{sim}{An object created by \link{simulate_ASL_table}.  If the default (\code{NULL})
is accepted, the simulation function will be called internally from the other
arguments provided.}

\item{case}{If a pre-determined case is to be used.  Allowed values are
\code{"stratified_witherror_lengthage"}, \code{"stratified_witherror_age"}, \code{"stratified_witherror_length"},
\code{"stratified_lengthage"}, \code{"stratified_age"}, \code{"stratified_length"},
\code{"stratified_Nunknown_lengthage"}, \code{"stratified_Nunknown_age"}, \code{"stratified_Nunknown_length"},
\code{"pooled_witherror_lengthage"}, \code{"pooled_witherror_age"}, \code{"pooled_witherror_length"}
\code{"pooled_lengthage"}, \code{"pooled_age"},  \code{"pooled_length"},
\code{"pooled_Nunknown_lengthage"}, \code{"pooled_Nunknown_length"}, or \code{"pooled_Nunknown_age"}. If the default
(\code{NULL}) is accepted, all simulation parameters below must be supplied by the user.}

\item{nstrata}{Number of sampling strata}

\item{nage}{Number of age categories}

\item{nt}{Sample size for each stratum}

\item{Nt}{Abundance for each stratum.  If abundance is completely unknown, sampling
weights may be used instead in the \verb{sampling_weights=} argument.}

\item{stratum_weights}{Optional vector of sampling weights for each stratum.  Defaults to \code{NULL}.}

\item{se_Nt}{Optional vector of standard errors of abundance for each stratum
Defaults to \code{NULL}, implying abundance is known without error.}

\item{mn_length}{Vector of mean lengths for each age category.  Defaults to
\code{NULL}, implying lengths were not considered.}

\item{sd_length}{Vector of standard deviations for each age category.  Defaults to
\code{NULL}, implying lengths were not considered.}

\item{ptz}{Matrix of probabilities of each age by stratum, with rows
corresponding to strata and columns corresponding to ages.  The probabilities
for each row will be normalized before simulation, therefore summing to one
is not required.  If a pooled (non-stratified) sample is to be taken and
age categories are to be considered, this should be supplied as a vector with
length equal to the number of ages.}

\item{nNA}{Number of NA values to randomly impute, to test robustness to NA.
Defaults to \code{0} in this case.}

\item{nsim}{Number of simulated replicates.  Defaults to \code{1000}, but more is recommended.}

\item{plot_pop}{Whether to make summary plots of the simulated population and
one representative sample, in addition to the plots produced in simulation.
Defaults to \code{TRUE}.}

\item{verbose}{Whether to print the parameters used in simulation to the console,
if one of the \code{case}s is accepted, as well as printing the method used within
\code{ASL_table()}.  Defaults to \code{TRUE}.}

\item{print_table}{Whether to print an example output table from \code{ASL_table()}
as an additional check.  Defaults to \code{FALSE}.}

\item{FPC}{Whether to incorporate the Finite Population Correction factor (FPC) in
variance calculations.  Allowed values are
\itemize{
\item \code{"ifknown"} (the default), which will use the FPC only if abundance is known without error
\item \code{"always"}, which will use the FPC wherever possible (i.e. if there is an estimate of abundance)
\item \code{"never"}, which will always ignore the FPC.
used if the abundance is considered to be known without error.
}}
}
\value{
\code{NULL}
}
\description{
This function may be used with a simulation object created by
\link{simulate_ASL_table}, or will call it internally if none is supplied.

The simulation function draws many samples from a population, given
population and sample characteristics supplied by the user, or else one of
twelve possible cases designed to cover all methods used in \code{ASL_table()}.

This function uses the simulation results to produce a sequence of plots
displaying the relationships between relative precision (how close we want
estimates to be to the true value(s)) and confidence (how often we can expect
estimates to be within the desired relative precision of the truth).

In terms of a typical Operational Plan Objective: "Estimate ... such that all
estimates are within (precision) of the true values (confidence) \% of the time".

Relative precision curves will be created for proportions (phat), abundance (Nhat),
and mean length (mn_length), depending on which are present.  Multiple curves
will be overlayed for each, if multiple categories are present for Age and/or Sex.
}
\examples{
# running rp_ using a pre-defined case
par(mfrow=c(2,2))
rp_ASL_table(case="stratified_witherror_lengthage",
             nsim=500, plot_pop=FALSE)

# running rp_ directly from simulation parameters
par(mfrow=c(2,2))
rp_ASL_table(nstrata = 4,
             nage = 5,

             # sample size for each stratum
             nt = c(100, 100, 100, 100),

             # abundance for each stratum
             Nt = c(10000, 20000, 30000, 40000),

             # (possible) SE for abundance by stratum
             se_Nt = 0.2*c(10000, 20000, 30000, 40000),

             # mean length FOR EACH AGE
             mn_length = c(150, 200, 250, 300, 350),

             # sd length FOR EACH AGE
             sd_length = c(40, 50, 60, 70, 80),

             # matrix of probabilities of each age BY stratum
             ptz = matrix(c(c(1,2,3,4,5),
                            c(1,2,5,5,2),
                            c(2,5,3,2,1),
                            c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5),
             nsim=500, plot_pop=FALSE)


## creating a simulation object first
simresults <- simulate_ASL_table(nstrata = 4,
                                 nage = 5,

                                 # sample size for each stratum
                                 nt = c(100, 100, 100, 100),

                                 # abundance for each stratum
                                 Nt = c(10000, 20000, 30000, 40000),

                                 # (possible) SE for abundance by stratum
                                 se_Nt = 0.2*c(10000, 20000, 30000, 40000),

                                 # mean length FOR EACH AGE
                                 mn_length = c(150, 200, 250, 300, 350),

                                 # sd length FOR EACH AGE
                                 sd_length = c(40, 50, 60, 70, 80),

                                 # matrix of probabilities of each age BY stratum
                                 ptz = matrix(c(c(1,2,3,4,5),
                                                c(1,2,5,5,2),
                                                c(2,5,3,2,1),
                                                c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5),
                                 plot_pop=FALSE)

## running rp_ with the object created
par(mfrow=c(2,2))
rp_ASL_table(sim=simresults)

## running rp_ again
par(mfrow=c(2,2))
rp_ASL_table(sim=simresults, conf_target=0.99)
}
\seealso{
\link{ASL_table}
}
\author{
Matt Tyers
}
