#' Simulation of many replicates of `ASL_table()`
#' @description This function is called internally within
#' \link{verify_ASL_table} and \link{rp_ASL_table}.  However, it may be advantageous
#' to first create simulation results if there will be multiple calls to
#' `verify_ASL_table()` or `rp_ASL_table()`.
#'
#' Eighteen pre-determined cases may be used in the `case=` argument, depending on
#' whether a stratified sampling scheme is used, whether abundance is known,
#' estimated with error, or completely unknown, and whether estimates pertaining to categorical (age)
#' and/or numeric (length) values is to be performed.  In general:
#' * Five age classes are simulated (if age is considered)
#' * Mean and standard deviation of length increase with age (if length and age are considered)
#' * Data are treated as belonging to four temporal strata (if stratified), in which:
#'   - Age class generally decreases with stratum (older fish arrive sooner)
#'   - Abundance and its standard error both increase with stratum (larger numbers of fish arrive later), but
#'   - An equal number of fish are sampled per temporal stratum.
#'
#' It should be noted that a small but non-trival number of entries for age,
#' length, and stratum are imputed as `NA` by default (see `nNA=` argument) in order
#' to test function robustness to missing data.
#'
#' A list will be returned, with matrices of simulation samples and associated
#' true values.  This will consist of some subset of the following, as appropriate:
#'  * Proportions
#'    - `$phat_sim`: A matrix of simulated proportion estimates, with columns corresponding
#'    to categories and rows corresponding to simulation samples.  This is calculated
#'    via the relevant methods used by `ASL_table()`.
#'    - `$phat_true`: A vector of true proportions.  This is calculated from simulation
#'    inputs `Nt` or `ptz` and `stratum_weights` as appropriate.
#'    - `$se_phat_sim`: A matrix of simulated proportion standard errors, with columns corresponding
#'    to categories and rows corresponding to simulation samples.  This is calculated
#'    via the relevant methods used by `ASL_table()`.
#'    - `$se_phat_true`: A vector of standard errors, calculated empirically from
#'    the simulated values `$phat_sim` thereby approximating the sampling distribution.
#'  * Abundances
#'    - `$Nhat_sim`: A matrix of simulated abundance estimates, with columns corresponding
#'    to categories and rows corresponding to simulation samples.  This is calculated
#'    via the relevant methods used by `ASL_table()`.
#'    - `$Nhat_true`: A vector of true abundances, given by user inputs `Nt`.
#'    - `$se_Nhat_sim`: A matrix of simulated abundance standard errors, with columns corresponding
#'    to categories and rows corresponding to simulation samples.  This is calculated
#'    via the relevant methods used by `ASL_table()`.
#'    - `$se_Nhat_true`: A vector of standard errors, calculated empirically from
#'    the simulated values `$Nhat_sim` thereby approximating the sampling distribution.
#'  * Mean Lengths
#'    - `$mn_length_sim`: A matrix of simulated mean length estimates, with columns corresponding
#'    to categories and rows corresponding to simulation samples.  This is calculated
#'    via the relevant methods used by `ASL_table()`.
#'    - `$mn_length_true`: A vector of true mean lengths, given by user inputs `mn_length`.
#'    - `$se_mn_length_sim`: A matrix of simulated mean length standard errors, with columns corresponding
#'    to categories and rows corresponding to simulation samples.  This is calculated
#'    via the relevant methods used by `ASL_table()`.
#'    - `$se_mn_length_true`: A vector of standard errors, calculated empirically from
#'    the simulated values `$mn_length_sim` thereby approximating the sampling distribution.
#' @param case If a pre-determined case is to be used.  Allowed values are
#' `"stratified_witherror_lengthage"`, `"stratified_witherror_age"`, `"stratified_witherror_length"`,
#' `"stratified_lengthage"`, `"stratified_age"`, `"stratified_length"`,
#' `"stratified_Nunknown_lengthage"`, `"stratified_Nunknown_age"`, `"stratified_Nunknown_length"`,
#' `"pooled_witherror_lengthage"`, `"pooled_witherror_age"`, `"pooled_witherror_length"`
#' `"pooled_lengthage"`, `"pooled_age"`,  `"pooled_length"`,
#' `"pooled_Nunknown_lengthage"`, `"pooled_Nunknown_length"`, or `"pooled_Nunknown_age"`. If the default
#' (`NULL`) is accepted, all simulation parameters below must be supplied by the user.
#' @param nstrata Number of sampling strata
#' @param nage Number of age categories
#' @param nt Sample size for each stratum
#' @param Nt Abundance for each stratum.  If abundance is completely unknown, sampling
#' weights may be used instead in the `sampling_weights=` argument.
#' @param stratum_weights Optional vector of sampling weights for each stratum.  Defaults to `NULL`.
#' @param se_Nt Optional vector of standard errors of abundance for each stratum
#' Defaults to `NULL`, implying abundance is known without error.
#' @param mn_length Vector of mean lengths for each age category.  Defaults to
#' `NULL`, implying lengths were not considered.
#' @param sd_length Vector of standard deviations for each age category.  Defaults to
#' `NULL`, implying lengths were not considered.
#' @param ptz Matrix of probabilities of each age by stratum, with rows
#' corresponding to strata and columns corresponding to ages.  The probabilities
#' for each row will be normalized before simulation, therefore summing to one
#' is not required.  If a pooled (non-stratified) sample is to be taken and
#' age categories are to be considered, this should be supplied as a vector with
#' length equal to the number of ages.
#' @param nNA Number of NA values to randomly impute, to test robustness to NA.  Defaults to `0`.
#' @param nsim Number of simulated replicates.  Defaults to `1000`, but more is recommended.
#' @param plot_pop Whether to make summary plots of the simulated population and
#' one representative sample, in addition to the plots produced in simulation.
#' Defaults to `TRUE`.
#' @param verbose Whether to print the parameters used in simulation to the console,
#' if one of the `case`s is accepted, as well as printing the method used within
#' `ASL_table()`.  Defaults to `TRUE`.
#' @param print_table Whether to print an example output table from `ASL_table()`
#' as an additional check.  Defaults to `FALSE`.
#' @return An object of class `"ASL_sim"`, which is a list of matrices of
#' simulation samples, and vectors of the associated
#' true values.  Details are described above.
#' @author Matt Tyers
#' @seealso [rp_ASL_table], [verify_ASL_table], [ASL_table]
#' @examples
#' ## creating a simulation object using case=
#' simresults <- simulate_ASL_table(case="stratified_witherror_lengthage",
#'                                  nsim=500, plot_pop=FALSE)
#'
#' ## creating a simulation object using simulation parameters
#' simresults <- simulate_ASL_table(nstrata = 4,
#'                                  nage = 5,
#'
#'                                  # sample size for each stratum
#'                                  nt = c(100, 100, 100, 100),
#'
#'                                  # abundance for each stratum
#'                                  Nt = c(10000, 20000, 30000, 40000),
#'
#'                                  # (possible) SE for abundance by stratum
#'                                  se_Nt = 0.2*c(10000, 20000, 30000, 40000),
#'
#'                                  # mean length FOR EACH AGE
#'                                  mn_length = c(150, 200, 250, 300, 350),
#'
#'                                  # sd length FOR EACH AGE
#'                                  sd_length = c(40, 50, 60, 70, 80),
#'
#'                                  # matrix of probabilities of each age BY stratum
#'                                  ptz = matrix(c(c(1,2,3,4,5),
#'                                                 c(1,2,5,5,2),
#'                                                 c(2,5,3,2,1),
#'                                                 c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5),
#'                                  nsim=500, plot_pop=FALSE)
#'
#' ## running rp_ with the object created
#' par(mfrow=c(2,2))
#' rp_ASL_table(sim=simresults)
#'
#' ## running rp_ again
#' par(mfrow=c(2,2))
#' rp_ASL_table(sim=simresults, conf_target=0.99)
#'
#' ## running verify_
#' par(mfrow=c(3,2))
#' verify_ASL_table(sim=simresults)
#' @export
simulate_ASL_table <- function(case=NULL,   # should this default to NULL?
                               nstrata,
                               nage,
                               nt,   # c(100, 100, 100, 100), # sample size for each stratum
                               Nt,   # c(10000, 20000, 30000, 40000),  # abundance for each stratum
                               stratum_weights = NULL,
                               se_Nt = NULL,   # 0.2*Nt, # (possible) SE for abundance by stratum,
                               mn_length = NULL,   # c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
                               sd_length = NULL,   # c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
                               ptz = NULL,   # matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                               # c(1,2,5,5,2),
                               # c(2,5,3,2,1),
                               # c(5,4,3,1,1)),
                               # byrow=TRUE,
                               # nrow=4,
                               # ncol=5),
                               nNA=0,
                               nsim=1000,   # 1000,    # number of simulated replicates
                               plot_pop=TRUE,   # whether to make summary plots of pop & one sample
                               verbose=TRUE, # whether to output cases in sim function
                               print_table=FALSE) {

  ### pre-populate cases
  if(!is.null(case)) {
    if(case=="stratified_witherror_lengthage") {
      nstrata <- 4
      nage <- 5
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
      se_Nt <- 0.2*c(10000, 20000, 30000, 40000) # (possible) SE for abundance by stratum,
      mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
      sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
      ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                      c(1,2,5,5,2),
                      c(2,5,3,2,1),
                      c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
      if(verbose) {
        cat("
nstrata = 4,
nage = 5,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
se_Nt = 0.2*c(10000, 20000, 30000, 40000), # (possible) SE for abundance by stratum
mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
               c(1,2,5,5,2),
               c(2,5,3,2,1),
               c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)

")
      }
    }
    if(case=="stratified_witherror_age") {
      nstrata <- 4
      nage <- 5
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
      se_Nt <- 0.2*c(10000, 20000, 30000, 40000) # (possible) SE for abundance by stratum,
      mn_length <- NULL # mean length FOR EACH AGE
      sd_length <- NULL # sd length FOR EACH AGE
      ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                      c(1,2,5,5,2),
                      c(2,5,3,2,1),
                      c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
      # mn_length=NULL, sd_length=NULL
      if(verbose) {
        cat("
nstrata = 4,
nage = 5,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
se_Nt = 0.2*c(10000, 20000, 30000, 40000), # (possible) SE for abundance by stratum
mn_length = NULL, # mean length FOR EACH AGE
sd_length = NULL, # sd length FOR EACH AGE
ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                c(1,2,5,5,2),
                c(2,5,3,2,1),
                c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
")
      }
    }
    if(case=="stratified_witherror_length") {
      nstrata <- 4
      nage <- 1
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
      se_Nt <- 0.2*c(10000, 20000, 30000, 40000) # (possible) SE for abundance by stratum,
      mn_length <- 300 # mean length FOR EACH AGE
      sd_length <- 70 # sd length FOR EACH AGE
      ptz <- NULL
      # mn_length=300, sd_length=70, ptz=NULL
      if(verbose) {
        cat("
nstrata = 4,
nage = 1,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
se_Nt = 0.2*c(10000, 20000, 30000, 40000), # (possible) SE for abundance by stratum
mn_length = 300, # mean length FOR EACH AGE
sd_length = 70, # sd length FOR EACH AGE
ptz = NULL
")
      }
    }
    if(case=="stratified_lengthage") {
      nstrata <- 4
      nage <- 5
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
      sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
      ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                      c(1,2,5,5,2),
                      c(2,5,3,2,1),
                      c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
      if(verbose) {
        cat("
nstrata = 4,
nage = 5,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                c(1,2,5,5,2),
                c(2,5,3,2,1),
                c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
")
      }
      # se_Nt = NULL
    }
    if(case=="stratified_age") {
      nstrata <- 4
      nage <- 5
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- NULL # mean length FOR EACH AGE
      sd_length <- NULL # sd length FOR EACH AGE
      ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                      c(1,2,5,5,2),
                      c(2,5,3,2,1),
                      c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
      # se_Nt = NULL, mn_length=NULL, sd_length=NULL
      if(verbose) {
        cat("
nstrata = 4,
nage = 5,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = NULL, # mean length FOR EACH AGE
sd_length = NULL, # sd length FOR EACH AGE
ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                c(1,2,5,5,2),
                c(2,5,3,2,1),
                c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
")
      }
    }
    if(case=="stratified_length") {
      nstrata <- 4
      nage <- 1
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- c(10000, 20000, 30000, 40000)  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- 300 # mean length FOR EACH AGE
      sd_length <- 70 # sd length FOR EACH AGE
      ptz <- NULL
      # se_Nt = NULL, mn_length=300, sd_length=70, ptz=NULL
      if(verbose) {
        cat("
nstrata = 4,
nage = 1,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = c(10000, 20000, 30000, 40000),  # abundance for each stratum
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = 300, # mean length FOR EACH AGE
sd_length = 70, # sd length FOR EACH AGE
ptz = NULL
")
      }
    }
    if(case=="pooled_witherror_lengthage") {
      nstrata <- 1
      nage <- 5
      nt <- 400 # sample size for each stratum
      Nt <- 10000  # abundance for each stratum
      se_Nt <- 0.2*10000 # (possible) SE for abundance by stratum,
      mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
      sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
      ptz <- 1:5
      # ptz=1:5, Nt=10000, nt=100
      if(verbose) {
        cat("
nstrata = 1,
nage = 5,
nt = 400, # sample size for each stratum
Nt = 10000,  # abundance for each stratum
se_Nt = 0.2*10000, # (possible) SE for abundance by stratum
mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
ptz = 1:5
")
      }
    }
    if(case=="pooled_witherror_age") {
      nstrata <- 1
      nage <- 5
      nt <- 400 # sample size for each stratum
      Nt <- 10000  # abundance for each stratum
      se_Nt <- 0.2*10000 # (possible) SE for abundance by stratum,
      mn_length <- NULL # mean length FOR EACH AGE
      sd_length <- NULL # sd length FOR EACH AGE
      ptz <- 1:5
      # mn_length=NULL, sd_length=NULL, ptz=1:5, Nt=10000, nt=100
      if(verbose) {
        cat("
nstrata = 1,
nage = 5,
nt = 400, # sample size for each stratum
Nt = 10000,  # abundance for each stratum
se_Nt = 0.2*10000, # (possible) SE for abundance by stratum
mn_length = NULL, # mean length FOR EACH AGE
sd_length = NULL, # sd length FOR EACH AGE
ptz = 1:5
")
      }
    }
    if(case=="pooled_witherror_length") {
      nstrata <- 1
      nage <- 1
      nt <- 400 # sample size for each stratum
      Nt <- 10000  # abundance for each stratum
      se_Nt <- 0.2*10000 # (possible) SE for abundance by stratum,
      mn_length <- 300 # mean length FOR EACH AGE
      sd_length <- 70 # sd length FOR EACH AGE
      ptz <- NULL
      # mn_length=300, sd_length=70, ptz=NULL, Nt=10000, nt=100
      if(verbose) {
        cat("
nstrata = 1,
nage = 1,
nt = 400, # sample size for each stratum
Nt = 10000,  # abundance for each stratum
se_Nt = 0.2*10000, # (possible) SE for abundance by stratum
mn_length = 300, # mean length FOR EACH AGE
sd_length = 70, # sd length FOR EACH AGE
ptz = NULL
")
      }
    }
    if(case=="pooled_lengthage") {
      nstrata <- 1
      nage <- 5
      nt <- 400 # sample size for each stratum
      Nt <- 10000  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
      sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
      ptz <- 1:5
      # se_Nt = NULL, ptz=1:5, Nt=10000, nt=100
      if(verbose) {
        cat("
nstrata = 1,
nage = 5,
nt = 400, # sample size for each stratum
Nt = 10000,  # abundance for each stratum
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
ptz = 1:5
")
      }
    }
    if(case=="pooled_age") {
      nstrata <- 1
      nage <- 5
      nt <- 400 # sample size for each stratum
      Nt <- 10000  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- NULL # mean length FOR EACH AGE
      sd_length <- NULL # sd length FOR EACH AGE
      ptz <- 1:5
      # se_Nt = NULL, mn_length=NULL, sd_length=NULL, ptz=1:5, Nt=10000, nt=100
      if(verbose) {
        cat("
nstrata = 1,
nage = 5,
nt = 400, # sample size for each stratum
Nt = 10000,  # abundance for each stratum
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = NULL, # mean length FOR EACH AGE
sd_length = NULL, # sd length FOR EACH AGE
ptz = 1:5
")
      }
    }
    if(case=="pooled_length") {
      nstrata <- 1
      nage <- 1
      nt <- 400 # sample size for each stratum
      Nt <- 10000 # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- 300 # mean length FOR EACH AGE
      sd_length <- 70 # sd length FOR EACH AGE
      ptz <- NULL
      # se_Nt = NULL, mn_length=300, sd_length=70, ptz=NULL, Nt=10000, nt=100
      if(verbose) {
        cat("
nstrata = 1,
nage = 1,
nt = 400, # sample size for each stratum
Nt = 10000, # abundance for each stratum
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = 300, # mean length FOR EACH AGE
sd_length = 70, # sd length FOR EACH AGE
ptz = NULL
")
      }
    }
    if(case=="stratified_Nunknown_lengthage") {
      nstrata <- 4
      nage <- 5
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- NULL
      stratum_weights <- 1:4  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
      sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
      ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                      c(1,2,5,5,2),
                      c(2,5,3,2,1),
                      c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
      if(verbose) {
        cat("
nstrata = 4,
nage = 5,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = NULL,  # abundance for each stratum
stratum_weights = 1:4, # stratum weights instead of abundance
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
               c(1,2,5,5,2),
               c(2,5,3,2,1),
               c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)

")
      }
    }
    if(case=="stratified_Nunknown_length") {
      nstrata <- 4
      nage <- 1
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- NULL
      stratum_weights <- 1:4  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- 300 # mean length FOR EACH AGE
      sd_length <- 70 # sd length FOR EACH AGE
      ptz <- NULL  # matrix of probabilities of each age BY stratum

      if(verbose) {
        cat("
nstrata = 4,
nage = 1,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = NULL,  # abundance for each stratum
stratum_weights = 1:4, # stratum weights instead of abundance
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = 300, # mean length FOR EACH AGE
sd_length = 70, # sd length FOR EACH AGE
ptz = NULL  # matrix of probabilities of each age BY stratum

")
      }
    }
    if(case=="stratified_Nunknown_age") {
      nstrata <- 4
      nage <- 5
      nt <- c(100, 100, 100, 100) # sample size for each stratum
      Nt <- NULL
      stratum_weights <- 1:4  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- NULL # mean length FOR EACH AGE
      sd_length <- NULL # sd length FOR EACH AGE
      ptz <- matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                      c(1,2,5,5,2),
                      c(2,5,3,2,1),
                      c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)
      if(verbose) {
        cat("
nstrata = 4,
nage = 5,
nt = c(100, 100, 100, 100), # sample size for each stratum
Nt = NULL,  # abundance for each stratum
stratum_weights = 1:4, # stratum weights instead of abundance
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = NULL, # mean length FOR EACH AGE
sd_length = NULL, # sd length FOR EACH AGE
ptz = matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
               c(1,2,5,5,2),
               c(2,5,3,2,1),
               c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5)

")
      }
    }
    if(case=="pooled_Nunknown_lengthage") {
      nstrata <- 1
      nage <- 5
      nt <- 400 # sample size for each stratum
      Nt <- NULL
      stratum_weights <- 1  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- c(150, 200, 250, 300, 350) # mean length FOR EACH AGE
      sd_length <- c(40, 50, 60, 70, 80) # sd length FOR EACH AGE
      ptz <- 1:5
      if(verbose) {
        cat("
nstrata = 1,
nage = 5,
nt = 400, # sample size for each stratum
Nt = NULL,  # abundance for each stratum
stratum_weights = 1, # stratum weights instead of abundance
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
sd_length = c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
ptz = 1:5

")
      }
    }
    if(case=="pooled_Nunknown_length") {
      nstrata <- 1
      nage <- 1
      nt <- 400 # sample size for each stratum
      Nt <- NULL
      stratum_weights <- 1  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- 300 # mean length FOR EACH AGE
      sd_length <- 70 # sd length FOR EACH AGE
      ptz <- NULL  # matrix of probabilities of each age BY stratum

      if(verbose) {
        cat("
nstrata = 1,
nage = 1,
nt = 400, # sample size for each stratum
Nt = NULL,  # abundance for each stratum
stratum_weights = 1, # stratum weights instead of abundance
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = 300, # mean length FOR EACH AGE
sd_length = 70, # sd length FOR EACH AGE
ptz = NULL  # matrix of probabilities of each age BY stratum

")
      }
    }
    if(case=="pooled_Nunknown_age") {
      nstrata <- 1
      nage <- 5
      nt <- 400 # sample size for each stratum
      Nt <- NULL
      stratum_weights <- 1  # abundance for each stratum
      se_Nt <- NULL # (possible) SE for abundance by stratum,
      mn_length <- NULL # mean length FOR EACH AGE
      sd_length <- NULL # sd length FOR EACH AGE
      ptz <- 1:5
      if(verbose) {
        cat("
nstrata = 1,
nage = 5,
nt = 400, # sample size for each stratum
Nt = NULL,  # abundance for each stratum
stratum_weights = 1, # stratum weights instead of abundance
se_Nt = NULL, # (possible) SE for abundance by stratum
mn_length = NULL, # mean length FOR EACH AGE
sd_length = NULL, # sd length FOR EACH AGE
ptz = 1:5

")
      }
    }
  }

  # pre-fixing ptz
  if(is.null(dim(ptz)) & !is.null(ptz)) {
    ptz <- t(as.matrix(ptz))
  }

  # check inputs: length(nt) vs length(Nt) vs. length(se_Nt) vs nrow(ptz)
  if(length(nt)!=nstrata) stop("length(nt) should equal nstrata.")
  if(length(Nt)!=nstrata & length(stratum_weights)!=nstrata) stop("length(Nt) or length(stratum_weights) should equal nstrata.")
  if(!is.null(se_Nt) & (length(se_Nt)!=nstrata)) stop("length(se_Nt) should equal nstrata.")
  if(!is.null(ptz)) {
    if(nrow(ptz)!=nstrata) stop("nrow(ptz) should equal nstrata.")
  }

  # check inputs: length(mn_length) vs length(sd_length) vs ncol(ptz)
  if(!is.null(mn_length) & (length(mn_length)!=nage)) stop("length(mn_length) should equal nage.")
  if(!is.null(sd_length) & (length(sd_length)!=nage)) stop("length(sd_length) should equal nage.")
  if(length(sd_length) != length(mn_length)) stop("mn_length and sd_length should have corresponding length.")
  if(!is.null(ptz)) {
    if(nage>1 & ncol(ptz)!=nage) stop("ncol(ptz) or length(ptz) must equal nage.")
  }

  #### no age means mn_length and sd_length have length 1, and ptz will be NULL
  #### no length means mn_length and sd_length will be NULL
  #### no strata means nt and Nt will be length 1 and ptz will be a vector

  # nstrata <- length(nt)


  ## constructing vectors at the POPULATION level
  if(!is.null(Nt)) {
    # stratum
    t <- rep(1,Nt[1])
    if(length(Nt) > 1) {
      for(i in 2:length(Nt)) t <- c(t, rep(i, Nt[i]))
    }

    # age
    # nage <- ifelse(!is.null(mn_length), length(mn_length), ncol(ptz))
    if(!is.null(ptz)) {
      age_sim <- NA*t
      for(i_t in 1:nstrata) {
        age_sim[t==i_t] <- sample(x=1:nage, size=sum(t==i_t), prob=ptz[i_t,],replace=TRUE)
      }
      # age_sim[t==1] <- sample(x=1:5, size=sum(t==1), prob=c(1,2,3,4,5), replace=T)
    } else {
      age_sim <- rep(1, length(t))
    }

    # length
    if(!is.null(mn_length)) {
      length_sim <- NA*t
      for(i_z in 1:nage) {
        length_sim[age_sim==i_z] <- rnorm(sum(age_sim==i_z), mn_length[i_z], sd_length[i_z])
      }
    } else {
      length_sim <- NULL
    }



    # # storing graphics state to save
    # parmfrow <- par("mfrow")
    # on.exit(par(mfrow=parmfrow))  # making sure to re-set graphics state

    # # plotting simulated population
    # mosaicplot(table(t,age_sim), xlab="Stratum", ylab="Age")
    # boxplot(length_sim~age_sim, xlab="Age", ylab="Length")
    # boxplot(length_sim~t, xlab="Stratum", ylab="Length")



    # simulate once to set up output array (yes this is a hack)
    thesample <- sample(seq_along(t)[t==1], size=nt[1])
    if(length(nt) > 1) {
      for(i_t in 2:nstrata) {
        thesample <- c(thesample, sample(seq_along(t)[t==i_t], size=nt[i_t]))
      }
    }
    # table(t[thesample])  # making sure it worked
  } else {
    # t is stratum at pop level
    # age from stratum
    age_sim <- NULL
    t <- NULL
    for(i_nt in 1:length(nt)) {
      if(!is.null(ptz)) {
        age_sim <- c(age_sim, sample(1:nage, size=nt[i_nt], prob=ptz[i_nt, ], replace=TRUE))
      } else {
        age_sim <- c(age_sim, rep(1, nt[i_nt]))
      }
      t <- c(t, rep(i_nt, nt[i_nt]))
    }

    # length from age
    if(!is.null(mn_length)) {
      length_sim <- rnorm(length(t), mean=mn_length[age_sim], sd=sd_length[age_sim])
    } else {
      length_sim <- NULL
    }

    # dummy vector to maintain compatibility with existing code
    thesample <- seq_along(t)
  }

  # plotting population & one sample (optionally)
  if(plot_pop) {
    if(prod(dim(table(t,age_sim))) > 1) {  #
      if(!is.null(Nt)) mosaicplot(table(t,age_sim), xlab="Stratum", ylab="Age", main="Population", col=grey.colors(nage, rev=TRUE))
      mosaicplot(table(t[thesample],age_sim[thesample]), xlab="Stratum", ylab="Age", main="Sample", col=grey.colors(nage, rev=TRUE))
    }
    if(!is.null(length_sim)) {
      if(!is.null(Nt)) boxplot(length_sim~age_sim, xlab="Age", ylab="Length", main="Population")
      boxplot(length_sim[thesample]~age_sim[thesample], xlab="Age", ylab="Length", main="Sample")
      if(!is.null(Nt)) boxplot(length_sim~t, xlab="Stratum", ylab="Length", main="Population")
      boxplot(length_sim[thesample]~t[thesample], xlab="Stratum", ylab="Length", main="Sample")
    }
  }

  if(is.null(se_Nt)) {
    Nhat_sim <- Nt
  } else {
    Nhat_sim <- rnorm(length(Nt), mean=Nt, sd=se_Nt) # could move this into function args
  }

  ### LEAVE IT OPEN FOR age_sim=NULL and length_sim=NULL
  if(nage == 1) age_sim <- NULL
  # if(nstrata == 1) t <- NULL

  thestratum <- as.integer(t[thesample])
  if(all(thestratum==1)) thestratum <- NULL

  theage <- age_sim[thesample]
  thelength <- length_sim[thesample]


  ## imputing some NA values to make sure the function is robust to NA!!
  if(!is.null(theage)) theage[sample(seq_along(theage), nNA)] <- NA
  if(!is.null(thelength)) thelength[sample(seq_along(thelength), nNA)] <- NA
  if(!is.null(thestratum)) thestratum[sample(seq_along(thestratum), nNA)] <- NA


  thetable <- ASL_table(age=theage,
                        length=thelength,
                        stratum=thestratum,
                        Nhat=Nhat_sim,
                        stratum_weights = stratum_weights,
                        se_Nhat=se_Nt,
                        verbose=verbose)

  # thetable <- ASL_table(age=age_sim[thesample],
  #                       length=length_sim[thesample],
  #                       stratum=thestratum,
  #                       Nhat=Nhat_sim,
  #                       se_Nhat=se_Nt,
  #                       verbose=verbose)  # find a way to add stratum_weights???
  if(print_table) print(thetable)

  # initiate all stuff to store
  # results <- array(dim=c(nrow(thetable), ncol(thetable), nsim))
  results <- array(dim=c(nsim, nrow(thetable), ncol(thetable)))

  ############### then loop nsim times!!! ###############
  for(i_sim in 1:nsim) {
    if(!is.null(Nt)) {
      # first take a sample of fish
      thesample <- sample(seq_along(t)[t==1], size=nt[1])
      if(length(nt) > 1) {
        for(i_t in 2:nstrata) {
          thesample <- c(thesample, sample(seq_along(t)[t==i_t], size=nt[i_t]))
        }
      }
      # table(t[thesample])  # making sure it worked
    } else {
      # t is stratum at pop level
      # age from stratum
      age_sim <- NULL
      t <- NULL
      for(i_nt in 1:length(nt)) {
        if(!is.null(ptz)) {
          age_sim <- c(age_sim, sample(1:nage, size=nt[i_nt], prob=ptz[i_nt, ], replace=TRUE))
        } else {
          age_sim <- c(age_sim, rep(1, nt[i_nt]))
        }
        t <- c(t, rep(i_nt, nt[i_nt]))
      }

      # length from age
      if(!is.null(mn_length)) {
        length_sim <- rnorm(length(t), mean=mn_length[age_sim], sd=sd_length[age_sim])
      } else {
        length_sim <- NULL
      }

      # dummy vector to maintain compatibility with existing code
      thesample <- seq_along(t)
      if(nage==1) age_sim <- NULL
    }

    if(is.null(se_Nt)) {
      Nhat_sim <- Nt
    } else {
      Nhat_sim <- rnorm(length(Nt), mean=Nt, sd=se_Nt) # could move this into function args
    }

    thestratum <- as.integer(t[thesample])    #### this is new
    if(all(thestratum==1)) thestratum <- NULL    #### this is new

    theage <- age_sim[thesample]
    thelength <- length_sim[thesample]


    ## imputing some NA values to make sure the function is robust to NA!!
    if(!is.null(theage)) theage[sample(seq_along(theage), nNA)] <- NA
    if(!is.null(thelength)) thelength[sample(seq_along(thelength), nNA)] <- NA
    if(!is.null(thestratum)) thestratum[sample(seq_along(thestratum), nNA)] <- NA


    thetable <- ASL_table(age=theage,
                          length=thelength,
                          stratum=thestratum,
                          Nhat=Nhat_sim,
                          stratum_weights=stratum_weights,
                          se_Nhat=se_Nt)  # find a way to add stratum_weights???
    # results[,,i_sim] <- as.matrix(thetable)
    results[i_sim,,] <- as.matrix(thetable)
  }
  # dimnames(results)[1:2] <- dimnames(thetable)
  dimnames(results)[2:3] <- dimnames(thetable)

  # return(results)

  ### bundling simulation results
  outlist <- list()
  if("phat" %in% dimnames(results)[[3]]) {
    outlist$phat_sim <- results[,,dimnames(results)[[3]]=="phat"]
    if(!is.null(Nt)) {
      outlist$phat_true <- as.numeric(table(age_sim)/sum(Nt))
    } else {
      outlist$phat_true <- colSums(ptz*stratum_weights)/sum(ptz*stratum_weights)
    }


    outlist$se_phat_sim <- results[,,dimnames(results)[[3]]=="se_phat"]
    outlist$se_phat_true <- apply(outlist$phat_sim, 2, sd, na.rm=TRUE)

  }

  if("Nhat" %in% dimnames(results)[[3]]) {
    outlist$Nhat_sim <- results[,,dimnames(results)[[3]]=="Nhat"]
    outlist$Nhat_true <- as.numeric(table(age_sim))

    outlist$se_Nhat_sim <- results[,,dimnames(results)[[3]]=="se_Nhat"]
    outlist$se_Nhat_true <- apply(outlist$Nhat_sim, 2, sd, na.rm=TRUE)

  }

  if("mn_length" %in% dimnames(results)[[3]]) {
    outlist$mn_length_sim <- as.matrix(results[,,dimnames(results)[[3]]=="mn_length"])
    outlist$mn_length_true <- mn_length

    outlist$se_mn_length_sim <- results[,,dimnames(results)[[3]]=="se_length"]
    outlist$se_mn_length_true <- apply(outlist$mn_length_sim, 2, sd, na.rm=TRUE)

  }

  class(outlist) <- "ASL_sim"
  return(outlist)
}





#' Empirical Verification of the Methods Used in `ASL_table()`
#' @description This function may be used with a simulation object created by
#' \link{simulate_ASL_table}, or will call it internally if none is supplied.
#'
#' The simulation function draws many samples from a population, given
#' population and sample characteristics supplied by the user, or else one of
#' twelve possible cases designed to cover all methods used in `ASL_table()`.
#'
#' The sampling distributions are approximated by simulation for all estimators of age
#' proportions, abundances, and mean lengths, and their respective standard
#' errors. The sampling distributions are then graphically compared to true
#' values (in the case of simulated point estimates) or empirical standard deviations of
#' simulated point estimates (in the case of simulated standard errors.)
#'
#' It should be noted that all simulated estimation is performed using the
#' current version of the `ASL_table()` function; therefore any bias or errors
#' that are identified should indicate a bug or incorrect derivation!
#'
#' Eighteen pre-determined cases may be used in the `case=` argument, depending on
#' whether a stratified sampling scheme is used, whether abundance is known,
#' estimated with error, or completely unknown, and whether estimates pertaining to categorical (age)
#' and/or numeric (length) values is to be performed.  In general:
#' * Five age classes are simulated (if age is considered)
#' * Mean and standard deviation of length increase with age (if length and age are considered)
#' * Data are treated as belonging to four temporal strata (if stratified), in which:
#'   - Age class generally decreases with stratum (older fish arrive sooner)
#'   - Abundance and its standard error both increase with stratum (larger numbers of fish arrive later), but
#'   - An equal number of fish are sampled per temporal stratum.
#'
#' It should be noted that a small but non-trival number of entries for age,
#' length, and stratum are imputed as `NA` by default (see `nNA=` argument) in order
#' to test function robustness to missing data.
#' @param sim An object created by \link{simulate_ASL_table}.  If the default (`NULL`)
#' is accepted, the simulation function will be called internally from the other
#' arguments provided.
#' @param case If a pre-determined case is to be used.  Allowed values are
#' `"stratified_witherror_lengthage"`, `"stratified_witherror_age"`, `"stratified_witherror_length"`,
#' `"stratified_lengthage"`, `"stratified_age"`, `"stratified_length"`,
#' `"stratified_Nunknown_lengthage"`, `"stratified_Nunknown_age"`, `"stratified_Nunknown_length"`,
#' `"pooled_witherror_lengthage"`, `"pooled_witherror_age"`, `"pooled_witherror_length"`
#' `"pooled_lengthage"`, `"pooled_age"`,  `"pooled_length"`,
#' `"pooled_Nunknown_lengthage"`, `"pooled_Nunknown_length"`, or `"pooled_Nunknown_age"`. If the default
#' (`NULL`) is accepted, all simulation parameters below must be supplied by the user.
#' @param nstrata Number of sampling strata
#' @param nage Number of age categories
#' @param nt Sample size for each stratum
#' @param Nt Abundance for each stratum.  If abundance is completely unknown, sampling
#' weights may be used instead in the `sampling_weights=` argument.
#' @param stratum_weights Optional vector of sampling weights for each stratum.  Defaults to `NULL`.
#' @param se_Nt Optional vector of standard errors of abundance for each stratum
#' Defaults to `NULL`, implying abundance is known without error.
#' @param mn_length Vector of mean lengths for each age category.  Defaults to
#' `NULL`, implying lengths were not considered.
#' @param sd_length Vector of standard deviations for each age category.  Defaults to
#' `NULL`, implying lengths were not considered.
#' @param ptz Matrix of probabilities of each age by stratum, with rows
#' corresponding to strata and columns corresponding to ages.  The probabilities
#' for each row will be normalized before simulation, therefore summing to one
#' is not required.  If a pooled (non-stratified) sample is to be taken and
#' age categories are to be considered, this should be supplied as a vector with
#' length equal to the number of ages.
#' @param nNA Number of NA values to randomly impute, to test robustness to NA.  Defaults to `10`.
#' @param nsim Number of simulated replicates.  Defaults to `1000`, but more is recommended.
#' @param plot_pop Whether to make summary plots of the simulated population and
#' one representative sample, in addition to the plots produced in simulation.
#' Defaults to `TRUE`.
#' @param verbose Whether to print the parameters used in simulation to the console,
#' if one of the `case`s is accepted, as well as printing the method used within
#' `ASL_table()`.  Defaults to `TRUE`.
#' @param print_table Whether to print an example output table from `ASL_table()`
#' as an additional check.  Defaults to `FALSE`.
#' @return `NULL`
#' @author Matt Tyers
#' @seealso [ASL_table]
#' @examples
#' # running using a pre-defined case
#' par(mfrow=c(2,2))
#' verify_ASL_table(case="stratified_witherror_lengthage", nsim=500)
#'
#' # running using simulation parameters directly
#' par(mfrow=c(2,2))
#' verify_ASL_table(nstrata = 4,
#'                  nage = 5,
#'
#'                  # sample size for each stratum
#'                  nt = c(100, 100, 100, 100),
#'
#'                  # abundance for each stratum
#'                  Nt = c(10000, 20000, 30000, 40000),
#'
#'                  # (possible) SE for abundance by stratum
#'                  se_Nt = 0.2*c(10000, 20000, 30000, 40000),
#'
#'                  # mean length FOR EACH AGE
#'                  mn_length = c(150, 200, 250, 300, 350),
#'
#'                  # sd length FOR EACH AGE
#'                  sd_length = c(40, 50, 60, 70, 80),
#'
#'                  # matrix of probabilities of each age BY stratum
#'                  ptz = matrix(c(c(1,2,3,4,5),
#'                                 c(1,2,5,5,2),
#'                                 c(2,5,3,2,1),
#'                                 c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5),
#'                  nsim=1000)
#'
#' # simulating first, then running
#' par(mfrow=c(2,2))
#' simresults <- simulate_ASL_table(case="stratified_witherror_lengthage", nsim=500)
#' verify_ASL_table(simresults)
#'
#'
#' # all pre-defined cases
#' \dontrun{
#' nsim <- 5000
#' cases <- c("stratified_witherror_lengthage", "stratified_witherror_age",
#' "stratified_witherror_length", "stratified_lengthage", "stratified_age",
#' "stratified_length", "stratified_Nunknown_lengthage", "stratified_Nunknown_age",
#' "stratified_Nunknown_length", "pooled_witherror_lengthage", "pooled_witherror_age",
#' "pooled_witherror_length", "pooled_lengthage", "pooled_age", "pooled_length",
#' "pooled_Nunknown_lengthage", "pooled_Nunknown_age", "pooled_Nunknown_length")
#' for(case_i in cases) {
#'   par(mfrow=c(3,2))
#'   verify_ASL_table(case=case_i, nsim=nsim)
#' }
#' }
#' @export
verify_ASL_table <- function(sim=NULL,
                              case=NULL,   # should this default to NULL?
                              nstrata,
                              nage,
                              nt,   # c(100, 100, 100, 100), # sample size for each stratum
                              Nt,   # c(10000, 20000, 30000, 40000),  # abundance for each stratum
                              stratum_weights = NULL,
                              se_Nt = NULL,   # 0.2*Nt, # (possible) SE for abundance by stratum,
                              mn_length = NULL,   # c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
                              sd_length = NULL,   # c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
                              ptz = NULL,   # matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                              # c(1,2,5,5,2),
                              # c(2,5,3,2,1),
                              # c(5,4,3,1,1)),
                              # byrow=TRUE,
                              # nrow=4,
                              # ncol=5),
                              nNA=10,
                              nsim=1000,   # 1000,    # number of simulated replicates
                              plot_pop=TRUE,   # whether to make summary plots of pop & one sample
                              verbose=TRUE, # whether to output cases in sim function
                              print_table=FALSE) {

  ## if sim= is supplied, make sure it is actually a simulation object!
  if(!is.null(sim)) {
    if(!inherits(sim, "ASL_sim")) {
      stop("Argument sim= must be an object returned from simulate_ASL_table()")
    }
  }

  ## if a pre-run simulation is not supplied, then run the simulation
  if(is.null(sim)) {

    # validate (class match?) that sim is useable

    sim <- simulate_ASL_table(case=case,
                              nstrata=nstrata,
                              nage=nage,
                              nt=nt,
                              Nt=Nt,
                              stratum_weights = stratum_weights,
                              se_Nt=se_Nt,
                              mn_length=mn_length,
                              sd_length=sd_length,
                              ptz=ptz,
                              nNA=nNA,
                              nsim=nsim,
                              plot_pop=plot_pop,
                              verbose=verbose,
                              print_table=print_table)
  }

  ### I'm sure there's a more efficient way to use all arguments!!

  ### plotting simulation results, overlayed with true values
  with(sim, {  # this is a hack, curious if it works
    if(!is.null(sim$phat_sim)) {

      boxplot(phat_sim, ylim=range(phat_true, phat_sim, na.rm=TRUE),
              main="p_hat", xlab="Age",
              border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
      points(phat_true, col=2, pch=16)
      legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
             fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))

      boxplot(se_phat_sim, ylim=range(se_phat_true, se_phat_sim, na.rm=TRUE),
              main="se(p_hat)", xlab="Age",
              border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
      points(se_phat_true, col=2, pch=16)
      legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
             fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
    }

    if(!is.null(sim$Nhat_sim)) {
      boxplot(Nhat_sim, ylim=range(Nhat_true, Nhat_sim, na.rm=TRUE),
              main="N_hat", xlab="Age",
              border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
      points(Nhat_true, col=2, pch=16)
      legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
             fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))

      boxplot(se_Nhat_sim, ylim=range(se_Nhat_true, se_Nhat_sim, na.rm=TRUE),
              main="se(N_hat)", xlab="Age",
              border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
      points(se_Nhat_true, col=2, pch=16)
      legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
             fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
    }

    if(!is.null(sim$mn_length_sim)) {
      boxplot(mn_length_sim, ylim=range(mn_length_true, mn_length_sim, na.rm=TRUE),
              main="mn_length", xlab="Age",
              border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
      points(mn_length_true, col=2, pch=16)
      legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
             fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))

      boxplot(se_mn_length_sim, ylim=range(se_mn_length_true, se_mn_length_sim, na.rm=TRUE),
              main="se(mn_length)", xlab="Age",
              border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
      points(se_mn_length_true, col=2, pch=16)
      legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
             fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
    }
  })
}





#' Relative precision for the Methods Used in `ASL_table()`
#' @description This function may be used with a simulation object created by
#' \link{simulate_ASL_table}, or will call it internally if none is supplied.
#'
#' The simulation function draws many samples from a population, given
#' population and sample characteristics supplied by the user, or else one of
#' twelve possible cases designed to cover all methods used in `ASL_table()`.
#'
#' This function uses the simulation results to produce a sequence of plots
#' displaying the relationships between relative precision (how close we want
#' estimates to be to the true value(s)) and confidence (how often we can expect
#' estimates to be within the desired relative precision of the truth).
#'
#' In terms of a typical Operational Plan Objective: "Estimate ... such that all
#' estimates are within (precision) of the true values (confidence) % of the time".
#'
#' Relative precision curves will be created for proportions (phat), abundance (Nhat),
#' and mean length (mn_length), depending on which are present.  Multiple curves
#' will be overlayed for each, if multiple categories are present for Age and/or Sex.
#' @param conf_target The target level of confidence desired.  Defaults to `0.9`.
#' @param sim An object created by \link{simulate_ASL_table}.  If the default (`NULL`)
#' is accepted, the simulation function will be called internally from the other
#' arguments provided.
#' @param case If a pre-determined case is to be used.  Allowed values are
#' `"stratified_witherror_lengthage"`, `"stratified_witherror_age"`, `"stratified_witherror_length"`,
#' `"stratified_lengthage"`, `"stratified_age"`, `"stratified_length"`,
#' `"stratified_Nunknown_lengthage"`, `"stratified_Nunknown_age"`, `"stratified_Nunknown_length"`,
#' `"pooled_witherror_lengthage"`, `"pooled_witherror_age"`, `"pooled_witherror_length"`
#' `"pooled_lengthage"`, `"pooled_age"`,  `"pooled_length"`,
#' `"pooled_Nunknown_lengthage"`, `"pooled_Nunknown_length"`, or `"pooled_Nunknown_age"`. If the default
#' (`NULL`) is accepted, all simulation parameters below must be supplied by the user.
#' @param nstrata Number of sampling strata
#' @param nage Number of age categories
#' @param nt Sample size for each stratum
#' @param Nt Abundance for each stratum.  If abundance is completely unknown, sampling
#' weights may be used instead in the `sampling_weights=` argument.
#' @param stratum_weights Optional vector of sampling weights for each stratum.  Defaults to `NULL`.
#' @param se_Nt Optional vector of standard errors of abundance for each stratum
#' Defaults to `NULL`, implying abundance is known without error.
#' @param mn_length Vector of mean lengths for each age category.  Defaults to
#' `NULL`, implying lengths were not considered.
#' @param sd_length Vector of standard deviations for each age category.  Defaults to
#' `NULL`, implying lengths were not considered.
#' @param ptz Matrix of probabilities of each age by stratum, with rows
#' corresponding to strata and columns corresponding to ages.  The probabilities
#' for each row will be normalized before simulation, therefore summing to one
#' is not required.  If a pooled (non-stratified) sample is to be taken and
#' age categories are to be considered, this should be supplied as a vector with
#' length equal to the number of ages.
#' @param nNA Number of NA values to randomly impute, to test robustness to NA.
#' Defaults to `0` in this case.
#' @param nsim Number of simulated replicates.  Defaults to `1000`, but more is recommended.
#' @param plot_pop Whether to make summary plots of the simulated population and
#' one representative sample, in addition to the plots produced in simulation.
#' Defaults to `TRUE`.
#' @param verbose Whether to print the parameters used in simulation to the console,
#' if one of the `case`s is accepted, as well as printing the method used within
#' `ASL_table()`.  Defaults to `TRUE`.
#' @param print_table Whether to print an example output table from `ASL_table()`
#' as an additional check.  Defaults to `FALSE`.
#' @return `NULL`
#' @author Matt Tyers
#' @seealso [ASL_table]
#' @importFrom graphics abline lines
#' @examples
#' # running rp_ using a pre-defined case
#' par(mfrow=c(2,2))
#' rp_ASL_table(case="stratified_witherror_lengthage",
#'              nsim=500, plot_pop=FALSE)
#'
#' # running rp_ directly from simulation parameters
#' par(mfrow=c(2,2))
#' rp_ASL_table(nstrata = 4,
#'              nage = 5,
#'
#'              # sample size for each stratum
#'              nt = c(100, 100, 100, 100),
#'
#'              # abundance for each stratum
#'              Nt = c(10000, 20000, 30000, 40000),
#'
#'              # (possible) SE for abundance by stratum
#'              se_Nt = 0.2*c(10000, 20000, 30000, 40000),
#'
#'              # mean length FOR EACH AGE
#'              mn_length = c(150, 200, 250, 300, 350),
#'
#'              # sd length FOR EACH AGE
#'              sd_length = c(40, 50, 60, 70, 80),
#'
#'              # matrix of probabilities of each age BY stratum
#'              ptz = matrix(c(c(1,2,3,4,5),
#'                             c(1,2,5,5,2),
#'                             c(2,5,3,2,1),
#'                             c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5),
#'              nsim=500, plot_pop=FALSE)
#'
#'
#' ## creating a simulation object first
#' simresults <- simulate_ASL_table(nstrata = 4,
#'                                  nage = 5,
#'
#'                                  # sample size for each stratum
#'                                  nt = c(100, 100, 100, 100),
#'
#'                                  # abundance for each stratum
#'                                  Nt = c(10000, 20000, 30000, 40000),
#'
#'                                  # (possible) SE for abundance by stratum
#'                                  se_Nt = 0.2*c(10000, 20000, 30000, 40000),
#'
#'                                  # mean length FOR EACH AGE
#'                                  mn_length = c(150, 200, 250, 300, 350),
#'
#'                                  # sd length FOR EACH AGE
#'                                  sd_length = c(40, 50, 60, 70, 80),
#'
#'                                  # matrix of probabilities of each age BY stratum
#'                                  ptz = matrix(c(c(1,2,3,4,5),
#'                                                 c(1,2,5,5,2),
#'                                                 c(2,5,3,2,1),
#'                                                 c(5,4,3,1,1)), byrow=TRUE, nrow=4, ncol=5),
#'                                  plot_pop=FALSE)
#'
#' ## running rp_ with the object created
#' par(mfrow=c(2,2))
#' rp_ASL_table(sim=simresults)
#'
#' ## running rp_ again
#' par(mfrow=c(2,2))
#' rp_ASL_table(sim=simresults, conf_target=0.99)
#' @export
rp_ASL_table <- function(conf_target=0.9,
                         sim=NULL,
                         case=NULL,   # should this default to NULL?
                         nstrata,
                         nage,
                         nt,   # c(100, 100, 100, 100), # sample size for each stratum
                         Nt,   # c(10000, 20000, 30000, 40000),  # abundance for each stratum
                         stratum_weights = NULL,
                         se_Nt = NULL,   # 0.2*Nt, # (possible) SE for abundance by stratum,
                         mn_length = NULL,   # c(150, 200, 250, 300, 350), # mean length FOR EACH AGE
                         sd_length = NULL,   # c(40, 50, 60, 70, 80), # sd length FOR EACH AGE
                         ptz = NULL,   # matrix(c(c(1,2,3,4,5),  # matrix of probabilities of each age BY stratum
                         # c(1,2,5,5,2),
                         # c(2,5,3,2,1),
                         # c(5,4,3,1,1)),
                         # byrow=TRUE,
                         # nrow=4,
                         # ncol=5),
                         nNA=0,
                         nsim=1000,   # 1000,    # number of simulated replicates
                         plot_pop=TRUE,   # whether to make summary plots of pop & one sample
                         verbose=TRUE, # whether to output cases in sim function
                         print_table=FALSE) {

  ## if sim= is supplied, make sure it is actually a simulation object!
  if(!is.null(sim)) {
    if(!inherits(sim, "ASL_sim")) {
      stop("Argument sim= must be an object returned from simulate_ASL_table()")
    }
  }

  ## if a pre-run simulation is not supplied, then run the simulation
  if(is.null(sim)) {

    # validate (class match?) that sim is useable

    sim <- simulate_ASL_table(case=case,
                              nstrata=nstrata,
                              nage=nage,
                              nt=nt,
                              Nt=Nt,
                              stratum_weights = stratum_weights,
                              se_Nt=se_Nt,
                              mn_length=mn_length,
                              sd_length=sd_length,
                              ptz=ptz,
                              nNA=nNA,
                              nsim=nsim,
                              plot_pop=plot_pop,
                              verbose=verbose,
                              print_table=print_table)
  }

  # if phat is part of the simulation results
  if(!is.null(sim$phat_sim)) {
    # calculate max d (x-axis, relative precision) for plotting
    dmax <- max(abs(sim$phat_sim - matrix(sim$phat_true,
                                          byrow=TRUE,
                                          nrow=nrow(sim$phat_sim),
                                          ncol=ncol(sim$phat_sim))))*100
    # x 100 to get on percentage points scale

    dd <- seq(0, dmax, length.out=1000)  # vector of x-axis values to consider

    # calculating a vector of confidence (how often) to accompany relative prec values
    aa_phat <- matrix(nrow=length(dd), ncol=length(sim$phat_true))
    for(i in 1:length(sim$phat_true)) {
      aa_phat[,i] <- colMeans(outer(abs(sim$phat_sim[,i] - sim$phat_true[i]), dd/100, "<="))
    }

    # calculating a vector of relative precisions to go with target confidence level (one for each category)
    dd_target <- dd[apply(aa_phat, 2, function(x) which.max(x>=conf_target))]
    percmain <- ifelse(all(dd_target > 0),
                       paste0(100*conf_target,
                              "% within ",
                              round(max(dd_target), 1),
                              " perc points of true"),"")
    plot(NA, xlim=c(0,dmax), ylim=0:1,
         main=c("- phat -", percmain),
         ylab="Confidence", xlab="Rel precision (perc points of true value)")
    for(i in 1:length(sim$phat_true)) {
      lines(dd, aa_phat[,i], col=4)
    }
    abline(v=dd_target, h=conf_target, col=2)
  }

  # if Nhat is part of the simulation results
  if(!is.null(sim$Nhat_sim)) {
    Nhat_true_mat <- matrix(sim$Nhat_true,
                            byrow=TRUE,
                            nrow=nrow(sim$Nhat_sim),
                            ncol=ncol(sim$Nhat_sim))
    # calculate max d (x-axis, relative precision) for plotting
    dmax <- max(abs(sim$Nhat_sim - Nhat_true_mat)/Nhat_true_mat)
    dd <- seq(0, dmax, length.out=1000)  # vector of x-axis values to consider

    aa_Nhat <- matrix(nrow=length(dd), ncol=length(sim$Nhat_true))
    for(i in 1:length(sim$Nhat_true)) {
      aa_Nhat[,i] <- colMeans(outer(abs(sim$Nhat_sim[,i] - sim$Nhat_true[i])/sim$Nhat_true[i], dd, "<="))
    }

    # calculating a vector of relative precisions to go with target confidence level (one for each category)
    dd_target <- dd[apply(aa_Nhat, 2, function(x) which.max(x>=conf_target))]
    percmain <- ifelse(all(dd_target > 0),
                       paste0(100*conf_target,
                              "% within ",
                              round(100*max(dd_target), 1),
                              "% of true"),"")
    plot(NA, xlim=c(0,dmax), ylim=0:1,
         main=c("- Nhat -", percmain),
         ylab="Confidence", xlab="Rel precision (prop of true value)")
    for(i in 1:length(sim$Nhat_true)) {
      lines(dd, aa_Nhat[,i], col=4)
    }
    abline(v=dd_target, h=conf_target, col=2)
  }


  # if mn_length is part of the simulation results
  if(!is.null(sim$mn_length_sim)) {
    mn_length_true_mat <- matrix(sim$mn_length_true,
                                 byrow=TRUE,
                                 nrow=nrow(sim$mn_length_sim),
                                 ncol=ncol(sim$mn_length_sim))
    # calculate max d (x-axis, relative precision) for plotting
    dmax <- max(abs(sim$mn_length_sim - mn_length_true_mat)/mn_length_true_mat)
    dd <- seq(0, dmax, length.out=1000)  # vector of x-axis values to consider

    aa_mn_length <- matrix(nrow=length(dd), ncol=length(sim$mn_length_true))
    for(i in 1:length(sim$mn_length_true)) {
      aa_mn_length[,i] <- colMeans(outer(abs(sim$mn_length_sim[,i] - sim$mn_length_true[i])/sim$mn_length_true[i], dd, "<="))
    }

    # calculating a vector of relative precisions to go with target confidence level (one for each category)
    dd_target <- dd[apply(aa_mn_length, 2, function(x) which.max(x>=conf_target))]
    percmain <- ifelse(all(dd_target > 0),
                       paste0(100*conf_target,
                              "% within ",
                              round(100*max(dd_target), 1),
                              "% of true"),"")
    plot(NA, xlim=c(0,dmax), ylim=0:1,
         main=c("- mean length -", percmain),
         ylab="Confidence", xlab="Rel precision (prop of true value)")
    for(i in 1:length(sim$mn_length_true)) {
      lines(dd, aa_mn_length[,i], col=4)
    }
    abline(v=dd_target, h=conf_target, col=2)
  }
}

