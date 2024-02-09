#' Construct an ASL Table from Vectors of Data
#' @description This function is intended to perform all the necessary calculations
#' to construct tables summarizing ASL (Age, Sex, and Length) for several data
#' scenarios:
#'
#' * Stratified vs. non-stratified
#' * Abundance known without error vs. estimated with error vs. unknown
#' * Data available for Age, Sex, and Length vs. some subset.
#'
#' The function will return a table with rows for each combination of Age and
#' Sex (as available), and will summarize the respective proportions, estimated
#' abundance, and associated Lengths (as available).
#'
#' Generally, `NULL` values in a given data vector indicate that the vector will
#' not be used.
#'
#' @param age Vector of ages.  If the default (`NULL`) is accepted, age will be
#' omitted from the resulting table.
#' @param sex Vector of sex.  If the default (`NULL`) is accepted, sex will be
#' omitted from the resulting table.
#' @param length Vector of length.  If the default (`NULL`) is accepted, length
#' will be omitted from the resulting table.
#' @param stratum Optional vector of stratum, if stratified estimators will be
#' used.  This must be formatted as positive whole numbers (1, 2, 3, etc).
#' Defaults to `NULL`.
#' @param Nhat Optional vector of abundance for each stratum, if stratified
#' estimators will be used.  The length of this vector must correspond to the
#' maximum value in the `stratum=` argument.  Defaults to `NULL`.
#' @param se_Nhat Optional vector of the standard error of abundance for each
#' stratum, if stratified estimators will be used and by-stratum abundance is
#' considered to have measurement error.    The length of this vector must
#' correspond to the maximum value in the `stratum=` argument.  Defaults to `NULL`.
#' @param stratum_weights Optional vector of weights for each stratum, if
#' relative weights are known but abundance is not.  Using this argument rather
#' than `Nhat` will also omit the Finite Population Correction factor (FPC) from
#' calculations, but weights are treated as constant (without error), likely
#' underestimating variance.  Defaults to `NULL`.
#' @param verbose Whether to print messages corresponding to the method used.
#' Defaults to `FALSE`.
#' @param FPC Whether to incorporate the Finite Population Correction factor (FPC) in
#' variance calculations.  Allowed values are
#' * `"ifknown"` (the default), which will use the FPC only if abundance is known without error
#' * `"always"`, which will use the FPC wherever possible (i.e. if there is an estimate of abundance)
#' * `"never"`, which will always ignore the FPC.
#' used if the abundance is considered to be known without error.
#' @return A data.frame with rows corresponding to categories of age and/or sex,
#' depending on data inputs.
#'
#' ## Stratified - If abundance is known without error
#'
#' ### If there are proportions
#'
#' The proportion of each age and/or sex category *z* will be estimated for each sampling stratum *t* as follows:
#'
#'  \deqn{\hat{p}_{tz}=\frac{n_{tz}}{n_t}}
#'
#' in which \eqn{n_{tz}} equals the number of fish sampled during sampling stratum \eqn{t} classified as age and/or sex category \eqn{z}, and \eqn{n_t} equals the number of fish sampled for age and/or sex determination within sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\hat{p}_{tz}} will be estimated as the following (Cochran 1977):
#'
#' \deqn{\hat{var}[\hat{p}_{tz}]=\frac{\hat{p}_{tz}(1-\hat{p}_{tz})}{n_t-1}\left(\frac{N_t-n_t}{N_t-1}\right)}
#'
#' if \eqn{N_t}, the total abundance of fish in sampling stratum \eqn{t}, is known and the finite population correction factor (FPC) is used; otherwise, as the following:
#'
#' \deqn{\hat{var}[\hat{p}_{tz}]=\frac{\hat{p}_{tz}(1-\hat{p}_{tz})}{n_t-1}}
#'
#' The total abundance by age and/or sex category in each sampling stratum will be estimated as follows:
#'
#' \deqn{\hat{N}_{tz}=N_t\hat{p}_{tz}}
#'
#' with variance estimated as
#'
#' \deqn{\hat{var}[\hat{N}_{tz}]=N_t^2\hat{var}[\hat{p}_{tz}]}
#'
#' The total abundance by age and/or sex category and its variance will then be estimated by summation as follows:
#'
#' \deqn{\hat{N}_z=\sum_{t=1}^{L}\hat{N}_{tz}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{N}_{z}]=\sum_{t=1}^{L}\hat{var}[\hat{N}_{tz}]}
#'
#' where \eqn{L} equals the number of sampling strata.
#'
#' Finally, the overall proportion by age and/or sex category and its variance will be estimated as follows:
#'
#' \deqn{\hat{p}_z=\frac{\hat{N}_z}{N}}
#'
#' and
#'
#' \deqn{\hat{var}[\hat{p}_z]=\frac{\hat{var}[\hat{N}_z]}{N^2}}
#'
#' where \eqn{N} is the total abundance across all sampling periods.
#'
#' The mean length by age and/or sex for each sampling stratum will be estimated as follows:
#'
#' \deqn{\bar{x}_{tz}=\frac{\sum_{i=1}^{n_{tz}}x_{tzi}}{n_{tz}}}
#'
#' where \eqn{x_{tzi}} is the length of the *i*th fish sampled of age and/or sex \eqn{z} during sampling stratum \eqn{t}.
#'
#' The sampling variance of \eqn{\bar{x}_{tz}} will be estimated as
#'
#' \deqn{\hat{var}[\bar{x}_{tz}]=\frac{\sum_{i=1}^{n_{tz}}(x_{tzi}-\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}\left(\frac{
#'     \hat{N}_{tz}-n_{tz}}{\hat{N}_{tz}-1}\right)}
#'
#' if the finite population correction factor (FPC) will be used; otherwise, as the following:
#'
#' \deqn{\hat{var}[\bar{x}_{tz}]=\frac{\sum_{i=1}^{n_{tz}}(x_{tzi}-\bar{x}_{tz})^2}{n_{tz}(n_{tz}-1)}}
#'
#' The mean length by age and/or sex category will then be estimated as follows:
#'
#' \deqn{\bar{x}_z=\sum_{t=1}^{L}\frac{\hat{N}_{tz}}{\hat{N}_z}\bar{x}_{tz}}
#'
#' with its variance approximated using a Taylor's series expansion (Mood et al. 1974):
#'
#' \deqn{\hat{var}[\bar{x}_z]\approx\sum_{t=1}^{L}\frac{\hat{N}_{tz}^2}{\hat{N}_z^2}\hat{var}[\bar{x}_{tz}]+\sum_{t=1}^{L}\frac{\left(\bar{x}_{tz}\hat{N}_z-\left(\sum_{u=1}^{L}\bar{x}_{uz}\hat{N}_{uz}\right)\right)^2}{\hat{N}_z^4}\hat{var}[\hat{N}_{tz}]}
#'
#' ### If there are no proportions to estimate
#'
#' The mean length for each sampling stratum will be estimated as follows, where \eqn{x_{ti}} is the length of the *i*th fish sampled within sampling stratum \eqn{t}, and \eqn{n_t} is the number of fish in stratum *t* sampled for length:
#'
#' \deqn{\bar{x}_{t}=\frac{\sum_{i=1}^{n_{t}}x_{ti}}{n_{t}}}
#'
#' The sampling variance of \eqn{\bar{x}_{t}} will be estimated as
#'
#' \deqn{\hat{var}[\bar{x}_{t}]=\frac{\sum_{i=1}^{n_{t}}(x_{ti}-\bar{x}_{t})^2}{n_{t}(n_{t}-1)}\left(\frac{
#' N_{t}-n_{t}}{N_{t}-1}\right)}
#'
#' if abundance per stratum \eqn{N_t} is known and if the finite population correction factor is used, otherwise as:
#'
#' \deqn{\hat{var}[\bar{x}_{t}]=\frac{\sum_{i=1}^{n_{t}}(x_{ti}-\bar{x}_{t})^2}{n_{t}(n_{t}-1)}}
#'
#' Stratified estimates of mean length will be calculated as follows, in which \eqn{N_t} and \eqn{\bar{x}_t} represent the abundance and average length associated with sampling stratum \eqn{t}, respectively, and \eqn{N} represents the total abundance:
#'
#' \deqn{\bar{x}=\frac{1}{N}\sum_{t=1}^L N_t\bar{x}_t}
#'
#' and
#'
#' \deqn{\hat{var}[\bar{x}]=\sum_{t=1}^L\left(\frac{N_t}{N}\right)^2 \hat{var}[\bar{x}_t]}
#'
#' @author Matt Tyers
#' @seealso [verify_ASL_table]
#' @examples
#' ### a look at input
#' head(sim_data$data)
#' sim_data$abundance
#'
#' ### tables for possible scenarios
#'
#' # Stratified with error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           stratum=sim_data$data$stratum,
#'           Nhat=sim_data$abundance$Nhat,
#'           se_Nhat=sim_data$abundance$se_Nhat)
#'
#' # Stratified without error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           stratum=sim_data$data$stratum,
#'           Nhat=sim_data$abundance$Nhat)
#'
#' # Pooled (not stratified) with error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           Nhat=100*1000,
#'           se_Nhat=10*1000)
#'
#' # Pooled (not stratified) without error in Nhat
#' ASL_table(age=sim_data$data$age,
#'           length=sim_data$data$length,
#'           Nhat=100*1000)
#' @importFrom grDevices adjustcolor grey.colors
#' @importFrom graphics boxplot legend mosaicplot par points
#' @importFrom stats rnorm sd var
#' @export
ASL_table <- function(age=NULL,
                      sex=NULL,
                      length=NULL,
                      stratum=NULL,
                      Nhat=NULL,
                      se_Nhat=NULL,
                      stratum_weights = NULL,
                      verbose=FALSE,
                      FPC=c("ifknown", "always", "never")) {

  # -------- globally dealing with inputs ----------
  # combining age/sex categories as available
  if(!is.null(sex) & !is.null(age)) {
    cats <- paste(sex, age)
  } else {
    if(!is.null(sex)) cats <- sex
    if(!is.null(age)) cats <- age
    if(is.null(sex) & is.null(age)) cats <- rep("Total", length(length))
  }

  if(!is.null(stratum)) {
    # taking out data rows where strata value is NA
    length <- length[!is.na(stratum)]
    cats <- cats[!is.na(stratum)]
    stratum <- stratum[!is.na(stratum)]
  }

  ### need to constrain FPC to c("ifknown","always","never")
  FPC <- match.arg(FPC)

  # this will become output dataframe
  out <- data.frame(n=tapply(X=seq_along(cats), INDEX=cats, FUN=base::length))

  # --------- Stratified ----------
  if(!is.null(stratum)) {
    if(verbose) cat("\n", "Stratified","\n")
    # checking that inputs make sense

    # if(!inherits(stratum, "integer")) {
    #   stop("Stratum values must be positive whole numbers (1, 2, 3, etc.)")  ############### this fails
    # }


    if(length(stratum_weights) != max(stratum) & length(Nhat) != max(stratum)) {
      stop("Length of stratum weights does not equal max stratum value")
    }

    # -- Stratified without error in Nhats --
    if(is.null(se_Nhat)) {
      if(verbose) cat("\n", "without error in Nhat", "\n")
      # if(is.na(FPC)) {
      #   FPC <- !is.null(Nhat) #& is.null(stratum_weights)
      # }
      if(is.null(Nhat) & !is.null(stratum_weights)) {
        Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
        # FPC <- FALSE
      } else {
        Nt <- Nhat # consistency with report eqns
      }

      # pulling these equations from the Jim Creek report
      if(length(unique(cats)) > 1) {

        # --- if there are proportions to estimate ---
        ntz <- table(stratum, cats, useNA="no")
        nt <- rowSums(ntz)
        if((is.null(Nhat) & !is.null(stratum_weights)) | FPC=="never") {
          FPC_vec <- 1#rep(1, length(nt))  # this will ignore FPC in further calculations
          if(verbose) cat("\n","Finite Population Correction factor NOT used for proportions","\n")
        } else {
          FPC_vec <- (Nt-nt)/(Nt-1)
          if(verbose) cat("\n","Finite Population Correction factor USED for proportions","\n")
        }
        ptz <- ntz/nt
        vptz <- ptz*(1-ptz)/(nt-1)*FPC_vec
        Ntz <- Nt*ptz
        vNtz <- (Nt^2)*vptz
        Nz <- colSums(Ntz)
        vNz <- colSums(vNtz)
        pz <- Nz/sum(Nt)
        vpz <- vNz/(sum(Nt)^2)
        out$phat <- pz
        out$se_phat <- sqrt(vpz)
        if(!is.null(Nhat)) {
          out$Nhat <- Nz
          out$se_Nhat <- sqrt(vNz)
        }

        # summarizing length
        if(!is.null(length)) {
          xbartz <- tapply(length, list(stratum, cats), mean, na.rm=TRUE)   ### this should have a na.rm=TRUE !!!
          vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          if((!is.null(Nhat) & is.null(stratum_weights)) & FPC=="always") {  ############
            FPC_vec2 <- (Ntz-vxbartz_denom)/(Ntz-1)  ####### new FPC here
            if(verbose) cat("\n","Finite Population Correction factor USED for means","\n")
          } else {
            FPC_vec2 <- 1#rep(1, length(nt))
            if(verbose) cat("\n","Finite Population Correction factor NOT used for means","\n")
          }
          vxbartz <- tapply(length, list(stratum, cats), var, na.rm=TRUE)/vxbartz_denom * FPC_vec2

          Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=T)

          # Taylor series expansion (Mood et al. 1974)
          prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
            colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)

          # ## estimate xbartz using my ugly derivation
          # vsumNtzx <- colSums((Ntz^2)*vxbartz + (xbartz^2)*vNtz - vNtz*vxbartz)
          # #vNz is calculated above
          # covthing <- colSums(xbartz*vNtz)
          # vxbarz <- ((colSums(Ntz*xbartz)/Nz)^2)*(vsumNtzx/(colSums(Ntz*xbartz)^2) + vNz/(Nz^2) - 2*covthing/(colSums(Ntz*xbartz)*Nz))

          out$n_length <- tapply(!is.na(length), cats, sum)
          out$mn_length <- xbarz
          out$se_length <- sqrt(vxbarz)
          out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          out$max_length <- tapply(length, cats, max, na.rm=TRUE)
        }
      } else {
        # --- if there are NO proportions to estimate ---
        nt <- table(stratum, is.na(length))[,1] # table(stratum)   ### now excludes NA values in length
        if((is.null(Nhat) | !is.null(stratum_weights)) | FPC=="never") {
          if(verbose) cat("\n","Finite Population Correction factor NOT used for means", "\n")
          FPC_vec <- 1#rep(1, length(nt))  # this will ignore FPC in further calculations
        } else {
          if(verbose) cat("\n","Finite Population Correction factor USED for means", "\n")
          FPC_vec <- (Nt-nt)/(Nt-1)
        }
        if(!is.null(stratum_weights)) {
          Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
        } else {
          Nt <- Nhat # consistency with report eqns
        }
        xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt*FPC_vec    #### moved FPC_vec to here
        out$n_length <- sum(!is.na(length))
        out$mn_length <- sum(Nt*xbart/sum(Nt))
        out$se_length <- sqrt(sum(((Nt/sum(Nt))^2)*vxbart))    ####### from here
        out$min_length <- min(length, na.rm=TRUE)
        out$max_length <- max(length, na.rm=TRUE)
      }
    } else {
      # -- Stratified with error in Nhats --
      # estimating proportions..
      # estimating abundance
      # summarizing length
      if(verbose) cat("\n", "WITH error in Nhat", "\n")
      # if(is.na(FPC)) {
      #   FPC <- FALSE
      # }
      Nt <- Nhat # consistency with report eqns
      # modifying these equations from the Jim Creek report

      if(length(unique(cats)) > 1) {
        # --- if there are proportions to estimate ---
        ntz <- table(stratum, cats, useNA="no")
        nt <- rowSums(ntz)  # this will be robust to the presence of NA
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          stop("Nhat must be supplied when se_Nhat is used.")
        } else {
          Nt <- Nhat # consistency with report eqns
          # FPC_vec <- (Nt-nt)/(Nt-1)
        }
        if(is.null(stratum_weights) & FPC=="always") {
          FPC_vec <- (Nt-nt)/(Nt-1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for proportions", "\n")
        } else {
          FPC_vec <- 1#rep(1, length(Nt))
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for proportions", "\n")
        }
        ptz <- ntz/nt
        vptz <- ptz*(1-ptz)/(nt-1)*FPC_vec
        Ntz <- Nt*ptz
        vNtz <- ((Nt^2)*vptz) + ((ptz^2)*(se_Nhat^2)) - (vptz*(se_Nhat^2)) ### Goodman
        Nz <- colSums(Ntz)
        vNz <- colSums(vNtz)
        pz <- Nz/sum(Nt)

        # delta method
        # covNzSum <- vNz  ## this is a minimum, there should be more
        covNzSum <- colSums((se_Nhat^2)*ptz, na.rm=TRUE)  # this seems unbiased!!

        vpz <- ((Nz/sum(Nz))^2)*(vNz/(Nz^2) + sum(se_Nhat^2)/((sum(Nz))^2) - 2*covNzSum/(Nz*sum(Nz)))
        out$phat <- pz
        out$se_phat <- sqrt(vpz)
        # if(any(is.na(out$se_phat))) {
        #   print(ntz)
        #   print(vpz)
        # }

        out$Nhat <- Nz
        out$se_Nhat <- sqrt(vNz)   ## this one seems unbiased still

        # summarizing length
        if(!is.null(length)) {
          xbartz <- tapply(length, list(stratum, cats), mean, na.rm=TRUE)   ### this should have a na.rm=TRUE !!!
          vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          if(is.null(stratum_weights) & FPC=="always") {
            FPC_vec2 <- (Ntz-vxbartz_denom)/(Ntz-1)  ####### new FPC here
            if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
          } else {
            FPC_vec2 <- 1#rep(1, length(nt))
            if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
          }
          vxbartz <- tapply(length, list(stratum, cats), var, na.rm=TRUE)/vxbartz_denom * FPC_vec2

          Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=TRUE)
          # Taylor series expansion (Mood et al. 1974)
          prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
            colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)
          out$n_length <- tapply(!is.na(length), cats, sum)
          out$mn_length <- xbarz
          out$se_length <- sqrt(vxbarz)
          out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          out$max_length <- tapply(length, cats, max, na.rm=TRUE)
        }
      } else {
        # --- if there are NO proportions to estimate ---
        nt <- table(stratum, is.na(length))[,1] # table(stratum)   ### now excludes NA values in length
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          stop("Nhat must be supplied when se_Nhat is used.")
        }
        if(is.null(stratum_weights) & FPC=="always") {
          FPC_vec <- (Nt-nt)/(Nt-1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
        } else {
          FPC_vec <- 1#rep(1, length(Nt))
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
        }
        Nt <- Nhat # consistency with report eqns
        xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt*FPC_vec ##################### include FPC_vec???
        out$n_length <- sum(nt)
        out$mn_length <- sum(Nt*xbart/sum(Nt))
        # NoverSumN <- Nt/sum(Nt)
        # # from delta method  -- validated by simulation   # currently no fpc anywhere
        # vNoverSumN <- (NoverSumN^2)*(((se_Nhat^2)/(Nt^2)) +
        #                                (sum(se_Nhat^2)/(sum(Nt)^2)) -
        #                                (2*(se_Nhat^2)/(Nt*sum(Nt))))
        # # from Goodman
        # out$se_length <- sqrt(sum(((NoverSumN^2)*vxbart) +   # vxbart validated by simulation
        #                             ((xbart^2)*vNoverSumN) -
        #                             (vxbart*vNoverSumN))) #### there is huge covariance in the sum!!
        # from covariance identities
        covSumNtSum <- sum(xbart*(se_Nhat^2))
        # from Goodman
        vSumNtXbart <- sum(((Nt^2)*vxbart) +
                             (xbart^2)*(se_Nhat^2) -
                             ((se_Nhat^2)*vxbart))
        # from delta method
        out$se_length <- sqrt(((sum(Nhat*xbart)/sum(Nhat))^2)*((vSumNtXbart/((sum(Nhat*xbart))^2)) +
                                                                 (sum(se_Nhat^2)/((sum(Nhat))^2)) -
                                                                 ((2*covSumNtSum)/(sum(Nhat)*sum(Nhat*xbart)))))
        out$min_length <- min(length, na.rm=TRUE)
        out$max_length <- max(length, na.rm=TRUE)
      }
    }
  }
  # --------- NOT Stratified---------
  if(is.null(stratum)) {
    if(verbose) cat("\n", "not stratified", "\n")

    # if(is.na(FPC)) {
    #   FPC <- !is.null(Nhat) & is.null(se_Nhat)
    # }
    # if(FPC & (is.null(Nhat) )) {  #WRONG  | !is.null(se_Nhat)
    #   FPC <- FALSE
    # }
    # FPC <- ifelse(!is.null(Nhat), (Nhat-sum(out$n))/(Nhat-1), 1)
    # FPC_mn <- FPC   # FPC used for variance of mean length

    # checking that inputs make sense
    if(length(Nhat) > 1 | length(se_Nhat) > 1 | length(stratum_weights) > 1) {
      stop("Stratum totals or weights are given, but no strata")
    }

    # estimating proportions..
    if(length(unique(cats)) > 1) {

      if(!is.null(Nhat) & (FPC=="always" | (FPC!="never" & is.null(se_Nhat)))) {
        FPC_prop <- (Nhat-sum(out$n))/(Nhat-1)
        # FPC_mn <- (Nhat-sum(out$n))/(Nhat-1)
        if(verbose) cat("\n", "Finite Population Correction factor USED for proportions", "\n")
      } else {
        FPC_prop <- 1
        # FPC_mn <- 1
        if(verbose) cat("\n", "Finite Population Correction factor NOT used for proportions", "\n")
      }

      out_ntot <- sum(out$n)
      out$phat <- out$n/out_ntot
      out$se_phat <- sqrt(out$phat*(1-out$phat)/(out_ntot-1)*FPC_prop)     ### added fpc

      # estimating abundance..
      if(!is.null(Nhat)) {

        out$Nhat <- Nhat * out$phat
        if(is.null(se_Nhat)) {
          if(verbose) cat("\n", "without error in Nhat", "\n")
          out$se_Nhat <- Nhat * out$se_phat   # verify this??
        } else {
          # Goodman's formula
          if(verbose) cat("\n", "with error in Nhat", "\n")
          out$se_Nhat <- sqrt(((Nhat^2)*(out$se_phat^2)) +
                                ((out$phat^2)*(se_Nhat^2)) -
                                ((out$se_phat^2)*(se_Nhat^2)))
        }

        # updating FPC for mean if there are categories AND est abundance
        # if(FPC=="always") {
        #   FPC_mn <- (out$Nhat - out$n)/(out$Nhat - 1)
        # }
      }
    }

    # and summarizing lengths..
    if(!is.null(length)) {

      # populating FPC
      if(length(unique(cats)) > 1) {  # if there are proportions
        if(!is.null(Nhat) & FPC=="always") {
          FPC_mn <- (out$Nhat - out$n)/(out$Nhat - 1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
        } else {
          FPC_mn <- 1
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
        }
      } else {  # no proportions
        if(!is.null(Nhat) & ((FPC!="never" & is.null(se_Nhat)) | (FPC=="always"))) {
          FPC_mn <- (Nhat-sum(out$n))/(Nhat-1)
          if(verbose) cat("\n", "Finite Population Correction factor USED for means", "\n")
        } else {
          FPC_mn <- 1
          if(verbose) cat("\n", "Finite Population Correction factor NOT used for means", "\n")
        }
      }

      n_notNA <- tapply(!is.na(length), cats, sum)
      out$n_length <- n_notNA
      out$mn_length <- tapply(length, cats, mean, na.rm=TRUE)
      out_sd_length <- tapply(length, cats, sd, na.rm=TRUE) * sqrt(FPC_mn)
        # sqrt((out$Nhat-n_notNA)/(out$Nhat-1)) ## NEW FPC    # sqrt(FPC)     ### added fpc
      out$se_length <- out_sd_length/sqrt(n_notNA)
      out$min_length <- tapply(length, cats, min, na.rm=TRUE)
      out$max_length <- tapply(length, cats, max, na.rm=TRUE)
    }
  }
  return(out)
}

# # simulating some fake data..
# nn <- 100
# sex <- sample(c("F","M"), nn, replace=T)
# age <- sample(c(1.1,1.2,1.3,2.1,2.2,2.3), nn, replace=T)
# length <- rnorm(nn, mean=700, sd=100)
#
# length[1] <- NA
#
# stratum <- sample(1:5, nn, replace=T)
# Nhat <- c(1000, 2000, 3000, 4000, 5000)
# ASL_table(
#   age=age,
#   # sex=sex,
#   length=length,
#   # Nhat=10000,
#   # se_Nhat=1000,
#   Nhat=Nhat,
#   se_Nhat=500,
#   # se_Nhat=Nhat/10,
#   # stratum_weights=Nhat/10,
#   stratum=stratum
# )




#' Empirical Verification of the Methods Used in `ASL_table()`
#' @description This function draws many samples from a population, given
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
#' Twelve pre-determined cases may be used in the `case=` argument, depending on
#' whether a stratified sampling scheme is used, whether abundance is known or
#' estimated with error, and whether estimates pertaining to categorical (age)
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
#' @param case If a pre-determined case is to be used.  Allowed values are
#' `"stratified_witherror_lengthage"`, `"stratified_witherror_age"`, `"stratified_witherror_length"`,
#' `"stratified_lengthage"`, `"stratified_age"`, `"stratified_length"`,
#' `"pooled_witherror_lengthage"`, `"pooled_witherror_age"`, `"pooled_witherror_length"`
#' `"pooled_lengthage"`, `"pooled_age"`, or `"pooled_length"`. If the default
#' (`NULL`) is accepted, all simulation parameters below must be supplied by the user.
#' @param nstrata Number of sampling strata
#' @param nage Number of age categories
#' @param nt Sample size for each stratum
#' @param Nt Abundance for each stratum
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
#' par(mfrow=c(2,2))
#' verify_ASL_table(case="stratified_witherror_lengthage", nsim=1000)
#'
#'
#' \dontrun{
#' nsim <- 5000
#' cases <- c("stratified_witherror_lengthage", "stratified_witherror_age",
#' "stratified_witherror_length", "stratified_lengthage", "stratified_age",
#' "stratified_length", "pooled_witherror_lengthage", "pooled_witherror_age",
#' "pooled_witherror_length", "pooled_lengthage", "pooled_age", "pooled_length")
#' for(case_i in cases) {
#'   par(mfrow=c(3,2))
#'   verify_ASL_table(case=case_i, nsim=5000)
#' }
#' }
#' @export
verify_ASL_table <- function(case=NULL,   # should this default to NULL?
                             nstrata,
                             nage,
                             nt,   # c(100, 100, 100, 100), # sample size for each stratum
                             Nt,   # c(10000, 20000, 30000, 40000),  # abundance for each stratum
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
  }

  # pre-fixing ptz
  if(is.null(dim(ptz)) & !is.null(ptz)) {
    ptz <- t(as.matrix(ptz))
  }

  # check inputs: length(nt) vs length(Nt) vs. length(se_Nt) vs nrow(ptz)
  if(length(nt)!=nstrata) stop("length(nt) should equal nstrata.")
  if(length(Nt)!=nstrata) stop("length(Nt) should equal nstrata.")
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

  # plotting population & one sample (optionally)
  if(plot_pop) {
    if(prod(dim(table(t,age_sim))) > 1) {  #
      mosaicplot(table(t,age_sim), xlab="Stratum", ylab="Age", main="Population", col=grey.colors(nage, rev=TRUE))
      mosaicplot(table(t[thesample],age_sim[thesample]), xlab="Stratum", ylab="Age", main="Sample", col=grey.colors(nage, rev=TRUE))
    }
    if(!is.null(length_sim)) {
      boxplot(length_sim~age_sim, xlab="Age", ylab="Length", main="Population")
      boxplot(length_sim[thesample]~age_sim[thesample], xlab="Age", ylab="Length", main="Sample")
      boxplot(length_sim~t, xlab="Stratum", ylab="Length", main="Population")
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
    # first take a sample of fish
    thesample <- sample(seq_along(t)[t==1], size=nt[1])
    if(length(nt) > 1) {
      for(i_t in 2:nstrata) {
        thesample <- c(thesample, sample(seq_along(t)[t==i_t], size=nt[i_t]))
      }
    }
    # table(t[thesample])  # making sure it worked

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
                          se_Nhat=se_Nt)  # find a way to add stratum_weights???
    # results[,,i_sim] <- as.matrix(thetable)
    results[i_sim,,] <- as.matrix(thetable)
  }
  # dimnames(results)[1:2] <- dimnames(thetable)
  dimnames(results)[2:3] <- dimnames(thetable)

  ### plotting simulation results, overlayed with true values
  if("phat" %in% dimnames(results)[[3]]) {
    phat_sim <- results[,,dimnames(results)[[3]]=="phat"]
    phat_true <- as.numeric(table(age_sim)/sum(Nt))
    boxplot(phat_sim, ylim=range(phat_true, phat_sim, na.rm=TRUE),
            main="p_hat", xlab="Age",
            border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
    points(phat_true, col=2, pch=16)
    legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
           fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))

    se_phat_sim <- results[,,dimnames(results)[[3]]=="se_phat"]
    se_phat_true <- apply(phat_sim, 2, sd, na.rm=TRUE)
    boxplot(se_phat_sim, ylim=range(se_phat_true, se_phat_sim, na.rm=TRUE),
            main="se(p_hat)", xlab="Age",
            border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
    points(se_phat_true, col=2, pch=16)
    legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
           fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
  }

  if("Nhat" %in% dimnames(results)[[3]]) {
    Nhat_sim <- results[,,dimnames(results)[[3]]=="Nhat"]
    Nhat_true <- as.numeric(table(age_sim))
    boxplot(Nhat_sim, ylim=range(Nhat_true, Nhat_sim, na.rm=TRUE),
            main="N_hat", xlab="Age",
            border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
    points(Nhat_true, col=2, pch=16)
    legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
           fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))

    se_Nhat_sim <- results[,,dimnames(results)[[3]]=="se_Nhat"]
    se_Nhat_true <- apply(Nhat_sim, 2, sd, na.rm=TRUE)
    boxplot(se_Nhat_sim, ylim=range(se_Nhat_true, se_Nhat_sim, na.rm=TRUE),
            main="se(N_hat)", xlab="Age",
            border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
    points(se_Nhat_true, col=2, pch=16)
    legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
           fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
  }

  if("mn_length" %in% dimnames(results)[[3]]) {
    mn_length_sim <- as.matrix(results[,,dimnames(results)[[3]]=="mn_length"])
    mn_length_true <- mn_length
    boxplot(mn_length_sim, ylim=range(mn_length_true, mn_length_sim, na.rm=TRUE),
            main="mn_length", xlab="Age",
            border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
    points(mn_length_true, col=2, pch=16)
    legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
           fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))

    se_mn_length_sim <- results[,,dimnames(results)[[3]]=="se_length"]
    se_mn_length_true <- apply(mn_length_sim, 2, sd, na.rm=TRUE)
    boxplot(se_mn_length_sim, ylim=range(se_mn_length_true, se_mn_length_sim, na.rm=TRUE),
            main="se(mn_length)", xlab="Age",
            border=adjustcolor(4,alpha.f=.6), col=adjustcolor(4,alpha.f=.2))
    points(se_mn_length_true, col=2, pch=16)
    legend("topright",legend=c("sim","true"),pch=c(NA,16),col=c(NA,2),
           fill=c(adjustcolor(4,alpha.f=.2),NA),border=c(adjustcolor(4,alpha.f=.6),NA))
  }
}


# par(mfrow=c(3,2))
# verify_ASL_table(nsim=50, plot_pop = T)   # stratified with error, age and length
# verify_ASL_table(nsim=50, plot_pop = T, mn_length=NULL, sd_length=NULL)   # stratified with error, just age
# verify_ASL_table(nsim=50, plot_pop = T, mn_length=300, sd_length=70, ptz=NULL)# stratified with error, just length
#
# verify_ASL_table(se_Nt = NULL, nsim=50, plot_pop = T)  # stratified, no error, age and length
# verify_ASL_table(se_Nt = NULL, nsim=50, plot_pop = T, mn_length=NULL, sd_length=NULL)   # stratified no error, just age
# verify_ASL_table(se_Nt = NULL, nsim=50, plot_pop = T, mn_length=300, sd_length=70, ptz=NULL)# stratified no error, just length
#
# verify_ASL_table(ptz=1:5, Nt=10000, nt=100, nsim=50, plot_pop = T)   # pooled with error, age and length
# verify_ASL_table(nsim=50, plot_pop = T, mn_length=NULL, sd_length=NULL, ptz=1:5, Nt=10000, nt=100)  # pooled with error, just age
# verify_ASL_table(nsim=50, plot_pop = T, mn_length=300, sd_length=70, ptz=NULL, Nt=10000, nt=100)# pooled with error, just length
#
# verify_ASL_table(se_Nt = NULL, ptz=1:5, Nt=10000, nt=100, nsim=50, plot_pop = T)   # pooled no error, age and length
# verify_ASL_table(se_Nt = NULL, nsim=50, plot_pop = T, mn_length=NULL, sd_length=NULL, ptz=1:5, Nt=10000, nt=100)  # pooled no error, just age
# verify_ASL_table(se_Nt = NULL, nsim=50, plot_pop = T, mn_length=300, sd_length=70, ptz=NULL, Nt=10000, nt=100)# pooled no error, just length

# nsim <- 5000
# par(mfrow=c(3,2))
# verify_ASL_table(case="stratified_witherror_lengthage", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="stratified_witherror_age", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="stratified_witherror_length", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="stratified_lengthage", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="stratified_age", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="stratified_length", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="pooled_witherror_lengthage", nsim=nsim, plot_pop = T)   # fails??
# par(mfrow=c(3,2))
# verify_ASL_table(case="pooled_witherror_age", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="pooled_witherror_length", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="pooled_lengthage", nsim=nsim, plot_pop = T)
# par(mfrow=c(3,2))
# verify_ASL_table(case="pooled_age", nsim=nsim, plot_pop = T)  # fails??  I think there are missing ages
# par(mfrow=c(3,2))
# verify_ASL_table(case="pooled_length", nsim=nsim, plot_pop = T)  # warnings again


#' Example data: sim_data
#'
#' A simulated dataset intended to illustrate \link{ASL_table}.  This dataset
#' is formatted as a list with two components:
#' * `sim_data$data` is a data frame with 400 observations of `$age`, `$length`,
#' and `$stratum`.
#' * `sim_data$abundance` is a data frame with 4 rows of `$Nhat` and `$se_Nhat`,
#' corresponding to each stratum.
"sim_data"


ASL_boilerplate <- function(age=NULL,
                      sex=NULL,
                      length=NULL,
                      stratum=NULL,
                      Nhat=NULL,
                      se_Nhat=NULL,
                      stratum_weights = NULL,
                      verbose=FALSE,
                      FPC=NA) {

  # -------- globally dealing with inputs ----------
  # combining age/sex categories as available
  if(!is.null(sex) & !is.null(age)) {
    cats <- paste(sex, age)
  } else {
    if(!is.null(sex)) cats <- sex
    if(!is.null(age)) cats <- age
    if(is.null(sex) & is.null(age)) cats <- rep("Total", length(length))
  }

  # if(!is.null(stratum)) {
  #   # taking out data rows where strata value is NA
  #   length <- length[!is.na(stratum)]
  #   cats <- cats[!is.na(stratum)]
  #   stratum <- stratum[!is.na(stratum)]
  # }
  #
  # # this will become output dataframe
  # out <- data.frame(n=tapply(X=seq_along(cats), INDEX=cats, FUN=base::length))

  # --------- Stratified ----------
  if(!is.null(stratum)) {
    # if(verbose) cat("\n", "Stratified","\n")
    # # checking that inputs make sense
    #
    # # if(!inherits(stratum, "integer")) {
    # #   stop("Stratum values must be positive whole numbers (1, 2, 3, etc.)")  ############### this fails
    # # }
    #
    #
    # if(length(stratum_weights) != max(stratum) & length(Nhat) != max(stratum)) {
    #   stop("Length of stratum weights does not equal max stratum value")
    # }

    # -- Stratified without error in Nhats --
    if(is.null(se_Nhat)) {
      # if(verbose) cat("\n", "without error in Nhat", "\n")
      if(is.na(FPC)) {
        FPC <- !is.null(Nhat) #& is.null(stratum_weights)
      }
      if(is.null(Nhat) & !is.null(stratum_weights)) {
        Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
        FPC <- FALSE
      } else {
        Nt <- Nhat # consistency with report eqns
      }

      # pulling these equations from the Jim Creek report
      if(length(unique(cats)) > 1) {

        # --- if there are proportions to estimate ---

        cat("\n","Estimate proportions, add FPC if present","\n")

        # ntz <- table(stratum, cats, useNA="no")
        # nt <- rowSums(ntz)
        # if(!FPC) {
        #   FPC_vec <- 1#rep(1, length(nt))  # this will ignore FPC in further calculations
        #   if(verbose) cat("\n","Finite Population Correction factor NOT used","\n")
        # } else {
        #   FPC_vec <- (Nt-nt)/(Nt-1)
        #   if(verbose) cat("\n","Finite Population Correction factor USED","\n")
        # }
        # ptz <- ntz/nt
        # vptz <- ptz*(1-ptz)/(nt-1)*FPC_vec
        # Ntz <- Nt*ptz
        # vNtz <- (Nt^2)*vptz
        # Nz <- colSums(Ntz)
        # vNz <- colSums(vNtz)
        # pz <- Nz/sum(Nt)
        # vpz <- vNz/(sum(Nt)^2)
        # out$phat <- pz
        # out$se_phat <- sqrt(vpz)
        if(!is.null(Nhat)) {

          cat("\n","Estimate abundance, add FPC if present","\n")


        #   out$Nhat <- Nz
        #   out$se_Nhat <- sqrt(vNz)
        }

        # summarizing length
        if(!is.null(length)) {

          cat("\n","Estimate mean length, add FPC if present","\n")

          # xbartz <- tapply(length, list(stratum, cats), mean, na.rm=TRUE)   ### this should have a na.rm=TRUE !!!
          # vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          # if(FPC) {
          #   FPC_vec2 <- (Ntz-vxbartz_denom)/(Ntz-1)  ####### new FPC here
          # } else {
          #   FPC_vec2 <- 1#rep(1, length(nt))
          # }
          # vxbartz <- tapply(length, list(stratum, cats), var, na.rm=TRUE)/vxbartz_denom * FPC_vec2
          #
          # Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          # xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=T)
          #
          # # Taylor series expansion (Mood et al. 1974)
          # prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          # vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
          #   colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)
          #
          # # ## estimate xbartz using my ugly derivation
          # # vsumNtzx <- colSums((Ntz^2)*vxbartz + (xbartz^2)*vNtz - vNtz*vxbartz)
          # # #vNz is calculated above
          # # covthing <- colSums(xbartz*vNtz)
          # # vxbarz <- ((colSums(Ntz*xbartz)/Nz)^2)*(vsumNtzx/(colSums(Ntz*xbartz)^2) + vNz/(Nz^2) - 2*covthing/(colSums(Ntz*xbartz)*Nz))
          #
          # out$n_length <- tapply(!is.na(length), cats, sum)
          # out$mn_length <- xbarz
          # out$se_length <- sqrt(vxbarz)
          # out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          # out$max_length <- tapply(length, cats, max, na.rm=TRUE)
        }
      } else {
        # --- if there are NO proportions to estimate ---

        cat("\n","Estimate mean length, add FPC if present","\n")

        # nt <- table(stratum, is.na(length))[,1] # table(stratum)   ### now excludes NA values in length
        # if(!FPC) {
        #   if(verbose) cat("\n","Finite Population Correction factor NOT used", "\n")
        #   FPC_vec <- 1#rep(1, length(nt))  # this will ignore FPC in further calculations
        # } else {
        #   if(verbose) cat("\n","Finite Population Correction factor USED", "\n")
        #   FPC_vec <- (Nt-nt)/(Nt-1)
        # }
        # if(!is.null(stratum_weights)) {
        #   Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
        # } else {
        #   Nt <- Nhat # consistency with report eqns
        # }
        # xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        # vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt
        # out$n_length <- sum(!is.na(length))
        # out$mn_length <- sum(Nt*xbart/sum(Nt))
        # out$se_length <- sqrt(sum(((Nt/sum(Nt))^2)*vxbart*FPC_vec))
        # out$min_length <- min(length, na.rm=TRUE)
        # out$max_length <- max(length, na.rm=TRUE)
      }
    } else {
      # -- Stratified with error in Nhats --
      # estimating proportions..
      # estimating abundance
      # summarizing length
      # if(verbose) cat("\n", "WITH error in Nhat", "\n")
      if(is.na(FPC)) {
        FPC <- FALSE
      }
      Nt <- Nhat # consistency with report eqns
      # modifying these equations from the Jim Creek report

      if(length(unique(cats)) > 1) {
        # --- if there are proportions to estimate ---

        cat("\n","Estimate proportions, add FPC if present","\n")

        # ntz <- table(stratum, cats, useNA="no")
        # nt <- rowSums(ntz)  # this will be robust to the presence of NA
        # if(is.null(Nhat) & !is.null(stratum_weights)) {
        #   stop("Nhat must be supplied when se_Nhat is used.")
        # } else {
        #   Nt <- Nhat # consistency with report eqns
        #   # FPC_vec <- (Nt-nt)/(Nt-1)
        # }
        # if(FPC) {
        #   FPC_vec <- (Nt-nt)/(Nt-1)
        #   if(verbose) cat("\n", "Finite Population Correction factor USED", "\n")
        # } else {
        #   FPC_vec <- 1#rep(1, length(Nt))
        #   if(verbose) cat("\n", "Finite Population Correction factor NOT used", "\n")
        # }
        # ptz <- ntz/nt
        # vptz <- ptz*(1-ptz)/(nt-1)*FPC_vec
        # Ntz <- Nt*ptz
        # vNtz <- ((Nt^2)*vptz) + ((ptz^2)*(se_Nhat^2)) - (vptz*(se_Nhat^2)) ### Goodman
        # Nz <- colSums(Ntz)
        # vNz <- colSums(vNtz)
        # pz <- Nz/sum(Nt)
        #
        # # delta method
        # # covNzSum <- vNz  ## this is a minimum, there should be more
        # covNzSum <- colSums((se_Nhat^2)*ptz, na.rm=TRUE)  # this seems unbiased!!
        #
        # vpz <- ((Nz/sum(Nz))^2)*(vNz/(Nz^2) + sum(se_Nhat^2)/((sum(Nz))^2) - 2*covNzSum/(Nz*sum(Nz)))
        # out$phat <- pz
        # out$se_phat <- sqrt(vpz)
        # # if(any(is.na(out$se_phat))) {
        # #   print(ntz)
        # #   print(vpz)
        # # }

        # out$Nhat <- Nz
        # out$se_Nhat <- sqrt(vNz)   ## this one seems unbiased still

        cat("\n","Estimate abundance, add FPC if present","\n")

        # summarizing length
        if(!is.null(length)) {
          # xbartz <- tapply(length, list(stratum, cats), mean, na.rm=TRUE)   ### this should have a na.rm=TRUE !!!
          # vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          # if(FPC) {
          #   FPC_vec2 <- (Ntz-vxbartz_denom)/(Ntz-1)  ####### new FPC here
          # } else {
          #   FPC_vec2 <- 1#rep(1, length(nt))
          # }
          # vxbartz <- tapply(length, list(stratum, cats), var, na.rm=TRUE)/vxbartz_denom * FPC_vec2
          #
          # Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          # xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=TRUE)
          # # Taylor series expansion (Mood et al. 1974)
          # prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          # vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
          #   colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)
          # out$n_length <- tapply(!is.na(length), cats, sum)
          # out$mn_length <- xbarz
          # out$se_length <- sqrt(vxbarz)
          # out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          # out$max_length <- tapply(length, cats, max, na.rm=TRUE)

          cat("\n","Estimate mean length, add FPC if present","\n")
        }
      } else {
        # --- if there are NO proportions to estimate ---

        cat("\n","Estimate mean length, add FPC if present","\n")
        # nt <- table(stratum, is.na(length))[,1] # table(stratum)   ### now excludes NA values in length
        # if(is.null(Nhat) & !is.null(stratum_weights)) {
        #   stop("Nhat must be supplied when se_Nhat is used.")
        # }
        # if(FPC) {
        #   FPC_vec <- (Nt-nt)/(Nt-1)
        #   if(verbose) cat("\n", "Finite Population Correction factor USED ------ except actually not", "\n")
        # } else {
        #   FPC_vec <- 1#rep(1, length(Nt))
        #   if(verbose) cat("\n", "Finite Population Correction factor NOT used", "\n")
        # }
        # Nt <- Nhat # consistency with report eqns
        # xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        # vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt
        # out$n_length <- sum(nt)
        # out$mn_length <- sum(Nt*xbart/sum(Nt))
        # # NoverSumN <- Nt/sum(Nt)
        # # # from delta method  -- validated by simulation   # currently no fpc anywhere
        # # vNoverSumN <- (NoverSumN^2)*(((se_Nhat^2)/(Nt^2)) +
        # #                                (sum(se_Nhat^2)/(sum(Nt)^2)) -
        # #                                (2*(se_Nhat^2)/(Nt*sum(Nt))))
        # # # from Goodman
        # # out$se_length <- sqrt(sum(((NoverSumN^2)*vxbart) +   # vxbart validated by simulation
        # #                             ((xbart^2)*vNoverSumN) -
        # #                             (vxbart*vNoverSumN))) #### there is huge covariance in the sum!!
        # # from covariance identities
        # covSumNtSum <- sum(xbart*(se_Nhat^2))
        # # from Goodman
        # vSumNtXbart <- sum(((Nt^2)*vxbart) +
        #                      (xbart^2)*(se_Nhat^2) -
        #                      ((se_Nhat^2)*vxbart))
        # # from delta method
        # out$se_length <- sqrt(((sum(Nhat*xbart)/sum(Nhat))^2)*((vSumNtXbart/((sum(Nhat*xbart))^2)) +
        #                                                          (sum(se_Nhat^2)/((sum(Nhat))^2)) -
        #                                                          ((2*covSumNtSum)/(sum(Nhat)*sum(Nhat*xbart)))))
        # out$min_length <- min(length, na.rm=TRUE)
        # out$max_length <- max(length, na.rm=TRUE)
      }
    }
  }
  # --------- NOT Stratified---------
  if(is.null(stratum)) {
    # if(verbose) cat("\n", "not stratified", "\n")

    if(is.na(FPC)) {
      FPC <- !is.null(Nhat) & is.null(se_Nhat)
    }
    if(FPC & (is.null(Nhat) )) {  #WRONG  | !is.null(se_Nhat)
      FPC <- FALSE
    }
    # if(FPC) {
    #   FPC_prop <- (Nhat-sum(out$n))/(Nhat-1)
    #   FPC_mn <- (Nhat-sum(out$n))/(Nhat-1)
    #   if(verbose) cat("\n", "Finite Population Correction factor USED", "\n")
    # } else {
    #   FPC_prop <- 1
    #   FPC_mn <- 1
    #   if(verbose) cat("\n", "Finite Population Correction factor NOT used", "\n")
    # }
    # # FPC <- ifelse(!is.null(Nhat), (Nhat-sum(out$n))/(Nhat-1), 1)
    # # FPC_mn <- FPC   # FPC used for variance of mean length
    #
    # # checking that inputs make sense
    # if(length(Nhat) > 1 | length(se_Nhat) > 1 | length(stratum_weights) > 1) {
    #   stop("Stratum totals or weights are given, but no strata")
    # }

    # estimating proportions..
    if(length(unique(cats)) > 1) {

      cat("\n","Estimate proportions, add FPC if present","\n")

      # out_ntot <- sum(out$n)
      # out$phat <- out$n/out_ntot
      # out$se_phat <- sqrt(out$phat*(1-out$phat)/(out_ntot-1)*FPC_prop)     ### added fpc

      # estimating abundance..
      if(!is.null(Nhat)) {


        # out$Nhat <- Nhat * out$phat
        if(is.null(se_Nhat)) {

          cat("\n","Estimate abundance, add FPC if present","\n")

          # if(verbose) cat("\n", "without error in Nhat", "\n")
          # out$se_Nhat <- Nhat * out$se_phat   # verify this??
        } else {

          cat("\n","Estimate abundance, add FPC if present","\n")

          # # Goodman's formula
          # if(verbose) cat("\n", "with error in Nhat", "\n")
          # out$se_Nhat <- sqrt(((Nhat^2)*(out$se_phat^2)) +
          #                       ((out$phat^2)*(se_Nhat^2)) -
          #                       ((out$se_phat^2)*(se_Nhat^2)))
        }

        # # updating FPC for mean if there are categories AND est abundance
        # if(FPC) {
        #   FPC_mn <- (out$Nhat - out$n)/(out$Nhat - 1)
        # }
      }
    }

    # and summarizing lengths..
    if(!is.null(length)) {

      cat("\n","Estimate mean length, add FPC if present","\n")

      # n_notNA <- tapply(!is.na(length), cats, sum)
      # out$n_length <- n_notNA
      # out$mn_length <- tapply(length, cats, mean, na.rm=TRUE)
      # out_sd_length <- tapply(length, cats, sd, na.rm=TRUE) * sqrt(FPC_mn)
      # # sqrt((out$Nhat-n_notNA)/(out$Nhat-1)) ## NEW FPC    # sqrt(FPC)     ### added fpc
      # out$se_length <- out_sd_length/sqrt(n_notNA)
      # out$min_length <- tapply(length, cats, min, na.rm=TRUE)
      # out$max_length <- tapply(length, cats, max, na.rm=TRUE)
    }
  }
  # return(out)
}

# ASL_boilerplate(age=sim_data$data$age,
#           length=sim_data$data$length,
#           stratum=sim_data$data$stratum,
#           Nhat=sim_data$abundance$Nhat,
#           se_Nhat=sim_data$abundance$se_Nhat)
#
# # Stratified without error in Nhat
# ASL_boilerplate(age=sim_data$data$age,
#           length=sim_data$data$length,
#           stratum=sim_data$data$stratum,
#           Nhat=sim_data$abundance$Nhat)
#
# # Pooled (not stratified) with error in Nhat
# ASL_boilerplate(age=sim_data$data$age,
#           length=sim_data$data$length,
#           Nhat=100*1000,
#           se_Nhat=10*1000)
#
# # Pooled (not stratified) without error in Nhat
# ASL_boilerplate(age=sim_data$data$age,
#           length=sim_data$data$length,
#           Nhat=100*1000)
