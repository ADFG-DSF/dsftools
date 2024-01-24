#' Construct an ASL Table from Vectors of Data
#' @description FILL IN DESCRIPTION HERE
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
#' @return A data.frame with rows corresponding to categories of age and/or sex,
#' depending on data inputs.
#' @author Matt Tyers
#' @examples
#'
#' @export
ASL_table <- function(age=NULL,
                      sex=NULL,
                      length=NULL,
                      stratum=NULL,
                      Nhat=NULL,
                      se_Nhat=NULL,
                      stratum_weights = NULL,
                      verbose=FALSE) {
  # -------- globally dealing with inputs ----------
  # combining age/sex categories as available
  if(!is.null(sex) & !is.null(age)) {
    cats <- paste(sex, age)
  } else {
    if(!is.null(sex)) cats <- sex
    if(!is.null(age)) cats <- age
    if(is.null(sex) & is.null(age)) cats <- rep("Total", length(length))
  }
  # this will become output dataframe
  out <- data.frame(n=tapply(X=seq_along(cats), INDEX=cats, FUN=base::length))
  # --------- Stratified ----------
  if(!is.null(stratum)) {
    if(verbose) print("stratified")
    # checking that inputs make sense

    # if(!inherits(stratum, "integer")) {
    #   stop("Stratum values must be positive whole numbers (1, 2, 3, etc.)")  ############### this fails
    # }
    if(length(stratum_weights) != max(stratum) & length(Nhat) != max(stratum)) {
      stop("Length of stratum weights does not equal max stratum value")
    }

    # -- Stratified without error in Nhats --
    if(is.null(se_Nhat)) {
      if(verbose) print("without error in Nhat")
      # pulling these equations from the Jim Creek report
      if(length(unique(cats)) > 1) {

        # --- if there are proportions to estimate ---
        ntz <- table(stratum, cats, useNA="no")
        nt <- rowSums(ntz)
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
          FPC_vec <- rep(1, length(nt))  # this will ignore FPC in further calculations
        } else {
          Nt <- Nhat # consistency with report eqns
          FPC_vec <- (Nt-nt)/(Nt-1)
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
        out$Nhat <- Nz
        out$se_Nhat <- sqrt(vNz)

        # summarizing length
        if(!is.null(length)) {
          xbartz <- tapply(length, list(stratum, cats), mean)
          vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          vxbartz <- tapply(length, list(stratum, cats), var)/vxbartz_denom
          Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=T)
          # Taylor series expansion (Mood et al. 1974)
          prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
            colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)
          out$mn_length <- xbarz
          out$se_length <- sqrt(vxbarz)
          out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          out$max_length <- tapply(length, cats, max, na.rm=TRUE)
        }
      } else {
        # --- if there are NO proportions to estimate ---
        nt <- table(stratum, is.na(length))[,1] # table(stratum)   ### now excludes NA values in length
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          FPC_vec <- rep(1, length(nt))  # this will ignore FPC in further calculations
          Nt <- stratum_weights   # weights will be used instead of true abundance (ok because they are assumed proportional)
        } else {
          FPC_vec <- (Nt-nt)/(Nt-1)
          Nt <- Nhat # consistency with report eqns
        }
        xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt
        out$mn_length <- sum(Nt*xbart/sum(Nt))
        out$se_length <- sqrt(sum(((Nt/sum(Nt))^2)*vxbart*FPC_vec))
        out$min_length <- min(length, na.rm=TRUE)
        out$max_length <- max(length, na.rm=TRUE)
      }
    } else {
      # -- Stratified with error in Nhats --
      # estimating proportions..
      # estimating abundance
      # summarizing length
      if(verbose) print("WITH error in Nhat")
      # modifying these equations from the Jim Creek report

      if(length(unique(cats)) > 1) {
        # --- if there are proportions to estimate ---
        ntz <- table(stratum, cats, useNA="no")
        nt <- rowSums(ntz)
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          stop("what the heck are you doing, boss?")
        } else {
          Nt <- Nhat # consistency with report eqns
          FPC_vec <- (Nt-nt)/(Nt-1)
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
        if(any(is.na(out$se_phat))) {
          print(ntz)
          print(vpz)
        }

        out$Nhat <- Nz
        out$se_Nhat <- sqrt(vNz)   ## this one seems unbiased still

        # summarizing length
        if(!is.null(length)) {
          xbartz <- tapply(length, list(stratum, cats), mean)
          vxbartz_denom <- table(stratum, cats, is.na(length))[,,1]  # this is like table(stratum, cats, useNA="no") but also excludes NA in length
          vxbartz <- tapply(length, list(stratum, cats), var)/vxbartz_denom
          Nz_mat <- matrix(Nz, nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          xbarz <- colSums(xbartz*Ntz/Nz_mat, na.rm=T)
          # Taylor series expansion (Mood et al. 1974)
          prod_mat <- matrix(colSums(xbartz*Ntz, na.rm=T), nrow=nrow(Ntz), ncol=ncol(Ntz), byrow=T)
          vxbarz <- colSums((Ntz^2)/(Nz_mat^2)*vxbartz, na.rm=T) +
            colSums(((xbartz*Nz_mat - prod_mat)^2)/(Nz_mat^4)*vNtz, na.rm=T)
          out$mn_length <- xbarz
          out$se_length <- sqrt(vxbarz)
          out$min_length <- tapply(length, cats, min, na.rm=TRUE)
          out$max_length <- tapply(length, cats, max, na.rm=TRUE)
        }
      } else {
        # --- if there are NO proportions to estimate ---
        nt <- table(stratum, is.na(length))[,1] # table(stratum)   ### now excludes NA values in length
        if(is.null(Nhat) & !is.null(stratum_weights)) {
          stop("what the heck are you doing, boss?")
        } else {
          FPC_vec <- (Nt-nt)/(Nt-1)
          Nt <- Nhat # consistency with report eqns
        }
        xbart <- tapply(length, stratum, mean, na.rm=TRUE)
        vxbart <- tapply(length, stratum, var, na.rm=TRUE)/nt
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
    if(verbose) print("not stratified")
    FPC <- ifelse(!is.null(Nhat), (Nhat-sum(out$n))/(Nhat-1), 1)

    # checking that inputs make sense
    if(length(Nhat) > 1 | length(se_Nhat) > 1 | length(stratum_weights) > 1) {
      stop("Stratum totals or weights are given, but no strata")
    }

    # estimating proportions..
    if(length(unique(cats)) > 1) {
      out_ntot <- sum(out$n)
      out$phat <- out$n/out_ntot
      out$se_phat <- sqrt(out$phat*(1-out$phat)/(out_ntot-1)*FPC)     ### added fpc

      # estimating abundance..
      if(!is.null(Nhat)) {
        out$Nhat <- Nhat * out$phat
        if(is.null(se_Nhat)) {
          if(verbose) print("without error in Nhat")
          out$se_Nhat <- Nhat * out$se_phat   # verify this??
        } else {
          # Goodman's formula
          if(verbose) print("with error in Nhat")
          out$se_Nhat <- sqrt(((Nhat^2)*(out$se_phat^2)) +
                                ((out$phat^2)*(se_Nhat^2)) -
                                ((out$se_phat^2)*(se_Nhat^2)))
        }
      }
    }

    # and summarizing lengths..
    if(!is.null(length)) {
      out$mn_length <- tapply(length, cats, mean, na.rm=TRUE)
      out_sd_length <- tapply(length, cats, sd, na.rm=TRUE)*sqrt(FPC)     ### added fpc
      n_notNA <- tapply(!is.na(length), cats, sum)
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
#   # age=age,
#   # sex=sex,
#   length=length,
#   # Nhat=10000,
#   # se_Nhat=5000
#   stratum_weights=Nhat/10,
#   stratum=stratum
# )

# include levels for all possible combinations (including those that don't occur?)
# incorporate fpc??
# think through what (should) happen with NA
# sampling_weights instead of Nhat - no fpc ..... calculations work with weights IN Nhat vec (except FPC), it all scales
