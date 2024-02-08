
### tests
library(testthat)
library(dsftools)

# need to check all possible combinations

# ASL_table(age=sim_data$data$age,
#           length=sim_data$data$length,
#           stratum=sim_data$data$stratum,
#           Nhat=sim_data$abundance$Nhat,
#           se_Nhat=sim_data$abundance$se_Nhat, FPC=NA)

# stratified vs pooled
# lengthage vs length vs age
# N known vs estimated
# FPC TRUE vs FALSE vs NA
countup <- 0
thelist <- list() # thesums <- NA
phat_arr <- se_phat_arr <- Nhat_arr <- se_Nhat_arr <- xbar_arr <- se_xbar_arr <-
  array(dim=c(2,3,2,3),
  dimnames=list(c("stratified","pooled"),
               c("lengthage","length","age"),
               c("N known","N estimated"),
               c("FPC", "no FPC","FPC NA")))
for(i_strat in 1:2) {
  for(i_which in 1:3) {
    for(i_known in 1:2) {
      # print("=============================================")
      for(i_fpc in 1:3) {
        countup <- countup + 1
        if(i_which==1) {
          age <- sim_data$data$age
          length <- sim_data$data$length
        }
        if(i_which==2) {
          age <- NULL
          length <- sim_data$data$length
        }
        if(i_which==3) {
          age <- sim_data$data$age
          length <- NULL
        }

        if(i_strat==1) {
          stratum <- sim_data$data$stratum
        }
        if(i_strat==2) {
          stratum <- NULL
        }

        if(i_known==1 & i_strat==1) {
          Nhat <- sim_data$abundance$Nhat
          se_Nhat <- NULL
        }
        if(i_known==1 & i_strat==2) {
          Nhat <- sum(sim_data$abundance$Nhat)
          se_Nhat <- NULL
        }
        if(i_known==2 & i_strat==1) {
          Nhat <- sim_data$abundance$Nhat
          se_Nhat <- sim_data$abundance$se_Nhat
        }
        if(i_known==2 & i_strat==2) {
          Nhat <- sum(sim_data$abundance$Nhat)
          se_Nhat <- sqrt(sum(sim_data$abundance$se_Nhat^2))
        }

        if(i_fpc==1) FPC <- "always" #TRUE
        if(i_fpc==2) FPC <- "never" #FALSE
        if(i_fpc==3) FPC <- "ifknown" #NA

        # print(ASL_table(age=age,
        #                 length=length,
        #                 stratum=stratum,
        #                 Nhat=Nhat,
        #                 se_Nhat=se_Nhat,
        #                 FPC=FPC,
        #                 verbose=TRUE))

        # thesums[countup] <- sum(ASL_table(age=age,
        #                                   length=length,
        #                                   stratum=stratum,
        #                                   Nhat=Nhat,
        #                                   se_Nhat=se_Nhat,
        #                                   FPC=FPC,
        #                                   verbose=FALSE))

        thelist[[countup]] <- ASL_table(age=age,
                                          length=length,
                                          stratum=stratum,
                                          Nhat=Nhat,
                                          se_Nhat=se_Nhat,
                                          FPC=FPC,
                                          verbose=FALSE)
        # if(abs(thesums[countup]-sums_check[countup])>0.0001) {
        #   print(c(i_strat,i_which,i_known,i_fpc))
        #   for(ii in 1:3) {
        #     print(ASL_table(age=age,
        #                     length=length,
        #                     stratum=stratum,
        #                     Nhat=Nhat,
        #                     se_Nhat=se_Nhat,
        #                     FPC=c("always","never","ifknown")[ii],
        #                     verbose=TRUE))
        #     print(thesums[countup]-sums_check[countup])
        #   }
        # }


        # phat_mn <- suppressWarnings(mean(thelist[[countup]]$phat))
        phat_mn <- suppressWarnings(thelist[[countup]]$phat[1])
        if(!is.null(phat_mn)) {
          phat_arr[i_strat,i_which,i_known,i_fpc] <- phat_mn
        }
        # se_phat_mn <- suppressWarnings(mean(thelist[[countup]]$se_phat))
        se_phat_mn <- suppressWarnings(thelist[[countup]]$se_phat[1])
        if(!is.null(se_phat_mn)) {
          se_phat_arr[i_strat,i_which,i_known,i_fpc] <- se_phat_mn
        }
        # Nhat_mn <- suppressWarnings(mean(thelist[[countup]]$Nhat))
        Nhat_mn <- suppressWarnings(thelist[[countup]]$Nhat[1])
        if(!is.null(Nhat_mn)) {
          Nhat_arr[i_strat,i_which,i_known,i_fpc] <- Nhat_mn
        }
        # se_Nhat_mn <- suppressWarnings(mean(thelist[[countup]]$se_Nhat))
        se_Nhat_mn <- suppressWarnings(thelist[[countup]]$se_Nhat[1])
        if(!is.null(se_Nhat_mn)) {
          se_Nhat_arr[i_strat,i_which,i_known,i_fpc] <- se_Nhat_mn
        }
        # xbar_mn <- suppressWarnings(mean(thelist[[countup]]$mn_length))
        xbar_mn <- suppressWarnings(thelist[[countup]]$mn_length[1])
        if(!is.null(xbar_mn)) {
          xbar_arr[i_strat,i_which,i_known,i_fpc] <- xbar_mn
        }
        # se_xbar_mn <- suppressWarnings(mean(thelist[[countup]]$se_mn_length))
        se_xbar_mn <- suppressWarnings(thelist[[countup]]$se_length[1])
        if(!is.null(se_xbar_mn)) {
          se_xbar_arr[i_strat,i_which,i_known,i_fpc] <- se_xbar_mn
        }
      }
    }
  }
}
# phat_arr[,1,,]
# phat_arr[,2,,]
# phat_arr[,3,,]
# se_phat_arr[,1,,]
# se_phat_arr[,2,,]
# se_phat_arr[,3,,]
#
# Nhat_arr[,1,,]
# Nhat_arr[,2,,]
# Nhat_arr[,3,,]
# se_Nhat_arr[,1,,]
# se_Nhat_arr[,2,,]
# se_Nhat_arr[,3,,]
#
# xbar_arr[,1,,]
# xbar_arr[,2,,]
# xbar_arr[,3,,]
# se_xbar_arr[,1,,]
# se_xbar_arr[,2,,]
# se_xbar_arr[,3,,]

# all(phat_arr[,1,1,] == phat_arr[,1,2,])  # N known same as N estimated
# all(phat_arr[,3,1,] == phat_arr[,3,2,])  # N known same as N estimated
# all(is.na(phat_arr[,2,,]))  # not estimated when no proportions
# all(phat_arr[,1,,1]==phat_arr[,1,,2] & phat_arr[,1,,3]==phat_arr[,1,,2]) # not affected by FPC
# all(phat_arr[,3,,1]==phat_arr[,3,,2] & phat_arr[,3,,3]==phat_arr[,3,,2]) # not affected by FPC
#
# all(se_phat_arr[1,1,1,] < se_phat_arr[1,1,2,])  # SE less when N known
# all(se_phat_arr[1,3,1,] < se_phat_arr[1,3,2,])
# all(se_phat_arr[2,1,1,1:2] == se_phat_arr[2,1,2,1:2])  # estimated N does not affect
# all(se_phat_arr[2,3,1,1:2] == se_phat_arr[2,3,2,1:2])
# all(se_phat_arr[1,1,,] > se_phat_arr[2,1,,])  # SE more when stratified
# all(se_phat_arr[1,3,,] > se_phat_arr[2,3,,])
# all(is.na(se_phat_arr[,2,,]))  # not estimated when no proportions
# all(se_phat_arr[,1,,1] < se_phat_arr[,1,,2])  # se less with fpc
# all(se_phat_arr[,3,,1] < se_phat_arr[,3,,2])  # se less with fpc
#
# all(Nhat_arr[,1,1,] == Nhat_arr[,1,2,])  # N known same as N estimated
# all(Nhat_arr[,3,1,] == Nhat_arr[,3,2,])  # N known same as N estimated
# all(is.na(Nhat_arr[,2,,]))  # not estimated when no proportions
# all(Nhat_arr[,1,,1]==Nhat_arr[,1,,2] & Nhat_arr[,1,,3]==Nhat_arr[,1,,2]) # not affected by FPC
# all(Nhat_arr[,3,,1]==Nhat_arr[,3,,2] & Nhat_arr[,3,,3]==Nhat_arr[,3,,2]) # not affected by FPC
#
# all(se_Nhat_arr[1,1,1,] < se_Nhat_arr[1,1,2,])  # SE less when N known
# all(se_Nhat_arr[1,3,1,] < se_Nhat_arr[1,3,2,])
# all(se_Nhat_arr[2,1,1,1:2] < se_Nhat_arr[2,1,2,1:2])  # less when N is known
# all(se_Nhat_arr[2,3,1,1:2] < se_Nhat_arr[2,3,2,1:2])
# all(se_Nhat_arr[1,1,,] > se_Nhat_arr[2,1,,])  # SE more when stratified
# all(se_Nhat_arr[1,3,,] > se_Nhat_arr[2,3,,])
# all(is.na(se_Nhat_arr[,2,,]))  # not estimated when no proportions
# all(se_Nhat_arr[,1,,1] < se_Nhat_arr[,1,,2])  # se less with fpc
# all(se_Nhat_arr[,3,,1] < se_Nhat_arr[,3,,2])  # se less with fpc
#
# all(xbar_arr[,1,1,] == xbar_arr[,1,2,])  # N known same as N estimated
# all(xbar_arr[,2,1,] == xbar_arr[,2,2,])  # N known same as N estimated
# all(is.na(xbar_arr[,3,,]))  # not estimated when no lengths
# all(xbar_arr[,1,,1]==xbar_arr[,1,,2] & xbar_arr[,1,,3]==xbar_arr[,1,,2]) # not affected by FPC
# all(xbar_arr[,2,,1]==xbar_arr[,2,,2] & xbar_arr[,2,,3]==xbar_arr[,2,,2]) # not affected by FPC
#
# all(se_xbar_arr[1,1,1,] < se_xbar_arr[1,1,2,])  # SE less when N known
# all(se_xbar_arr[1,2,1,] < se_xbar_arr[1,2,2,])
# all(se_xbar_arr[2,1,1,1:2] == se_xbar_arr[2,1,2,1:2])  # estimated N does not affect
# all(se_xbar_arr[2,2,1,1:2] == se_xbar_arr[2,2,2,1:2])
# all(se_xbar_arr[1,1,,] > se_xbar_arr[2,1,,])  # SE more when stratified
# all(se_xbar_arr[1,2,,] > se_xbar_arr[2,2,,])
# all(is.na(se_xbar_arr[,3,,]))  # not estimated when no proportions
# all(se_xbar_arr[,1,,1] < se_xbar_arr[,1,,2])  # se less with fpc
# all(se_xbar_arr[,2,1,1] < se_xbar_arr[,2,1,2])  # se less with fpc (subset)


# stratified with weights vs pooled and no Nhat given
# lengthage vs length vs age
# FPC TRUE vs FALSE vs NA
phat_arr1 <- se_phat_arr1 <- Nhat_arr1 <- se_Nhat_arr1 <- xbar_arr1 <- se_xbar_arr1 <-
  array(dim=c(2,3,3),
        dimnames=list(c("stratWts","pooled_noNhat"),
                      c("lengthage","length","age"),
                      c("FPC", "no FPC","FPC NA")))
for(i_strat in 1:2) {
  for(i_which in 1:3) {
    # print("=============================================")
    for(i_fpc in 1:3) {
      countup <- countup + 1
      if(i_which==1) {
        age <- sim_data$data$age
        length <- sim_data$data$length
      }
      if(i_which==2) {
        age <- NULL
        length <- sim_data$data$length
      }
      if(i_which==3) {
        age <- sim_data$data$age
        length <- NULL
      }

      if(i_strat==1) {
        stratum <- sim_data$data$stratum
        stratum_weights <- sim_data$abundance$Nhat/sum(sim_data$abundance$Nhat)
        Nhat <- NULL
      }
      if(i_strat==2) {
        stratum <- NULL
        Nhat <- NULL
        stratum_weights <- NULL
      }

      if(i_fpc==1) FPC <- "always" #TRUE
      if(i_fpc==2) FPC <- "never" #FALSE
      if(i_fpc==3) FPC <- "ifknown" #NA

      # print(ASL_table(age=age,
      #                 length=length,
      #                 stratum=stratum,
      # stratum_weights = stratum_weights,
      #                 Nhat=Nhat,
      #                 FPC=FPC,
      #                 verbose=TRUE))

      # thesums[countup] <- sum(ASL_table(age=age,
      #                                   length=length,
      #                                   stratum=stratum,
      #                                   stratum_weights = stratum_weights,
      #                                   Nhat=Nhat,
      #                                   FPC=FPC,
      #                                   verbose=FALSE))

      thelist[[countup]] <- ASL_table(age=age,
                                        length=length,
                                        stratum=stratum,
                                        stratum_weights = stratum_weights,
                                        Nhat=Nhat,
                                        FPC=FPC,
                                        verbose=FALSE)

      phat_mn <- suppressWarnings(thelist[[countup]]$phat[1])
      if(!is.null(phat_mn)) {
        phat_arr1[i_strat,i_which,i_fpc] <- phat_mn
      }
      # se_phat_mn <- suppressWarnings(mean(thelist[[countup]]$se_phat))
      se_phat_mn <- suppressWarnings(thelist[[countup]]$se_phat[1])
      if(!is.null(se_phat_mn)) {
        se_phat_arr1[i_strat,i_which,i_fpc] <- se_phat_mn
      }
      # Nhat_mn <- suppressWarnings(mean(thelist[[countup]]$Nhat))
      Nhat_mn <- suppressWarnings(thelist[[countup]]$Nhat[1])
      if(!is.null(Nhat_mn)) {
        Nhat_arr1[i_strat,i_which,i_fpc] <- Nhat_mn
      }
      # se_Nhat_mn <- suppressWarnings(mean(thelist[[countup]]$se_Nhat))
      se_Nhat_mn <- suppressWarnings(thelist[[countup]]$se_Nhat[1])
      if(!is.null(se_Nhat_mn)) {
        se_Nhat_arr1[i_strat,i_which,i_fpc] <- se_Nhat_mn
      }
      # xbar_mn <- suppressWarnings(mean(thelist[[countup]]$mn_length))
      xbar_mn <- suppressWarnings(thelist[[countup]]$mn_length[1])
      if(!is.null(xbar_mn)) {
        xbar_arr1[i_strat,i_which,i_fpc] <- xbar_mn
      }
      # se_xbar_mn <- suppressWarnings(mean(thelist[[countup]]$se_mn_length))
      se_xbar_mn <- suppressWarnings(thelist[[countup]]$se_length[1])
      if(!is.null(se_xbar_mn)) {
        se_xbar_arr1[i_strat,i_which,i_fpc] <- se_xbar_mn
      }
    }
  }
}

# # phat_arr1[,1,]
# # phat_arr1[,2,]
# # phat_arr1[,3,]
# # se_phat_arr1[,1,]
# # se_phat_arr1[,2,]
# # se_phat_arr1[,3,]
# #
# # Nhat_arr1[,1,]
# # Nhat_arr1[,2,]
# # Nhat_arr1[,3,]
# # se_Nhat_arr1[,1,]
# # se_Nhat_arr1[,2,]
# # se_Nhat_arr1[,3,]
# #
# # xbar_arr1[,1,]
# # xbar_arr1[,2,]
# # xbar_arr1[,3,]
# # se_xbar_arr1[,1,]
# # se_xbar_arr1[,2,]
# # se_xbar_arr1[,3,]
#
# # all(phat_arr1[,1,1,] == phat_arr1[,1,2,])  # N known same as N estimated
# # all(phat_arr1[,3,1,] == phat_arr1[,3,2,])  # N known same as N estimated
# all(is.na(phat_arr1[,2,]))  # not estimated when no proportions
# all(phat_arr1[,1,1]==phat_arr1[,1,2] & phat_arr1[,1,3]==phat_arr1[,1,2]) # not affected by FPC
# all(phat_arr1[,3,1]==phat_arr1[,3,2] & phat_arr1[,3,3]==phat_arr1[,3,2]) # not affected by FPC
#
# # all(se_phat_arr1[1,1,1,] < se_phat_arr1[1,1,2,])  # SE less when N known
# # all(se_phat_arr1[1,3,1,] < se_phat_arr1[1,3,2,])
# # all(se_phat_arr1[2,1,1,1:2] == se_phat_arr1[2,1,2,1:2])  # estimated N does not affect
# # all(se_phat_arr1[2,3,1,1:2] == se_phat_arr1[2,3,2,1:2])
# all(se_phat_arr1[1,1,] > se_phat_arr1[2,1,])  # SE more when stratified
# all(se_phat_arr1[1,3,] > se_phat_arr1[2,3,])
# all(is.na(se_phat_arr1[,2,]))  # not estimated when no proportions
# # all(se_phat_arr1[,1,1] < se_phat_arr1[,1,2])  # se less with fpc
# # all(se_phat_arr1[,3,1] < se_phat_arr1[,3,2])  # se less with fpc
#
# # all(Nhat_arr1[,1,1,] == Nhat_arr1[,1,2,])  # N known same as N estimated
# # all(Nhat_arr1[,3,1,] == Nhat_arr1[,3,2,])  # N known same as N estimated
# # all(is.na(Nhat_arr1[,2,,]))  # not estimated when no proportions
# # all(Nhat_arr1[,1,,1]==Nhat_arr1[,1,,2] & Nhat_arr1[,1,,3]==Nhat_arr1[,1,,2]) # not affected by FPC
# # all(Nhat_arr1[,3,,1]==Nhat_arr1[,3,,2] & Nhat_arr1[,3,,3]==Nhat_arr1[,3,,2]) # not affected by FPC
# #
# # all(se_Nhat_arr1[1,1,1,] < se_Nhat_arr1[1,1,2,])  # SE less when N known
# # all(se_Nhat_arr1[1,3,1,] < se_Nhat_arr1[1,3,2,])
# # all(se_Nhat_arr1[2,1,1,1:2] < se_Nhat_arr1[2,1,2,1:2])  # less when N is known
# # all(se_Nhat_arr1[2,3,1,1:2] < se_Nhat_arr1[2,3,2,1:2])
# # all(se_Nhat_arr1[1,1,,] > se_Nhat_arr1[2,1,,])  # SE more when stratified
# # all(se_Nhat_arr1[1,3,,] > se_Nhat_arr1[2,3,,])
# # all(is.na(se_Nhat_arr1[,2,,]))  # not estimated when no proportions
# # all(se_Nhat_arr1[,1,,1] < se_Nhat_arr1[,1,,2])  # se less with fpc
# # all(se_Nhat_arr1[,3,,1] < se_Nhat_arr1[,3,,2])  # se less with fpc
#
# # all(xbar_arr1[,1,1,] == xbar_arr1[,1,2,])  # N known same as N estimated
# # all(xbar_arr1[,2,1,] == xbar_arr1[,2,2,])  # N known same as N estimated
# all(is.na(xbar_arr1[,3,]))  # not estimated when no lengths
# all(xbar_arr1[,1,1]==xbar_arr1[,1,2] & xbar_arr1[,1,3]==xbar_arr1[,1,2]) # not affected by FPC
# all(xbar_arr1[,2,1]==xbar_arr1[,2,2] & xbar_arr1[,2,3]==xbar_arr1[,2,2]) # not affected by FPC
#
# # all(se_xbar_arr1[1,1,1,] < se_xbar_arr1[1,1,2,])  # SE less when N known
# # all(se_xbar_arr1[1,2,1,] < se_xbar_arr1[1,2,2,])
# # all(se_xbar_arr1[2,1,1,1:2] == se_xbar_arr1[2,1,2,1:2])  # estimated N does not affect
# # all(se_xbar_arr1[2,2,1,1:2] == se_xbar_arr1[2,2,2,1:2])
# all(se_xbar_arr1[1,1,] > se_xbar_arr1[2,1,])  # SE more when stratified
# all(se_xbar_arr1[1,2,] > se_xbar_arr1[2,2,])
# all(is.na(se_xbar_arr1[,3,]))  # not estimated when no proportions
# # all(se_xbar_arr1[,1,1] < se_xbar_arr1[,1,2])  # se less with fpc
# # all(se_xbar_arr1[,2,1,1] < se_xbar_arr1[,2,1,2])  # se less with fpc (subset)





therows <- sapply(thelist, nrow)
thecols <- sapply(thelist, ncol)
thesums <- sapply(thelist, sum)

rows_check <- c(5, 5, 5, 5, 5, 5, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                5, 5, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 1, 1, 5, 5,
                5, 5, 5, 5, 1, 1, 1, 5, 5, 5)
cols_check <- c(10, 10, 10, 10, 10, 10,  6,  6,  6,  6,  6,  6,  5,  5,  5,  5,
                5,  5, 10, 10, 10, 10, 10, 10, 6,  6,  6,  6,  6,  6,  5,  5,  5,
                5,  5, 5,  8,  8,  8,  6,  6,  6,  3,  3,  3,  8,  8,  8, 6,  6,
                6,  3,  3,  3)

sums_check <- c(120558.3029, 120575.2692, 120558.3712, 126112.8508, 126124.0163, 126124.0163,   1593.3601,
1593.3670,   1593.3601,   1593.8774,   1593.8774,   1593.8774, 116386.0509, 116402.9466,
116386.0509, 121940.2906, 121951.3860, 121951.3860, 120202.1631, 120221.7422, 120202.2332,
124932.9399, 124946.3919, 124946.3919,   1614.9671,   1614.9754,   1614.9671,   1614.9671,
1614.9754,   1614.9754, 116044.4955, 116064.0044, 116044.4955, 120775.2723, 120788.6541,
120788.6541,   4563.4247,   4563.4247,   4563.4247,   1593.3670,   1593.3670,   1593.3670,
391.1022,    391.1022,    391.1022,   4553.8367,   4553.8367,   4553.8367,   1614.9754,
1614.9754,   1614.9754,    396.0989,    396.0989,    396.0989)

# all(abs(thesums - thecheck) < 0.1)

test_that("ASL_table", {
  expect_true(all(phat_arr[,1,1,] == phat_arr[,1,2,]))  # N known same as N estimated
  expect_true(all(phat_arr[,3,1,] == phat_arr[,3,2,]))  # N known same as N estimated
  expect_true(all(is.na(phat_arr[,2,,])))  # not estimated when no proportions
  expect_true(all(phat_arr[,1,,1]==phat_arr[,1,,2] & phat_arr[,1,,3]==phat_arr[,1,,2])) # not affected by FPC
  expect_true(all(phat_arr[,3,,1]==phat_arr[,3,,2] & phat_arr[,3,,3]==phat_arr[,3,,2])) # not affected by FPC

  expect_true(all(se_phat_arr[1,1,1,] < se_phat_arr[1,1,2,]))  # SE less when N known
  expect_true(all(se_phat_arr[1,3,1,] < se_phat_arr[1,3,2,]))
  expect_true(all(se_phat_arr[2,1,1,1:2] == se_phat_arr[2,1,2,1:2]))  # estimated N does not affect
  expect_true(all(se_phat_arr[2,3,1,1:2] == se_phat_arr[2,3,2,1:2]))
  expect_true(all(se_phat_arr[1,1,,] > se_phat_arr[2,1,,]))  # SE more when stratified
  expect_true(all(se_phat_arr[1,3,,] > se_phat_arr[2,3,,]))
  expect_true(all(is.na(se_phat_arr[,2,,])))  # not estimated when no proportions
  expect_true(all(se_phat_arr[,1,,1] < se_phat_arr[,1,,2]))  # se less with fpc
  expect_true(all(se_phat_arr[,3,,1] < se_phat_arr[,3,,2]))  # se less with fpc

  expect_true(all(Nhat_arr[,1,1,] == Nhat_arr[,1,2,]))  # N known same as N estimated
  expect_true(all(Nhat_arr[,3,1,] == Nhat_arr[,3,2,]))  # N known same as N estimated
  expect_true(all(is.na(Nhat_arr[,2,,])))  # not estimated when no proportions
  expect_true(all(Nhat_arr[,1,,1]==Nhat_arr[,1,,2] & Nhat_arr[,1,,3]==Nhat_arr[,1,,2])) # not affected by FPC
  expect_true(all(Nhat_arr[,3,,1]==Nhat_arr[,3,,2] & Nhat_arr[,3,,3]==Nhat_arr[,3,,2])) # not affected by FPC

  expect_true(all(se_Nhat_arr[1,1,1,] < se_Nhat_arr[1,1,2,]))  # SE less when N known
  expect_true(all(se_Nhat_arr[1,3,1,] < se_Nhat_arr[1,3,2,]))
  expect_true(all(se_Nhat_arr[2,1,1,1:2] < se_Nhat_arr[2,1,2,1:2]))  # less when N is known
  expect_true(all(se_Nhat_arr[2,3,1,1:2] < se_Nhat_arr[2,3,2,1:2]))
  expect_true(all(se_Nhat_arr[1,1,,] > se_Nhat_arr[2,1,,]))  # SE more when stratified
  expect_true(all(se_Nhat_arr[1,3,,] > se_Nhat_arr[2,3,,]))
  expect_true(all(is.na(se_Nhat_arr[,2,,])))  # not estimated when no proportions
  expect_true(all(se_Nhat_arr[,1,,1] < se_Nhat_arr[,1,,2]))  # se less with fpc
  expect_true(all(se_Nhat_arr[,3,,1] < se_Nhat_arr[,3,,2]))  # se less with fpc

  expect_true(all(xbar_arr[,1,1,] == xbar_arr[,1,2,]))  # N known same as N estimated
  expect_true(all(xbar_arr[,2,1,] == xbar_arr[,2,2,]))  # N known same as N estimated
  expect_true(all(is.na(xbar_arr[,3,,])))  # not estimated when no lengths
  expect_true(all(xbar_arr[,1,,1]==xbar_arr[,1,,2] & xbar_arr[,1,,3]==xbar_arr[,1,,2])) # not affected by FPC
  expect_true(all(xbar_arr[,2,,1]==xbar_arr[,2,,2] & xbar_arr[,2,,3]==xbar_arr[,2,,2])) # not affected by FPC

  expect_true(all(se_xbar_arr[1,1,1,] < se_xbar_arr[1,1,2,]))  # SE less when N known
  expect_true(all(se_xbar_arr[1,2,1,] < se_xbar_arr[1,2,2,]))
  expect_true(all(se_xbar_arr[2,1,1,1:2] == se_xbar_arr[2,1,2,1:2]))  # estimated N does not affect
  expect_true(all(se_xbar_arr[2,2,1,1:2] == se_xbar_arr[2,2,2,1:2]))
  expect_true(all(se_xbar_arr[1,1,,] > se_xbar_arr[2,1,,]))  # SE more when stratified
  expect_true(all(se_xbar_arr[1,2,,] > se_xbar_arr[2,2,,]))
  expect_true(all(is.na(se_xbar_arr[,3,,])))  # not estimated when no proportions
  expect_true(all(se_xbar_arr[,1,,1] < se_xbar_arr[,1,,2]))  # se less with fpc
  expect_true(all(se_xbar_arr[,2,1,1] < se_xbar_arr[,2,1,2]))  # se less with fpc (subset)

  expect_true(all(is.na(phat_arr1[,2,])))  # not estimated when no proportions
  expect_true(all(phat_arr1[,1,1]==phat_arr1[,1,2] & phat_arr1[,1,3]==phat_arr1[,1,2])) # not affected by FPC
  expect_true(all(phat_arr1[,3,1]==phat_arr1[,3,2] & phat_arr1[,3,3]==phat_arr1[,3,2])) # not affected by FPC

  expect_true(all(se_phat_arr1[1,1,] > se_phat_arr1[2,1,]))  # SE more when stratified
  expect_true(all(se_phat_arr1[1,3,] > se_phat_arr1[2,3,]))
  expect_true(all(is.na(se_phat_arr1[,2,])))  # not estimated when no proportions

  expect_true(all(is.na(xbar_arr1[,3,])))  # not estimated when no lengths
  expect_true(all(xbar_arr1[,1,1]==xbar_arr1[,1,2] & xbar_arr1[,1,3]==xbar_arr1[,1,2])) # not affected by FPC
  expect_true(all(xbar_arr1[,2,1]==xbar_arr1[,2,2] & xbar_arr1[,2,3]==xbar_arr1[,2,2])) # not affected by FPC

  expect_true(all(se_xbar_arr1[1,1,] > se_xbar_arr1[2,1,]))  # SE more when stratified
  expect_true(all(se_xbar_arr1[1,2,] > se_xbar_arr1[2,2,]))
  expect_true(all(is.na(se_xbar_arr1[,3,])))  # not estimated when no proportions

  expect_true(all(therows == rows_check))
  expect_true(all(thecols == cols_check))
  expect_true(all(abs(thesums - sums_check) < 0.0001))
})

cases <- c("stratified_witherror_lengthage", "stratified_witherror_age",
           "stratified_witherror_length", "stratified_lengthage", "stratified_age",
           "stratified_length", "pooled_witherror_lengthage", "pooled_witherror_age",
           "pooled_witherror_length", "pooled_lengthage", "pooled_age", "pooled_length")
test_that("verify_ASL_table", {
  for(case_i in cases) {
    expect_silent(verify_ASL_table(case=case_i, nsim=10, verbose=FALSE))
  }
})

# for(case_i in cases) {
#   verify_ASL_table(case=case_i, nsim=10, verbose=FALSE)
# }
