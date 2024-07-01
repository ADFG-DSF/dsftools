
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
1593.3670,   1593.3601,   1593.8714,   1593.8774,   1593.8774, 116386.0509, 116402.9466,
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

## making sure verify_, simulate_, and rp_ all RUN with every pre-defined case
cases <- c("stratified_witherror_lengthage", "stratified_witherror_age",
           "stratified_witherror_length", "stratified_lengthage", "stratified_age",
           "stratified_length", "stratified_Nunknown_lengthage", "stratified_Nunknown_age",
           "stratified_Nunknown_length", "pooled_witherror_lengthage", "pooled_witherror_age",
           "pooled_witherror_length", "pooled_lengthage", "pooled_age", "pooled_length",
           "pooled_Nunknown_lengthage", "pooled_Nunknown_age", "pooled_Nunknown_length")
test_that("verify_ASL_table", {
  for(case_i in cases) {
    expect_silent(verify_ASL_table(case=case_i, nsim=10, verbose=FALSE))
  }
})
test_that("simulate_ASL_table", {
  for(case_i in cases) {
    expect_silent(simulate_ASL_table(case=case_i, nsim=10, verbose=FALSE))
  }
})
test_that("rp_ASL_table", {
  for(case_i in cases) {
    expect_silent(rp_ASL_table(case=case_i, nsim=10, verbose=FALSE))
  }
})

## testing behavior of simulate_ASL_table     500 sims in 12 sec
simslist <- list()
for(i in 1:length(cases)) {
  # print(i)
  simslist[[i]] <- simulate_ASL_table(case=cases[i], nsim=100, verbose=FALSE, plot_pop=FALSE)
}


simlengths <- c(12,  8,  4, 12,  8,  4,  8,  4,  4, 12,  8,  4, 12,  8,  4,  8,  4,  4)

simdims_phat <- sapply(simslist, function(x) dim(x$phat_sim))
simdims_Nhat <- sapply(simslist, function(x) dim(x$Nhat_sim))
simdims_mn_length <- sapply(simslist, function(x) dim(x$mn_length_sim))
simdims_se_phat <- sapply(simslist, function(x) dim(x$se_phat_sim))
simdims_se_Nhat <- sapply(simslist, function(x) dim(x$se_Nhat_sim))
simdims_se_mn_length <- sapply(simslist, function(x) dim(x$se_mn_length_sim))

truelengths_phat <- sapply(simslist, function(x) length(x$phat_true))
truelengths_Nhat <- sapply(simslist, function(x) length(x$Nhat_true))
truelengths_mn_length <- sapply(simslist, function(x) length(x$mn_length_true))
truelengths_se_phat <- sapply(simslist, function(x) length(x$se_phat_true))
truelengths_se_Nhat <- sapply(simslist, function(x) length(x$se_Nhat_true))
truelengths_se_mn_length <- sapply(simslist, function(x) length(x$se_mn_length_true))

# defining a custom function to make sure the true values are within
# the quartiles of the sims
check_quartiles <- function(simvals, truevals, lwr=0.1, upr=.9) {   ## lwr and upr used to be .25 and .75
  if(is.null(simvals) & is.null(truevals)) {
    out <- TRUE  # don't complain if the sim & true vals aren't there
  }
  if(xor(is.null(simvals), is.null(truevals))) {
    out <- FALSE  # complain if one is there and the other isn't
  }
  if(!is.null(simvals) & !is.null(truevals)) {
    out <- apply(simvals, 2, quantile, lwr) <= truevals &
      apply(simvals, 2, quantile, upr) >= truevals
  }
  return(out)
}

test_that("simulate_ASL_table intensive", {
  expect_equal(sapply(simslist, length), simlengths)

  expect_equal(sum(unlist(simdims_phat)), 1260)
  expect_equal(sum(unlist(simdims_Nhat)), 840)
  expect_equal(sum(unlist(simdims_mn_length)), 1236)
  expect_equal(sum(unlist(simdims_se_phat)), 1260)
  expect_equal(sum(unlist(simdims_se_Nhat)), 840)
  expect_equal(sum(unlist(simdims_se_mn_length)), 1236)

  expect_equal(sum(truelengths_phat), 60)
  expect_equal(sum(truelengths_Nhat), 40)
  expect_equal(sum(truelengths_mn_length), 36)
  expect_equal(sum(truelengths_se_phat), 60)
  expect_equal(sum(truelengths_se_Nhat), 40)
  expect_equal(sum(truelengths_se_mn_length), 36)

  for(i in 1:length(cases)) {
    expect_true(all(check_quartiles(simvals = simslist[[i]]$phat_sim,
                                    truevals = simslist[[i]]$phat_true)))
    expect_true(all(check_quartiles(simvals = simslist[[i]]$Nhat_sim,
                                    truevals = simslist[[i]]$Nhat_true)))
    expect_true(all(check_quartiles(simvals = simslist[[i]]$mn_length_sim,
                                    truevals = simslist[[i]]$mn_length_true)))
    ## these are less stable in pooled case (weird!!)  Need to figure out
    # expect_true(all(check_quartiles(simvals = simslist[[i]]$se_phat_sim,
    #                                 truevals = simslist[[i]]$se_phat_true)))
    # expect_true(all(check_quartiles(simvals = simslist[[i]]$se_Nhat_sim,
    #                                 truevals = simslist[[i]]$se_Nhat_true)))
    # expect_true(all(check_quartiles(simvals = simslist[[i]]$se_mn_length_sim,
    #                                 truevals = simslist[[i]]$se_mn_length_true)))
  }
})

# for(case_i in cases) {
#   verify_ASL_table(case=case_i, nsim=10, verbose=FALSE)
# }

test_that("logit, expit, se", {
  expect_equal(sum(logit(seq(0.5, 0.99, by=0.01))), 66.09069, tolerance = 0.00001)
  expect_equal(sum(expit(seq(0.5, 0.99, by=0.01))), 33.86414, tolerance = 0.00001)
  expect_equal(sum(is.na(logit(c(seq(from=0.01, to=0.99, by=0.01), NA)))), 1)
  expect_equal(sum(!is.na(logit(c(seq(from=0.01, to=0.99, by=0.01), NA)))), 99)
  expect_equal(sum(is.na(expit(c(seq(from=-5, to=5, by=0.1), NA)))), 1)
  expect_equal(sum(!is.na(expit(c(seq(from=-5, to=5, by=0.1), NA)))), 101)
  expect_equal(se(c(8, 6, 7, 5, 3, 0, 9, NA), na.rm=TRUE), 1.172241, tolerance = 0.00001)
  expect_true(is.na(se(c(8, 6, 7, 5, 3, 0, 9, NA))))
})

xx <- rnorm(n=100000, mean=10, sd=1)
check1 <- c(0.1276836, 0.1646694, 0.1961813, 0.2573959)
check2 <- c(0.68156, 0.95510, 0.99731)
check3 <- c(1.277260, 1.640136, 1.954461, 2.568028)
check4 <- c(0.07989, 0.15827, 0.23603)
names(check1) <- names(check3) <- c('80%', '90%', '95%', '99%')
names(check2) <- names(check4) <- c("0.1", "0.2", "0.3")
test_that("rp", {
  expect_equal(rp(sim_vec = xx, true_val = 10, confidence = c(0.8, 0.9, 0.95, 0.99)),
               check1,
               tolerance = 0.01)
  expect_equal(rp(sim_vec = xx, true_val = 10, accuracy = c(0.1, 0.2, 0.3)),
               check2,
               tolerance=0.01)
  expect_equal(rp(sim_vec = xx, true_val = 10, confidence = c(0.8, 0.9, 0.95, 0.99), relative=FALSE),
               check3,
               tolerance = 0.01)
  expect_equal(rp(sim_vec = xx, true_val = 10, accuracy = c(0.1, 0.2, 0.3), relative=FALSE),
               check4,
               tolerance=0.01)
})

n <- 100
xx <- rnorm(n, mean=10, sd=1)
dfsim <- data.frame(a = xx + rnorm(n, sd=0.5),
                    b = xx + rnorm(n, sd=0.5),
                    c = xx + rnorm(n, sd=2),
                    d = xx + rnorm(n, sd=4),
                    e = -xx + rnorm(n, sd=2),
                    f = -xx + rnorm(n, sd=0.5),
                    g = -xx + rnorm(n, sd=0.5))
dfsim2 <- dfsim3 <- dfsim
names(dfsim2)[1:4] <- paste0("aa[", 1:4, "]")
dfsim3[5, 5] <- NA

test_that("plotcor", {
  expect_silent(plotcor(cor(dfsim)))
  expect_silent(plotcor(cor(dfsim), mincor=.5, legend=FALSE))
  expect_silent(plotcor(cor(dfsim), colors=c("cornflowerblue","orange")))
  expect_silent(plotcor(cor(dfsim2), maxn=3))
  expect_silent(plotcor(cor(dfsim2)))
  expect_silent(plotcor(cor(dfsim3)))
  expect_silent(plotcor(cor(dfsim3, use="na.or.complete")))
})

xx1 <- rnorm(100)
xx2 <- c(xx1, NA)
test_that("inside", {
  expect_equal(xx1 %inside% 1:3, (xx1 >= 1 & xx1 <= 3))
  expect_equal(xx2 %inside% 1:3, ifelse(is.na(xx2), FALSE, (xx2 >= 1 & xx2 <= 3)))
  expect_equal(which(1:10 %inside% 1:3), 1:3)
  expect_equal(which(1:10 %inside()% 1:3), 2)
  expect_equal(which(1:10 %inside[)% 1:3), 1:2)
  expect_equal(which(1:10 %inside(]% 1:3), 2:3)
})

test_that("subsets", {
  expect_equal(xx1 %s_l% 0, xx1[xx1 < 0])
  expect_equal(xx1 %s_leq% 0, xx1[xx1 <= 0])
  expect_equal(xx1 %s_g% 0, xx1[xx1 > 0])
  expect_equal(xx1 %s_geq% 0, xx1[xx1 >= 0])
  expect_equal(xx2 %s_l% 0, xx2[xx2 < 0])
  expect_equal(1:10 %s_inside% 1:3, 1:3)
  expect_equal(1:10 %s_inside()% 1:3, 2)
  expect_equal(1:10 %s_inside[)% 1:3, 1:2)
  expect_equal(1:10 %s_inside(]% 1:3, 2:3)
})
