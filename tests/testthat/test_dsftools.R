
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
# error vs not
# lengthage vs length vs age
# N known vs estimated
# FPC TRUE vs FALSE vs NA
countup <- 0
thesums <- NA
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

        if(i_fpc==1) FPC <- TRUE
        if(i_fpc==2) FPC <- FALSE
        if(i_fpc==3) FPC <- NA

        # print(ASL_table(age=age,
        #                 length=length,
        #                 stratum=stratum,
        #                 Nhat=Nhat,
        #                 se_Nhat=se_Nhat,
        #                 FPC=FPC,
        #                 verbose=TRUE))

        thesums[countup] <- sum(ASL_table(age=age,
                        length=length,
                        stratum=stratum,
                        Nhat=Nhat,
                        se_Nhat=se_Nhat,
                        FPC=FPC,
                        verbose=FALSE))
        # if(abs(thesums[countup]-thecheck[countup])>0.01) {
        #   print(c(i_strat,i_which,i_known,i_fpc))
        # }
      }
    }
  }
}

# stratified with weights vs pooled and no Nhat given
# lengthage vs length vs age
# stratum_weights, no Nhat given
# FPC TRUE vs FALSE vs NA
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

      if(i_fpc==1) FPC <- TRUE
      if(i_fpc==2) FPC <- FALSE
      if(i_fpc==3) FPC <- NA

      # print(ASL_table(age=age,
      #                 length=length,
      #                 stratum=stratum,
      # stratum_weights = stratum_weights,
      #                 Nhat=Nhat,
      #                 FPC=FPC,
      #                 verbose=TRUE))

      thesums[countup] <- sum(ASL_table(age=age,
                                        length=length,
                                        stratum=stratum,
                                        stratum_weights = stratum_weights,
                                        Nhat=Nhat,
                                        FPC=FPC,
                                        verbose=FALSE))
    }
  }
}



thecheck <- c(120196.2650, 120213.2479, 120196.2650, 125756.9079, 125768.0804, 125768.0804,
1203.3601,   1203.3670,   1203.3601,   1593.8774,   1593.8774,   1593.8774,
116386.0509, 116402.9466, 116386.0509, 121940.2906, 121951.3860, 121951.3860,
120202.1631, 120221.7422, 120202.1631, 124946.3919, 124946.3919, 124946.3919,
1614.9671,   1614.9754,   1614.9671,   1614.9754,   1614.9754,   1614.9754,
116044.4955, 116064.0044, 116044.4955, 120788.6541, 120788.6541, 120788.6541,
4201.4035,   4201.4035,   4201.4035,   1203.3670,   1203.3670,   1203.3670,
391.1022,    391.1022,    391.1022,   4553.8367,   4553.8367,   4553.8367,
1614.9754,   1614.9754,   1614.9754,    396.0989,    396.0989,    396.0989)

# all(abs(thesums - thecheck) < 0.1)

test_that("ASL_table", {
  expect_true(all(abs(thesums - thecheck) < 0.001))
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
