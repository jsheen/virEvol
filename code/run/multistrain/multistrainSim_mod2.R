rm(list = ls())
# ------------------------------------------------------------------------------
# Model 2: SEIR with vaccination 
# ------------------------------------------------------------------------------
# Source functional scripts and set seed and nsim ------------------------------
source('~/virEvol/code/func/gen_parameters.R')
source('~/virEvol/code/func/multistrain_eqns.R')
set.seed(0)
nsim <- 100
interyear_input <- 10
maxyear_input <- 100

# Set model 2 specific parameters and functions --------------------------------
c1=0.1
c2=0.3
# Initial susceptible population
fS_init = pop_size - 1
# Initial strain 1 infectious population in farms
fI1_init = 1
# Percent of susceptible chickens vaccinated in each time period
perc_vax = 0.33 
# Time that perc_vax is vaccinated
inter_vax_time = 120
# Vaccination rate of chickens of farms per susceptible chicken of farm per day
v = perc_vax / inter_vax_time 
# Rate of loss of immunity due to vaccination per chicken per day
v_hat = (1 / 120)
# Threshold value for extinction
threshold_extinction = 1.1
# Assign model 2
eqn <- multi_eqn_mod2

# Create function that will output the virulence strategies every interyear ----
multistrainSim_mod2 <- function(interyear=interyear_input, maxyear=maxyear_input) {
  # Check that maxyear divisible by interyear
  if (maxyear %% interyear != 0) {
    stop('Maxyear must be divisible by interyear in order for while loop to work.')
  }
  
  # Randomly select first strain
  res_vir <- runif(1, min=2, max=100)
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  p_1 <- (res_vir) / 100
  
  # Parameters
  parameters <- c(fbet1=fbet1, fbet2=0, fbet3=0, fbet4=0, fbet5=0, fbet6=0, fbet7=0, fbet8=0, fbet9=0, fbet10=0,
                  p_1=p_1, p_2=0, p_3=0, p_4=0, p_5=0, p_6=0, p_7=0, p_8=0, p_9=0, p_10=0,
                  sig=sig, gamm=gamm,
                  mort=mort, b=b, nat_mort=nat_mort,
                  v=v, v_hat=v_hat)
  
  # Prepare result dataframe
  virstrats <- data.frame(matrix(NA, ncol=10, nrow=(maxyear / interyear) + 1))
  virstrats[1,1] <- res_vir
  
  # Run epidemic
  endyear <- interyear
  init <- c(fS=fS_init,
            fE1=0, fE2=0, fE3=0, fE4=0, fE5=0, fE6=0, fE7=0, fE8=0, fE9=0, fE10=0,
            fI1=fI1_init, fI2=0, fI3=0, fI4=0, fI5=0, fI6=0, fI7=0, fI8=0, fI9=0, fI10=0,
            fR=0, fV=0, 
            fVE1=0, fVE2=0, fVE3=0, fVE4=0, fVE5=0, fVE6=0, fVE7=0, fVE8=0, fVE9=0, fVE10=0,
            fVI1=0, fVI2=0, fVI3=0, fVI4=0, fVI5=0, fVI6=0, fVI7=0, fVI8=0, fVI9=0, fVI10=0)
  time <- seq(0, interyear*365, by = 0.1)
  out <- ode(y=init, times=time, func=eqn, parms=parameters)
  out.df <- as.data.frame(out)
  while (endyear < maxyear) {
    print(endyear)
    # Set to 0 for those that are extinct
    for (col_dex in 3:12) {
      if ((out.df[nrow(out.df), col_dex] + out.df[nrow(out.df), (col_dex + 10)] + out.df[nrow(out.df), (col_dex + 22)] + out.df[nrow(out.df), (col_dex + 32)]) < threshold_extinction) {
        out.df[nrow(out.df), col_dex] <- 0
        out.df[nrow(out.df), (col_dex + 10)] <- 0
        out.df[nrow(out.df), (col_dex + 22)] <- 0
        out.df[nrow(out.df), (col_dex + 32)] <- 0
        
        # Set parameters for this strain to 0 as well
        parameters[col_dex - 2] <- 0
        parameters[col_dex - 2 + 10] <- 0
      }
    }
    
    # Check if all extinct
    all_extinct <- TRUE
    for (col_dex in 3:12) {
      if ((out.df[nrow(out.df), col_dex] + out.df[nrow(out.df), col_dex + 10] + out.df[nrow(out.df), (col_dex + 22)] + out.df[nrow(out.df), (col_dex + 32)]) >= threshold_extinction) {
        all_extinct <- FALSE
      }
    }
    if (all_extinct) {
      break
    } else {
      endyear <- endyear + interyear
    }
    
    # If there are already ten strains and none are extinct
    if (length(which(out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 13:22] == 0)) == 10) {
      I_dex_setzero <- which((out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 13:22]) == min(out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 13:22]))[1] + 12
      VI_dex_setzero <- I_dex_setzero + 22
      out.df[nrow(out.df), I_dex_setzero] <- 0
      out.df[nrow(out.df), I_dex_setzero - 10] <- 0
      out.df[nrow(out.df), VI_dex_setzero] <- 0
      out.df[nrow(out.df), VI_dex_setzero - 10] <- 0
    }
    
    # Introduce new strain
    I_dex_tointroduce <- which((out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 13:22]) == min(out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 13:22]))[1] + 12
    tobase_newvir <- mean(unlist(virstrats[((endyear / interyear) - 1),]), na.rm=T)
    min_newvir <- ifelse((tobase_newvir - 10) < 2, 2, (tobase_newvir - 10))
    max_newvir <- ifelse((tobase_newvir + 10) > 100, 100, (tobase_newvir + 10))
    new_vir <- runif(1, min=min_newvir, max=max_newvir)
    fbet_new <- (c1 * (new_vir)^c2) / pop_size
    p_new <- (new_vir) / 100
    parameters[I_dex_tointroduce -10 -2] <- fbet_new
    parameters[I_dex_tointroduce -10 -2 +10] <- p_new
    init <- c(fS=out.df$fS[nrow(out.df)],
              fE1=out.df$fE1[nrow(out.df)], fE2=out.df$fE2[nrow(out.df)], fE3=out.df$fE3[nrow(out.df)], fE4=out.df$fE4[nrow(out.df)], fE5=out.df$fE5[nrow(out.df)], fE6=out.df$fE6[nrow(out.df)], fE7=out.df$fE7[nrow(out.df)], fE8=out.df$fE8[nrow(out.df)], fE9=out.df$fE9[nrow(out.df)], fE10=out.df$fE10[nrow(out.df)], 
              fI1=out.df$fI1[nrow(out.df)], fI2=out.df$fI2[nrow(out.df)], fI3=out.df$fI3[nrow(out.df)], fI4=out.df$fI4[nrow(out.df)], fI5=out.df$fI5[nrow(out.df)], fI6=out.df$fI6[nrow(out.df)], fI7=out.df$fI7[nrow(out.df)], fI8=out.df$fI8[nrow(out.df)], fI9=out.df$fI9[nrow(out.df)], fI10=out.df$fI10[nrow(out.df)], 
              fR=out.df$fR[nrow(out.df)], fV=out.df$fV[nrow(out.df)], 
              fVE1=out.df$fVE1[nrow(out.df)], fVE2=out.df$fVE2[nrow(out.df)], fVE3=out.df$fVE3[nrow(out.df)], fVE4=out.df$fVE4[nrow(out.df)], fVE5=out.df$fVE5[nrow(out.df)], fVE6=out.df$fVE6[nrow(out.df)], fVE7=out.df$fVE7[nrow(out.df)], fVE8=out.df$fVE8[nrow(out.df)], fVE9=out.df$fVE9[nrow(out.df)], fVE10=out.df$fVE10[nrow(out.df)],
              fVI1=out.df$fVI1[nrow(out.df)], fVI2=out.df$fVI2[nrow(out.df)], fVI3=out.df$fVI3[nrow(out.df)], fVI4=out.df$fVI4[nrow(out.df)], fVI5=out.df$fVI5[nrow(out.df)], fVI6=out.df$fVI6[nrow(out.df)], fVI7=out.df$fVI7[nrow(out.df)], fVI8=out.df$fVI8[nrow(out.df)], fVI9=out.df$fVI9[nrow(out.df)], fVI10=out.df$fVI10[nrow(out.df)])
    init[I_dex_tointroduce - 1] <- 1
    
    # Write down current strains, cycling through ps
    for (param_dex in 11:20) {
      if (parameters[param_dex] != 0) {
        virstrats[(endyear / interyear), param_dex - 10] <- parameters[param_dex] * 100
      } else {
        virstrats[(endyear / interyear), param_dex - 10] <- NA
      }
    }
    
    # Run sim
    time <- seq(0, interyear*365, by = 0.1)
    out <- ode(y=init, times=time, func=eqn, parms=parameters)
    out.df <- as.data.frame(out)
  }
  # Last iteration
  if (!all_extinct) {
    if (endyear != maxyear) {
      stop('Error in end of simulation.')
    }
    # Set to 0 for those that are extinct
    for (col_dex in 3:12) {
      if ((out.df[nrow(out.df), col_dex] + out.df[nrow(out.df), (col_dex + 10)] + out.df[nrow(out.df), (col_dex + 22)] + out.df[nrow(out.df), (col_dex + 32)]) < threshold_extinction) {
        out.df[nrow(out.df), col_dex] <- 0
        out.df[nrow(out.df), (col_dex + 10)] <- 0
        out.df[nrow(out.df), (col_dex + 22)] <- 0
        out.df[nrow(out.df), (col_dex + 32)] <- 0
        
        # Set parameters for this strain to 0 as well
        parameters[col_dex - 2] <- 0
        parameters[col_dex - 2 + 10] <- 0
      }
    }
    # Write down current strains, cycling through ps
    for (param_dex in 11:20) {
      if (parameters[param_dex] != 0) {
        virstrats[(endyear / interyear) + 1, param_dex - 10] <- parameters[param_dex] * 100
      } else {
        virstrats[(endyear / interyear) + 1, param_dex - 10] <- NA
      }
    }
  }
  return(virstrats)
}

# Run nsim number of simulations in parallel -----------------------------------
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix <- foreach(i=1:nsim, .combine=rbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = multistrainSim_mod2()
  # write.csv(tempMatrix, paste0("~/virEvol/scratch/", i, ".csv")) # for debugging purposes
  tempMatrix
}
stopCluster(cl)
write.csv(finalMatrix, paste0('~/virEvol/code_output/multistrain_res/mod2_interyear', interyear_input, '_maxyear', maxyear_input, '.csv'))
rm(list = ls())
