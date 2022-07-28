# ------------------------------------------------------------------------------
# Model 3: SEIR with vaccination and migration
# ------------------------------------------------------------------------------
# Source functional scripts and set seed and nsim ------------------------------
source('~/virEvol/code/func/gen_parameters.R')
source('~/virEvol/code/func/multistrain_eqns.R')
set.seed(0)
nsim <- 1000
interyear_input <- 1
maxyear_input <- 100

# Set model 3 specific parameters and functions --------------------------------
# Initial susceptible population in farms
fS_init = (pop_size * 1/2) - 1
# Initial susceptible population in markets
mS_init = (pop_size * 1/2) - 1
# Initial strain 1 infectious population in farms
fI1_init = 1
# Initial strain 1 infectious population in markets
mI1_init = 1 
# Percent of susceptible chickens vaccinated in each time period
perc_vax = 0.33 
# Time that perc_vax is vaccinated
inter_vax_time = 120 
# Vaccination rate of chickens of farms per susceptible chicken of farm per day
v = perc_vax / inter_vax_time 
# Rate of loss of immunity due to vaccination per chicken per day
v_hat = 1 / 120
# Percent sold in interval
perc_sold_per_farm = 0.33
# Days between successive sales of chickens of a farm
inter_sell_time_per_farm = 120 
# Migration rate of chickens from farms to markets per chicken per day, if unvaccinated
m_fm = perc_sold_per_farm / inter_sell_time_per_farm
# Frac of base migration rate for vaccinated chickens
diff_vax = 1
# Migration rate of chickens from farms to markets per chicken, if vaccinated
m_fm_vax = (perc_sold_per_farm * diff_vax) / inter_sell_time_per_farm
# Migration rate of chickens from markets to farms per chicken per day
m_mf = 1 / 60
# Ratio of contact rate in markets vs. farms
bet_mf_ratio = 10
# Threshold value for extinction
threshold_extinction = 3
# Assign model 3
eqn <- multi_eqn_mod3

# Create function that will output the virulence strategies every interyear ----
multistrainSim_mod3_mpatch <- function(interyear=interyear_input, maxyear=maxyear_input) {
  # Check that maxyear divisible by interyear
  if (maxyear %% interyear != 0) {
    stop('Maxyear must be divisible by interyear in order for while loop to work.')
  }
  
  # Randomly select first strain
  res_vir <- runif(1, min=2, max=100)
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  mbet1 <- fbet1 * bet_mf_ratio
  p_1 <- (res_vir) / 100
  
  # Parameters
  parameters <- c(fbet1=fbet1, fbet2=0, fbet3=0, fbet4=0, fbet5=0, fbet6=0, fbet7=0, fbet8=0, fbet9=0, fbet10=0,
                  p_1=p_1, p_2=0, p_3=0, p_4=0, p_5=0, p_6=0, p_7=0, p_8=0, p_9=0, p_10=0,
                  mbet1=mbet1, mbet2=0, mbet3=0, mbet4=0, mbet5=0, mbet6=0, mbet7=0, mbet8=0, mbet9=0, mbet10=0,
                  sig=sig, gamm=gamm,
                  mort=mort, b=b, nat_mort=nat_mort,
                  v=v, v_hat=v_hat,
                  m_fm=m_fm, m_fm_vax=m_fm_vax, m_mf=m_mf)
  
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
            fVI1=0, fVI2=0, fVI3=0, fVI4=0, fVI5=0, fVI6=0, fVI7=0, fVI8=0, fVI9=0, fVI10=0,
            mS=mS_init,
            mE1=0, mE2=0, mE3=0, mE4=0, mE5=0, mE6=0, mE7=0, mE8=0, mE9=0, mE10=0,
            mI1=mI1_init, mI2=0, mI3=0, mI4=0, mI5=0, mI6=0, mI7=0, mI8=0, mI9=0, mI10=0,
            mR=0, mV=0, 
            mVE1=0, mVE2=0, mVE3=0, mVE4=0, mVE5=0, mVE6=0, mVE7=0, mVE8=0, mVE9=0, mVE10=0,
            mVI1=0, mVI2=0, mVI3=0, mVI4=0, mVI5=0, mVI6=0, mVI7=0, mVI8=0, mVI9=0, mVI10=0)
  time <- seq(0, interyear*365, by = 0.1)
  out <- ode(y=init, times=time, func=eqn, parms=parameters)
  out.df <- as.data.frame(out)
  while (endyear < maxyear) {
    # Set to 0 for those that are extinct
    for (col_dex in 3:12) {
      if ((out.df[nrow(out.df), col_dex] + out.df[nrow(out.df), (col_dex + 10)] + 
           out.df[nrow(out.df), (col_dex + 22)] + out.df[nrow(out.df), (col_dex + 32)] +
           out.df[nrow(out.df), (col_dex + 43)] + out.df[nrow(out.df), (col_dex + 53)]) < threshold_extinction) {
        out.df[nrow(out.df), col_dex] <- 0
        out.df[nrow(out.df), (col_dex + 10)] <- 0
        out.df[nrow(out.df), (col_dex + 22)] <- 0
        out.df[nrow(out.df), (col_dex + 32)] <- 0
        out.df[nrow(out.df), (col_dex + 43)] <- 0
        out.df[nrow(out.df), (col_dex + 53)] <- 0
        
        # Set parameters for this strain to 0 as well (including market beta)
        parameters[col_dex - 2] <- 0
        parameters[col_dex - 2 + 10] <- 0
        parameters[col_dex - 2 + 20] <- 0
      }
    }
    
    # Check if all extinct
    all_extinct <- TRUE
    for (col_dex in 3:12) {
      if ((out.df[nrow(out.df), col_dex] + out.df[nrow(out.df), col_dex + 10] + 
           out.df[nrow(out.df), (col_dex + 22)] + out.df[nrow(out.df), (col_dex + 32)] +
           out.df[nrow(out.df), (col_dex + 43)] + out.df[nrow(out.df), (col_dex + 53)] +
           out.df[nrow(out.df), (col_dex + 65)] + out.df[nrow(out.df), (col_dex + 75)]) >= threshold_extinction) {
        all_extinct <- FALSE
      }
    }
    if (all_extinct) {
      break
    } else {
      endyear <- endyear + interyear
    }
    
    # If there are already ten strains and none are extinct (exclude exposed in this calculation)
    if (length(which(out.df[nrow(out.df), 13:22] + out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 56:65] + out.df[nrow(out.df), 78:87] == 0)) == 10) {
      I_dex_setzero <- which((out.df[nrow(out.df), 13:22] + out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 56:65] + out.df[nrow(out.df), 78:87]) == 
                               min(out.df[nrow(out.df), 13:22] + out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 56:65] + out.df[nrow(out.df), 78:87]))[1] + 12
      VI_dex_setzero <- I_dex_setzero + 22
      mI_dex_setzero <- I_dex_setzero + 43
      mVI_dex_setzero <- I_dex_setzero + 65
      out.df[nrow(out.df), I_dex_setzero] <- 0
      out.df[nrow(out.df), I_dex_setzero - 10] <- 0
      out.df[nrow(out.df), VI_dex_setzero] <- 0
      out.df[nrow(out.df), VI_dex_setzero - 10] <- 0
      out.df[nrow(out.df), mI_dex_setzero] <- 0
      out.df[nrow(out.df), mI_dex_setzero - 10] <- 0
      out.df[nrow(out.df), mVI_dex_setzero] <- 0
      out.df[nrow(out.df), mVI_dex_setzero - 10] <- 0
    }
    
    # Introduce new strain
    I_dex_tointroduce <- which((out.df[nrow(out.df), 13:22] + out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 56:65] + out.df[nrow(out.df), 78:87]) == 
                                 min(out.df[nrow(out.df), 13:22] + out.df[nrow(out.df), 35:44] + out.df[nrow(out.df), 56:65] + out.df[nrow(out.df), 78:87]))[1] + 12
    new_vir <- runif(1, min=2, max=100)
    fbet_new <- (c1 * (new_vir)^c2) / pop_size
    mbet_new <- fbet_new * bet_mf_ratio
    p_new <- (new_vir) / 100
    parameters[I_dex_tointroduce -10 -2] <- fbet_new
    parameters[I_dex_tointroduce -10 -2 +10] <- p_new
    parameters[I_dex_tointroduce -10 -2 +20] <- mbet_new
    init <- c(fS=out.df$fS[nrow(out.df)],
              fE1=out.df$fE1[nrow(out.df)], fE2=out.df$fE2[nrow(out.df)], fE3=out.df$fE3[nrow(out.df)], fE4=out.df$fE4[nrow(out.df)], fE5=out.df$fE5[nrow(out.df)], fE6=out.df$fE6[nrow(out.df)], fE7=out.df$fE7[nrow(out.df)], fE8=out.df$fE8[nrow(out.df)], fE9=out.df$fE9[nrow(out.df)], fE10=out.df$fE10[nrow(out.df)], 
              fI1=out.df$fI1[nrow(out.df)], fI2=out.df$fI2[nrow(out.df)], fI3=out.df$fI3[nrow(out.df)], fI4=out.df$fI4[nrow(out.df)], fI5=out.df$fI5[nrow(out.df)], fI6=out.df$fI6[nrow(out.df)], fI7=out.df$fI7[nrow(out.df)], fI8=out.df$fI8[nrow(out.df)], fI9=out.df$fI9[nrow(out.df)], fI10=out.df$fI10[nrow(out.df)], 
              fR=out.df$fR[nrow(out.df)], fV=out.df$fV[nrow(out.df)], 
              fVE1=out.df$fVE1[nrow(out.df)], fVE2=out.df$fVE2[nrow(out.df)], fVE3=out.df$fVE3[nrow(out.df)], fVE4=out.df$fVE4[nrow(out.df)], fVE5=out.df$fVE5[nrow(out.df)], fVE6=out.df$fVE6[nrow(out.df)], fVE7=out.df$fVE7[nrow(out.df)], fVE8=out.df$fVE8[nrow(out.df)], fVE9=out.df$fVE9[nrow(out.df)], fVE10=out.df$fVE10[nrow(out.df)],
              fVI1=out.df$fVI1[nrow(out.df)], fVI2=out.df$fVI2[nrow(out.df)], fVI3=out.df$fVI3[nrow(out.df)], fVI4=out.df$fVI4[nrow(out.df)], fVI5=out.df$fVI5[nrow(out.df)], fVI6=out.df$fVI6[nrow(out.df)], fVI7=out.df$fVI7[nrow(out.df)], fVI8=out.df$fVI8[nrow(out.df)], fVI9=out.df$fVI9[nrow(out.df)], fVI10=out.df$fVI10[nrow(out.df)],
              mS=out.df$mS[nrow(out.df)],
              mE1=out.df$mE1[nrow(out.df)], mE2=out.df$mE2[nrow(out.df)], mE3=out.df$mE3[nrow(out.df)], mE4=out.df$mE4[nrow(out.df)], mE5=out.df$mE5[nrow(out.df)], mE6=out.df$mE6[nrow(out.df)], mE7=out.df$mE7[nrow(out.df)], mE8=out.df$mE8[nrow(out.df)], mE9=out.df$mE9[nrow(out.df)], mE10=out.df$mE10[nrow(out.df)],
              mI1=out.df$mI1[nrow(out.df)], mI2=out.df$mI2[nrow(out.df)], mI3=out.df$mI3[nrow(out.df)], mI4=out.df$mI4[nrow(out.df)], mI5=out.df$mI5[nrow(out.df)], mI6=out.df$mI6[nrow(out.df)], mI7=out.df$mI7[nrow(out.df)], mI8=out.df$mI8[nrow(out.df)], mI9=out.df$mI9[nrow(out.df)], mI10=out.df$mI10[nrow(out.df)],
              mR=out.df$mR[nrow(out.df)], mV=out.df$mV[nrow(out.df)], 
              mVE1=out.df$mVE1[nrow(out.df)], mVE2=out.df$mVE2[nrow(out.df)], mVE3=out.df$mVE3[nrow(out.df)], mVE4=out.df$mVE4[nrow(out.df)], mVE5=out.df$mVE5[nrow(out.df)], mVE6=out.df$mVE6[nrow(out.df)], mVE7=out.df$mVE7[nrow(out.df)], mVE8=out.df$mVE8[nrow(out.df)], mVE9=out.df$mVE9[nrow(out.df)], mVE10=out.df$mVE10[nrow(out.df)],
              mVI1=out.df$mVI1[nrow(out.df)], mVI2=out.df$mVI2[nrow(out.df)], mVI3=out.df$mVI3[nrow(out.df)], mVI4=out.df$mVI4[nrow(out.df)], mVI5=out.df$mVI5[nrow(out.df)], mVI6=out.df$mVI6[nrow(out.df)], mVI7=out.df$mVI7[nrow(out.df)], mVI8=out.df$mVI8[nrow(out.df)], mVI9=out.df$mVI9[nrow(out.df)], mVI10=out.df$mVI10[nrow(out.df)])
    init[I_dex_tointroduce - 1] <- 1
    init[I_dex_tointroduce + 42] <- 1
    
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
  if (endyear != maxyear) {
    stop('Error in end of simulation.')
  }
  # Set to 0 for those that are extinct
  for (col_dex in 3:12) {
    if ((out.df[nrow(out.df), col_dex] + out.df[nrow(out.df), (col_dex + 10)] + 
         out.df[nrow(out.df), (col_dex + 22)] + out.df[nrow(out.df), (col_dex + 32)] +
         out.df[nrow(out.df), (col_dex + 43)] + out.df[nrow(out.df), (col_dex + 53)]) < threshold_extinction) {
      out.df[nrow(out.df), col_dex] <- 0
      out.df[nrow(out.df), (col_dex + 10)] <- 0
      out.df[nrow(out.df), (col_dex + 22)] <- 0
      out.df[nrow(out.df), (col_dex + 32)] <- 0
      out.df[nrow(out.df), (col_dex + 43)] <- 0
      out.df[nrow(out.df), (col_dex + 53)] <- 0
      
      # Set parameters for this strain to 0 as well (including market beta)
      parameters[col_dex - 2] <- 0
      parameters[col_dex - 2 + 10] <- 0
      parameters[col_dex - 2 + 20] <- 0
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
  tempMatrix = multistrainSim_mod3_mpatch()
  write.csv(tempMatrix, paste0("~/virEvol/scratch/", i, ".csv")) # for debugging purposes
  tempMatrix
}
stopCluster(cl)
write.csv(finalMatrix, paste0('~/virEvol/code_output/multistrain_res/mod3_mpatch_interyear', interyear_input, '_maxyear', maxyear_input, '.csv'))

