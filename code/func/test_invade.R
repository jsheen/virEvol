# Test invasion functions ------------------------------------------------------
# Legend: 0 = Resident wins (resident above extinction threshold, invader is not above extinction threshold)
#         1 = Invader wins (invader above extinction threshold, resident is not above extinction threshold)
#         2 = If there is extinction of the resident strain before invader introduced
#         3 = Coexistence (both resident and invader are above extinction threshold)
#         4 = If invader is extinct (both resident and invader are below extinction threshold)
#         5 = Error in code

# Helper function --------------------------------------------------------------
check_equil <- function(out.df, uniq_vals_threshold) {
  in_equil <- TRUE
  start_check <- which(out.df$time == (out.df$time[nrow(out.df)] - 10))
  for (col_dex in 2:ncol(out.df)) {
    if (length(unique(round(out.df[start_check:nrow(out.df),col_dex]))) > uniq_vals_threshold) {
      print(col_dex)
      print(length(unique(round(out.df[start_check:nrow(out.df),col_dex]))))
      print(unique(round(out.df[start_check:nrow(out.df),col_dex])))
      in_equil <- FALSE
    }
  }
  return(in_equil)
}

# Model 1: SEIR ----------------------------------------------------------------
test_invade_mod1 <- function(res_vir, invade_vir) {
  # Strain specific virulence parameters
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  fbet2 <- (c1 * (invade_vir)^c2) / pop_size
  p_1 <- ((res_vir) / 100)
  p_2 <- ((invade_vir) / 100)
  
  # Parameters
  parameters <- c(fbet1=fbet1, fbet2=fbet2, sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, nat_mort=nat_mort)
  
  # Run resident strain until equilibrium (no vaccination)
  init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR=0)
  time_eq1 <- seq(0, t_max_eq1, by = 0.1)
  out_eq1 <- ode(y=init, times=time_eq1, func=eqn, parms=parameters)
  out_eq1.df <- as.data.frame(out_eq1)
  
  # For debugging purposes
  plot.out.df.mod1(out_eq1.df)
  
  # Check equilibrium 1
  if (!check_equil(out_eq1.df, uniq_vals_threshold=1000)) {
    write.csv(list(), paste0('~/virEvol/scratch_equil/err_mod1_equil1_', res_vir, '_', invade_vir, '.csv'))
  }
  
  # Output result for res_vir and invade_vir combination
  res <- 5
  if ((out_eq1.df$fE1[nrow(out_eq1.df)] + 
       out_eq1.df$fI1[nrow(out_eq1.df)]) < threshold_extinction) {
    # If there is extinction of the resident strain before invader introduced
    res <- 2
  } else {
    eq2_init <- c(fS=out_eq1.df$fS[nrow(out_eq1.df)], 
                  fE1=out_eq1.df$fE1[nrow(out_eq1.df)], fE2=out_eq1.df$fE2[nrow(out_eq1.df)], 
                  fI1=out_eq1.df$fI1[nrow(out_eq1.df)], fI2=1, 
                  fR=out_eq1.df$fR[nrow(out_eq1.df)])
    
    # Run model until second equilibrium
    time_eq2 <- seq(0, t_max_eq2, by = 1)
    out_eq2 <- ode(y=eq2_init, times=time_eq2, func=eqn, parms=parameters)
    out_eq2.df <- as.data.frame(out_eq2)
    
    # For debugging purposes
    # plot.out.df.mod1(out_eq2.df)
    
    # Get number of resident and invader infected at time of second equilibrium
    num_EI_res <- out_eq2.df$fE1[nrow(out_eq2.df)] + out_eq2.df$fI1[nrow(out_eq2.df)]
    num_EI_invader <- out_eq2.df$fE2[nrow(out_eq2.df)] + out_eq2.df$fI2[nrow(out_eq2.df)]
    
    # Determine which result occurred at end of second equilibrium
    if (num_EI_res >= threshold_extinction) { # If resident is not extinct
      if (num_EI_invader < threshold_extinction) { # If invader is extinct
        # Resident wins (resident above extinction threshold, invader is not above extinction threshold)
        res <- 0
      } else { # If invader is not extinct
        # Coexistence (both resident and invader are above extinction threshold)
        res <- 3
      }
    } else { # If resident is extinct
      if (num_EI_invader >= threshold_extinction) { # If invader is not extinct
        # Invader wins (invader above extinction threshold, resident is not above extinction threshold)
        res <- 1
      } else { # If invader is extinct (both resident and invader are below extinction threshold)
        # Both extinct
        res <- 4
      }
    }
  }
  return(c(res))
}

# Model 2: SEIR with vaccination -----------------------------------------------
test_invade_mod2 <- function(res_vir, invade_vir) {
  # Strain specific virulence parameters
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  fbet2 <- (c1 * (invade_vir)^c2) / pop_size
  p_1 <- ((res_vir) / 100)
  p_2 <- ((invade_vir) / 100)
  
  # Parameters
  parameters <- c(fbet1=fbet1, fbet2=fbet2, sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, nat_mort=nat_mort, 
                  v=v, v_hat=v_hat)
  
  # Run resident strain until equilibrium (no vaccination)
  init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR=0, fV=0, fVE1=0, fVE2=0, fVI1=0, fVI2=0)
  time_eq1 <- seq(0, t_max_eq1, by = 0.1)
  out_eq1 <- ode(y=init, times=time_eq1, func=eqn, parms=parameters)
  out_eq1.df <- as.data.frame(out_eq1)
  
  # For debugging purposes
  # plot.out.df.mod2(out_eq1.df)
  
  # Check equilibrium 1
  if (!check_equil(out_eq1.df, uniq_vals_threshold=1000)) {
    write.csv(list(), paste0('~/virEvol/scratch_equil/err_mod2_equil1_', res_vir, '_', invade_vir, '.csv'))
  }
  
  # Output result for res_vir and invade_vir combination
  res <- 5
  if ((out_eq1.df$fE1[nrow(out_eq1.df)] + 
       out_eq1.df$fI1[nrow(out_eq1.df)] +
       out_eq1.df$fVE1[nrow(out_eq1.df)] + 
       out_eq1.df$fVI1[nrow(out_eq1.df)]) < threshold_extinction) {
    # If there is extinction of the resident strain before invader introduced
    res <- 2
  } else {
    eq2_init <- c(fS=out_eq1.df$fS[nrow(out_eq1.df)], 
                  fE1=out_eq1.df$fE1[nrow(out_eq1.df)], fE2=out_eq1.df$fE2[nrow(out_eq1.df)], 
                  fI1=out_eq1.df$fI1[nrow(out_eq1.df)], fI2=1, 
                  fR=out_eq1.df$fR[nrow(out_eq1.df)],
                  fV=out_eq1.df$fV[nrow(out_eq1.df)],
                  fVE1=out_eq1.df$fVE1[nrow(out_eq1.df)], fVE2=out_eq1.df$fVE2[nrow(out_eq1.df)],
                  fVI1=out_eq1.df$fVI1[nrow(out_eq1.df)], fVI2=out_eq1.df$fVI2[nrow(out_eq1.df)])
    
    # Run model until second equilibrium
    time_eq2 <- seq(0, t_max_eq2, by = 1)
    out_eq2 <- ode(y=eq2_init, times=time_eq2, func=eqn, parms=parameters)
    out_eq2.df <- as.data.frame(out_eq2)
    
    # For debugging purposes
    # plot.out.df.mod2(out_eq2.df)
    
    # Get number of resident and invader infected at time of second equilibrium
    num_EI_res <- out_eq2.df$fE1[nrow(out_eq2.df)] + out_eq2.df$fI1[nrow(out_eq2.df)] + out_eq2.df$fVE1[nrow(out_eq2.df)] + out_eq2.df$fVI1[nrow(out_eq2.df)]
    num_EI_invader <- out_eq2.df$fE2[nrow(out_eq2.df)] + out_eq2.df$fI2[nrow(out_eq2.df)] + out_eq2.df$fVE2[nrow(out_eq2.df)] + out_eq2.df$fVI2[nrow(out_eq2.df)]
    
    # Determine which result occurred at end of second equilibrium
    if (num_EI_res >= threshold_extinction) { # If resident is not extinct
      if (num_EI_invader < threshold_extinction) { # If invader is extinct
        # Resident wins (resident above extinction threshold, invader is not above extinction threshold)
        res <- 0
      } else { # If invader is not extinct
        # Coexistence (both resident and invader are above extinction threshold)
        res <- 3
      }
    } else { # If resident is extinct
      if (num_EI_invader >= threshold_extinction) { # If invader is not extinct
        # Invader wins (invader above extinction threshold, resident is not above extinction threshold)
        res <- 1
      } else { # If invader is extinct (both resident and invader are below extinction threshold)
        # Both extinct
        res <- 4
      }
    }
  }
  return(c(res))
}

# Model 3: SEIR with vaccination and migration ---------------------------------
test_invade_mod3 <- function(res_vir, invade_vir) {
  # Strain specific virulence parameters
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  mbet1 <- fbet1 * bet_mf_ratio
  fbet2 <- (c1 * (invade_vir)^c2) / pop_size
  mbet2 <- fbet2 * bet_mf_ratio
  p_1 <- ((res_vir) / 100)
  p_2 <- ((invade_vir) / 100)
  
  # Parameters
  parameters <- c(fbet1=fbet1, fbet2=fbet2, sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, nat_mort=nat_mort, 
                  v=v, v_hat=v_hat, m_fm=m_fm, m_fm_vax=m_fm_vax, m_mf=m_mf, mbet1=mbet1, mbet2=mbet2)
  
  # Run resident strain until equilibrium (no vaccination)
  init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR=0, fV=0, fVE1=0, fVE2=0, fVI1=0, fVI2=0,
            mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=0, mR=0, mV=0, mVE1=0, mVE2=0, mVI1=0, mVI2=0)
  time_eq1 <- seq(0, t_max_eq1, by = 0.1)
  out_eq1 <- ode(y=init, times=time_eq1, eqn, parms=parameters)
  out_eq1.df <- as.data.frame(out_eq1)
  
  # For debugging purposes
  # plot.out.df.mod3(out_eq1.df)
  
  # Check equilibrium 1
  if (!check_equil(out_eq1.df, uniq_vals_threshold=1000)) {
    write.csv(list(), paste0('~/virEvol/scratch_equil/err_mod3_equil1_', res_vir, '_', invade_vir, '.csv'))
  }
  
  # Output result for res_vir and invade_vir combination
  res <- 5
  if ((out_eq1.df$fE1[nrow(out_eq1.df)] + 
       out_eq1.df$fI1[nrow(out_eq1.df)] +
       out_eq1.df$fVE1[nrow(out_eq1.df)] + 
       out_eq1.df$fVI1[nrow(out_eq1.df)] +
       out_eq1.df$mE1[nrow(out_eq1.df)] + 
       out_eq1.df$mI1[nrow(out_eq1.df)] +
       out_eq1.df$mVE1[nrow(out_eq1.df)] + 
       out_eq1.df$mVI1[nrow(out_eq1.df)]) < threshold_extinction) {
    # If there is extinction of the resident strain before invader introduced
    res <- 2
  } else {
    eq2_init <- c(fS=out_eq1.df$fS[nrow(out_eq1.df)], 
                  fE1=out_eq1.df$fE1[nrow(out_eq1.df)], fE2=out_eq1.df$fE2[nrow(out_eq1.df)], 
                  fI1=out_eq1.df$fI1[nrow(out_eq1.df)], fI2=1,
                  fR=out_eq1.df$fR[nrow(out_eq1.df)], fV=out_eq1.df$fV[nrow(out_eq1.df)],
                  fVE1=out_eq1.df$fVE1[nrow(out_eq1.df)], fVE2=out_eq1.df$fVE2[nrow(out_eq1.df)], 
                  fVI1=out_eq1.df$fVI1[nrow(out_eq1.df)], fVI2=out_eq1.df$fVI2[nrow(out_eq1.df)],
                  mS=out_eq1.df$mS[nrow(out_eq1.df)], 
                  mE1=out_eq1.df$mE1[nrow(out_eq1.df)], mE2=out_eq1.df$mE2[nrow(out_eq1.df)], 
                  mI1=out_eq1.df$mI1[nrow(out_eq1.df)], mI2=1,
                  mR=out_eq1.df$mR[nrow(out_eq1.df)],
                  mV=out_eq1.df$mV[nrow(out_eq1.df)],
                  mVE1=out_eq1.df$mVE1[nrow(out_eq1.df)], mVE2=out_eq1.df$mVE2[nrow(out_eq1.df)], 
                  mVI1=out_eq1.df$mVI1[nrow(out_eq1.df)], mVI2=out_eq1.df$mVI2[nrow(out_eq1.df)])
    
    # Run model until second equilibrium
    time_eq2 <- seq(0, t_max_eq2, by = 1)
    out_eq2 <- ode(y=eq2_init, times=time_eq2, func=eqn, parms=parameters)
    out_eq2.df <- as.data.frame(out_eq2)
    
    # For debugging purposes
    # plot.out.df.mod3(out_eq2.df)
    
    # Get number of resident and invader infected at time of second equilibrium
    num_EI_res <- out_eq2.df$fE1[nrow(out_eq2.df)] + out_eq2.df$fI1[nrow(out_eq2.df)] + out_eq2.df$fVE1[nrow(out_eq2.df)] + out_eq2.df$fVI1[nrow(out_eq2.df)] +
      out_eq2.df$mE1[nrow(out_eq2.df)] + out_eq2.df$mI1[nrow(out_eq2.df)] + out_eq2.df$mVE1[nrow(out_eq2.df)] + out_eq2.df$mVI1[nrow(out_eq2.df)]
    num_EI_invader <- out_eq2.df$fE2[nrow(out_eq2.df)] + out_eq2.df$fI2[nrow(out_eq2.df)] + out_eq2.df$fVE2[nrow(out_eq2.df)] + out_eq2.df$fVI2[nrow(out_eq2.df)] +
      out_eq2.df$mE2[nrow(out_eq2.df)] + out_eq2.df$mI2[nrow(out_eq2.df)] + out_eq2.df$mVE2[nrow(out_eq2.df)] + out_eq2.df$mVI2[nrow(out_eq2.df)]
    
    # Determine which result occurred at end of second equilibrium
    if (num_EI_res >= threshold_extinction) { # If resident is not extinct
      if (num_EI_invader < threshold_extinction) { # If invader is extinct
        # Resident wins (resident above extinction threshold, invader is not above extinction threshold)
        res <- 0
      } else { # If invader is not extinct
        # Coexistence (both resident and invader are above extinction threshold)
        res <- 3
      }
    } else { # If resident is extinct
      if (num_EI_invader >= threshold_extinction) { # If invader is not extinct
        # Invader wins (invader above extinction threshold, resident is not above extinction threshold)
        res <- 1
      } else { # If invader is extinct (both resident and invader are below extinction threshold)
        # Both extinct
        res <- 4
      }
    }
  }
  return(c(res))
}
