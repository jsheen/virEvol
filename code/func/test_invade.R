# Test invasion function -------------------------------------------------------
# Legend: 0 = Resident wins (resident above extinction threshold, invader is not above extinction threshold)
#         1 = Invader wins (invader above extinction threshold, resident is not above extinction threshold)
#         2 = If there is extinction of the resident strain before invader introduced
#         3 = Coexistence (both resident and invader are above extinction threshold)
#         4 = If invader is extinct (both resident and invader are below extinction threshold)
#         5 = Error in code
test_invade <- function(res_vir, invade_vir, markets) {
  # Strain specific parameters
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  fbet2 <- (c1 * (invade_vir)^c2) / pop_size
  p_1 <- ((res_vir) / 100)
  p_2 <- ((invade_vir) / 100)
  
  # If markets are in model
  if (markets) {
    mbet1 <- fbet1 * mfbet_ratio
    mbet2 <- fbet2 * mfbet_ratio
  }
  
  # Run resident strain until equilibrium (no vaccination)
  time_eq1 <- seq(0, t_max_eq1, by = 1)
  out_eq1 <- ode(y=init, times=time_eq1, eqn, parms=parameters)
  out_eq1.df <- as.data.frame(out_eq1)
  
  # For debugging
  # plot.out.df(out_eq1.df)
  
  # Output result for res_vir and invade_vir combination
  res <- 5
  if ((out_eq1.df$fE1[nrow(out_eq1.df)] + 
       out_eq1.df$fI1[nrow(out_eq1.df)] + 
       out_eq1.df$mE1[nrow(out_eq1.df)] +
       out_eq1.df$mI1[nrow(out_eq1.df)] +
       out_eq1.df$fV_E1[nrow(out_eq1.df)] +
       out_eq1.df$fV_I1[nrow(out_eq1.df)] + 
       out_eq1.df$mV_E1[nrow(out_eq1.df)] +
       out_eq1.df$mV_I1[nrow(out_eq1.df)]) < threshold_extinction) {
    # If there is extinction of the resident strain before invader introduced
    res <- 2
  } else {
    eq2_init <- c(fS=out_eq1.df$fS[nrow(out_eq1.df)], fE1=out_eq1.df$fE1[nrow(out_eq1.df)], 
                  fE2=out_eq1.df$fE2[nrow(out_eq1.df)], fI1=out_eq1.df$fI1[nrow(out_eq1.df)], 
                  fI2=1, fR1=out_eq1.df$fR1[nrow(out_eq1.df)], 
                  fR2=out_eq1.df$fR2[nrow(out_eq1.df)], 
                  fV=out_eq1.df$fV[nrow(out_eq1.df)], fV_E1=out_eq1.df$fV_E1[nrow(out_eq1.df)], 
                  fV_I1=out_eq1.df$fV_I1[nrow(out_eq1.df)], fV_E2=out_eq1.df$fV_E2[nrow(out_eq1.df)], 
                  fV_I2=out_eq1.df$fV_I2[nrow(out_eq1.df)],
                  mS=out_eq1.df$mS[nrow(out_eq1.df)], mE1=out_eq1.df$mE1[nrow(out_eq1.df)], 
                  mE2=out_eq1.df$mE2[nrow(out_eq1.df)], mI1=out_eq1.df$mI1[nrow(out_eq1.df)], 
                  mI2=1, mR1=out_eq1.df$mR1[nrow(out_eq1.df)], 
                  mR2=out_eq1.df$mR2[nrow(out_eq1.df)], 
                  mV=out_eq1.df$mV[nrow(out_eq1.df)], mV_E1=out_eq1.df$mV_E1[nrow(out_eq1.df)], 
                  mV_I1=out_eq1.df$mV_I1[nrow(out_eq1.df)], mV_E2=out_eq1.df$mV_E2[nrow(out_eq1.df)], 
                  mV_I2=out_eq1.df$mV_I2[nrow(out_eq1.df)])
    time_eq2 <- seq(0, t_max_eq2, by = 1)
    out_eq2 <- ode(y=eq2_init, times=time_eq2, eqn, parms=parameters)
    out_eq2.df <- as.data.frame(out_eq2)
    
    # For debugging purposes
    # plot.out.df(out_eq2.df)
    
    num_EI_res <- out_eq2.df$fE1[nrow(out_eq2.df)] + 
      out_eq2.df$fI1[nrow(out_eq2.df)] + 
      out_eq2.df$fV_E1[nrow(out_eq2.df)] + 
      out_eq2.df$fV_I1[nrow(out_eq2.df)] + 
      out_eq2.df$mE1[nrow(out_eq2.df)] +
      out_eq2.df$mI1[nrow(out_eq2.df)] + 
      out_eq2.df$mV_E1[nrow(out_eq2.df)] +
      out_eq2.df$mV_I1[nrow(out_eq2.df)]
    num_EI_invader <- out_eq2.df$fE2[nrow(out_eq2.df)] + 
      out_eq2.df$fI2[nrow(out_eq2.df)] + 
      out_eq2.df$fV_E2[nrow(out_eq2.df)] + 
      out_eq2.df$fV_I2[nrow(out_eq2.df)] + 
      out_eq2.df$mE2[nrow(out_eq2.df)] +
      out_eq2.df$mI2[nrow(out_eq2.df)] + 
      out_eq2.df$mV_E2[nrow(out_eq2.df)] +
      out_eq2.df$mV_I2[nrow(out_eq2.df)]
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