# Function to get results of simulated invasion for a given resident strain ----
# Legend: 0 = res strain wins
#         1 = invade strain wins
#         2 = res strain extinct before invader introduced
#         3 = both strains are not extinct at equilibrium
#         4 = both strains are extinct at after invader introduced
test_invade <- function(res_vir, invade_vir) {
  #print(paste0("res vir: ", res_vir))
  #print(paste0("invade vir: ", invade_vir))
  res <- NA
  
  # Strain specific parameters
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  mbet1 <- fbet1 * mfbet_ratio
  fbet2 <- (c1 * (invade_vir)^c2) / pop_size
  mbet2 <- fbet2 * mfbet_ratio
  p_1 <- ((res_vir) / 100)
  p_2 <- ((invade_vir) / 100)
  
  # Run resident strain until equilibrium (no vaccination)
  parameters <- c(fbet1=fbet1, fbet2=fbet2, mbet1=mbet1, mbet2=mbet2,
                  sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, 
                  m_fm=m_fm, m_fm_vax=m_fm_vax, m_mf=m_mf, 
                  v=v, v_hat=v_hat, theta=theta)
  init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR1=0, fR2=0, 
            fV=0, fV_E1=0, fV_I1=0, fV_E2=0, fV_I2=0,
            mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=0, mR1=0, mR2=0, 
            mV=0, mV_E1=0, mV_I1=0, mV_E2=0, mV_I2=0)
  time_eq1 <- seq(0, t_max_eq1, by = 1)
  out_eq1 <- ode(y=init, times=time_eq1, eqn, parms=parameters)
  out_eq1.df <- as.data.frame(out_eq1)
  #plot.out.df(out_eq1.df)
  
  # If resident has not gone extinct
  if ((out_eq1.df$fE1[nrow(out_eq1.df)] + 
       out_eq1.df$fI1[nrow(out_eq1.df)] + 
       out_eq1.df$mE1[nrow(out_eq1.df)] +
       out_eq1.df$mI1[nrow(out_eq1.df)] +
       out_eq1.df$fV_E1[nrow(out_eq1.df)] +
       out_eq1.df$fV_I1[nrow(out_eq1.df)] + 
       out_eq1.df$mV_E1[nrow(out_eq1.df)] +
       out_eq1.df$mV_I1[nrow(out_eq1.df)]) < threshold_extinction) {
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
    #plot.out.df(out_eq2.df)
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
    if (num_EI_res >= threshold_extinction) { # if resident is not extinct
      if (num_EI_invader < threshold_extinction) { # if invader is extinct
        res <- 0 # resident wins
      } else { # if invader is not extinct
        res <- 3 # coexistence
      }
    } else { # if resident is extinct
      if (num_EI_invader >= threshold_extinction) { # if invader is not extinct
        res <- 1 # invader wins
      } else { # if invader is extinct
        res <- 4 # both extinct
      }
    }
  }
  return(c(res))
}