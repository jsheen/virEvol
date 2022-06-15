# ------------------------------------------------------------------------------
# Space modifies predictions of virulence evolution through differential immunization
# status and contact rate in a backyard poultry farm, live market system
#
# @author: Justin Sheen
# @description: a two-patch, deterministic model that shows the dependency of virulence 
#               evolution on space due to differential immunization status and contact rate.
#               With this model, we seek to understand whether, from a modeling standpoint,
#               a higher virulent form of the virus would be permitted, and whether this
#               higher virulent form will only be able to circulate within farms. The rate 
#               of migration may limit the amount of space the virulent virus is able to 
#               take over as well.
# ------------------------------------------------------------------------------

# Assumptions ------------------------------------------------------------------
# - We assume density independent transmission rates
# - We assume that there is no transmission of the vaccine, and that there are no other lower-virulent strains of the virus circulating
# - We assume that probability of death is the same between strain 1 and strain 2 after first infection
# - Vaccination completely protects chickens from mortality, and is only conducted in farms
# - We assume no super-infection in the system

# Load libraries ---------------------------------------------------------------
library(deSolve)
library(foreach)
library(doParallel)
library(pracma)
library(plot.matrix)

# Fixed parameters of the model ------------------------------------------------
t_max_eq1 = 2e3 # time of simulation to first equilibrium (days)
t_max_eq2 = 2e3 # time of simulation to second equilibrium (days)
pop_size = 1e6 # population size
fS_init = (pop_size * 1/2) - 1 # initial susceptible population in farms
mS_init = (pop_size * 1/2) - 1 # initial susceptible population in markets
fI1_init = 1 # initial strain 1 infectious population in farms
mI1_init = 1 # initial strain 1 infectious population in markets
fI2_init = 0 # initial strain 2 infectious population in farms
mI2_init = 0 # initial strain 2 infectious population in markets
sig = 1 / 5 # transition rate of infectiousness per chicken per day
gamm = 1 / 10 # transition rate of recovery per chicken per day
mort = 1 / 4 # disease mortality rate per chicken per day
nat_mort = 1 / 730 # natural mortality rate per chicken per day
b = 1 / 120 # average birth rate of chickens for raising
perc_sold_per_farm = 0.3 # percent sold in interval
inter_sell_time_per_farm = 120 # days between successive sales of chickens of a farm
m_fm = perc_sold_per_farm / inter_sell_time_per_farm # migration rate of chickens from farms to markets per chicken per day
m_mf = (1 / 7) # migration rate of chickens from markets to farms per chicken per day
perc_vax = 0.1 # percent vaccinated at each campaign
inter_vax_time = 120 # time that perc_vax is vaccinated
v = perc_vax / inter_vax_time # vaccination rate of chickens of farms per susceptible chicken of farm per day
v_hat = (1 / 126) # rate of loss of immunity due to vaccination per chicken per day
theta = (1 / 126) # rate of loss of immunity due to previous infection per chicken per day
vir_steps = seq(0.01, 95.01, 5)
mfbet_ratio = 5

# Plotting function ------------------------------------------------------------
plot.out.df <- function(out.df) {
  # Population size
  plot(out.df$time, rowSums(out.df[-c(1)]), main='Population size', ylab='N', xlab='days', type='l')
  
  # Susceptible pool divided into markets and farms
  plot(out.df$time, out.df$fS, col='black', type='l', ylim=c(0, max(out.df$fS, out.df$mS)), main='Susceptible divided into markets and farms', ylab='S', xlab='days')
  lines(out.df$time, out.df$mS, col='red')
  legend('topright', legend=c("farms", "markets"),
         col=c("black", "red"), lty=1, cex=0.8)
  
  # Infectious pool divided into markets and farms
  plot(out.df$time, out.df$fI1, type='l', ylim=c(0, max(out.df$fI1, out.df$mI1, out.df$fI2, out.df$mI2)),
       main='Incubation and Infectious divided into markets and farms', ylab='I', xlab='days', col='orange', lty=1)
  lines(out.df$time, out.df$mI1, col='orange', lty=2)
  lines(out.df$time, out.df$fI2, col='blue', lty=1)
  lines(out.df$time, out.df$mI2, col='blue', lty=2)
  legend('topright', legend=c("farms_strain1", "markets_strain1", "farms_strain2", "markets_strain2"),
         col=c("orange", "orange", "blue", "blue"), lty=c(1,2,1,2), cex=0.8)
  
  # Comparing strains 1 and 2
  plot(out.df$time, out.df$mI1 + out.df$fI1, col='orange', type='l', ylim=c(0, max((out.df$fI1 + out.df$mI1), (out.df$fI2 + out.df$mI2))),
       main='Strain 1 vs. Strain 2', ylab='I', xlab='days')
  lines(out.df$time, out.df$mI2 + out.df$fI2, col='blue')
  legend('topright', legend=c("strain 1", "strain 2"),
         col=c("orange", "blue"), lty=1, cex=0.8)
  
  # Plot recovered pool
  plot(out.df$time, out.df$mR1 + out.df$fR1, col='orange', type='l', ylim=c(0, max((out.df$fR1 + out.df$mR1), (out.df$fR2 + out.df$mR2))),
       main='Strain 1 vs. Strain 2', ylab='R', xlab='days')
  lines(out.df$time, out.df$mR2 + out.df$fR2, col='blue')
  legend('topright', legend=c("strain 1", "strain 2"),
         col=c("orange", "blue"), lty=1, cex=0.8)
  
  # Vaccination divided into markets and farms
  plot(out.df$time, out.df$mV, col='red', type='l', ylim=c(0, max(out.df$mV, out.df$fV)), main='Vaccination divided into markets and farms', ylab='V', xlab='days')
  lines(out.df$time, out.df$fV, col='black')
  lines(out.df$time, out.df$mV_I1 + out.df$mV_I2, col='red', lty=2)
  lines(out.df$time, out.df$fV_I1 + out.df$fV_I2, col='black', lty=2)
  legend('topright', legend=c("farms", "markets", 'farms (infectious)', 'markets (infectious)'),
         col=c("black", "red", 'black', 'red'), lty=c(1, 1, 2, 2), cex=0.8)
}

# # Plot transmission-mortality tradeoff curve -----------------------------------
virulences <- seq(0.01, 100, 0.1)
morts <- ((virulences * 0.5) / 100) + 0.5
betas <- (0.5 * (virulences)^0.33) + 0.1 # beta is, roughly, the number of chickens a single infectious chicken infects in a month
plot(morts, betas, type='l')

# Model equations --------------------------------------------------------------
eqn <- function(time, state, parameters){
  with(as.list(c(state, parameters)),{
    # Backyard poultry farms
    dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR1 + fR2 + fV + fV_E1 + fV_I1 + fV_E2 + fV_I2)*(1 - ((fS + fE1 + fE2 + fI1 + fI2 + fR1 + fR2 + fV + fV_E1 + fV_I1 + fV_E2 + fV_I2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (fbet1*fS*fI1) -(fbet2*fS*fI2) -
      (fbet1*fS*fV_I1) -(fbet2*fS*fV_I2) -
      v*fS +v_hat*fV -
      m_fm*fS +m_mf*mS -
      nat_mort*fS +
      theta*fR1 + theta*fR2
    dfE1 = (fbet1*fS*fI1) +(fbet1*fS*fV_I1) -
      sig*fE1 -
      m_fm*fE1 +m_mf*mE1 -
      nat_mort*fE1
    dfE2 = (fbet2*fS*fI2) +(fbet2*fS*fV_I2) -
      sig*fE2 -
      m_fm*fE2 +m_mf*mE2 -
      nat_mort*fE2
    dfI1 = sig*fE1 -
      gamm*(1-p_1)*fI1 -mort*p_1*fI1 -
      m_fm*fI1 +m_mf*mI1 -
      nat_mort*fI1
    dfI2 = sig*fE2 -
      gamm*(1-p_2)*fI2 -mort*p_2*fI2 -
      m_fm*fI2 + m_mf*mI2 -
      nat_mort*fI2
    dfR1 = gamm*(1-p_1)*fI1 -
      theta*fR1 -
      m_fm*fR1 +m_mf*mR1 -
      nat_mort*fR1
    dfR2 = gamm*(1-p_2)*fI2 -
      theta*fR2 -
      m_fm*fR2 +m_mf*mR2 -
      nat_mort*fR2
    dfV = v*(fS) -v_hat*fV -
      (fbet1*fV*fI1) -(fbet2*fV*fI2) -
      (fbet1*fV*fV_I1) -(fbet2*fV*fV_I2) +
      gamm*fV_I1 +gamm*fV_I2 -
      m_fm*fV + m_mf*mV -
      nat_mort*fV
    dfV_E1 = (fbet1*fV*fI1) +(fbet1*fV*fV_I1) -
      sig*fV_E1 -
      m_fm*fV_E1 +m_mf*mV_E1 -
      nat_mort*fV_E1
    dfV_I1 = sig*fV_E1 -
      gamm*fV_I1 -
      m_fm*fV_I1 +m_mf*mV_I1 -
      nat_mort*fV_I1
    dfV_E2 = (fbet2*fV*fI2) +(fbet2*fV*fV_I2) -
      sig*fV_E2 -
      m_fm*fV_E2 +m_mf*mV_E2 -
      nat_mort*fV_E2
    dfV_I2 = sig*fV_E2 -
      gamm*fV_I2 -
      m_fm*fV_I2 +m_mf*mV_I2 -
      nat_mort*fV_I2
    # Live bird markets
    dmS = b*(mS + mE1 + mE2 + mI1 + mI2 + mR1 + mR2 + mV + mV_E1 + mV_I1 + mV_E2 + mV_I2)*(1 - ((mS + mE1 + mE2 + mI1 + mI2 + mR1 + mR2 + mV + mV_E1 + mV_I1 + mV_E2 + mV_I2) / ((b / (b - nat_mort)) * (pop_size)))) -
      (mbet1*mS*mI1) -(mbet2*mS*mI2) -
      (mbet1*mS*mV_I1) -(mbet2*mS*mV_I2) -
      m_mf*mS +m_fm*fS -
      nat_mort*mS +
      theta*mR1 + theta*mR2
    dmE1 = (mbet1*mS*mI1) +(mbet1*mS*mV_I1) -
      sig*mE1 -
      m_mf*mE1 +m_fm*fE1 -
      nat_mort*mE1
    dmE2 = (mbet2*mS*mI2) +(mbet2*mS*mV_I2)  -
      sig*mE2 -
      m_mf*mE2 +m_fm*fE2 -
      nat_mort*mE2
    dmI1 = sig*mE1 -
      gamm*(1-p_1)*mI1 -mort*p_1*mI1 -
      m_mf*mI1 +m_fm*fI1 -
      nat_mort*mI1
    dmI2 = sig*mE2 -
      gamm*(1-p_2)*mI2 -mort*p_2*mI2 -
      m_mf*mI2 +m_fm*fI2 -
      nat_mort*mI2
    dmR1 = gamm*(1-p_1)*mI1 -
      theta*mR1 -
      m_mf*mR1 +m_fm*fR1 -
      nat_mort*mR1
    dmR2 = gamm*(1-p_2)*mI2 -
      theta*mR2 -
      m_mf*mR2 +m_fm*fR2 -
      nat_mort*mR2
    dmV = -v_hat*mV -
      (mbet1*mV*mI1) -(mbet2*mV*mI2) -
      (mbet1*mV*mV_I1) -(mbet2*mV*mV_I2) +
      gamm*mV_I1 + gamm*mV_I2 +
      m_fm*fV -m_mf*mV -
      nat_mort*mV
    dmV_E1 = (mbet1*mV*mI1) +(mbet1*mV*mV_I1) -
      sig*mV_E1 +
      m_fm*fV_E1 -m_mf*mV_E1 -
      nat_mort*mV_E1
    dmV_I1 = sig*mV_E1 -
      gamm*mV_I1 +
      m_fm*fV_I1 -m_mf*mV_I1 -
      nat_mort*mV_I1
    dmV_E2 = (mbet2*mV*mI2) +(mbet2*mV*mV_I2) -
      sig*mV_E2 +
      m_fm*fV_E2 -m_mf*mV_E2 -
      nat_mort*mV_E2
    dmV_I2 = sig*mV_E2 -
      gamm*mV_I2 +
      m_fm*fV_I2 -m_mf*mV_I2 -
      nat_mort*mV_I2
    return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR1, dfR2,
                  dfV, dfV_E1, dfV_I1, dfV_E2, dfV_I2, 
                  dmS, dmE1, dmE2, dmI1, dmI2, dmR1, dmR2, 
                  dmV, dmV_E1, dmV_I1, dmV_E2, dmV_I2)))})}

# Function to get results of simulated invasion for a given resident strain ----
# Legend: 0 = res strain wins; 
#         1 = invade strain wins; 
#         2 = res strain extinct before invader introduced; 
#         3 = both strains are extinct at after invader introduced; 
#         4 = both strains are not extinct at equilibrium
test_invade <- function(res_vir, invade_vir) {
  print(paste0("res vir: ", res_vir))
  print(paste0("invade vir: ", invade_vir))
  res <- NA
  
  # Strain specific parameters
  fbet1 <- ((0.5 * (res_vir)^0.33) + 0.1) / pop_size
  mbet1 <- fbet1 * mfbet_ratio
  fbet2 <- ((0.5 * (invade_vir)^0.33) + 0.1) / pop_size
  mbet2 <- fbet2 * mfbet_ratio
  p_1 <- ((res_vir * 0.5) / 100) + 0.5
  p_2 <- ((invade_vir * 0.5) / 100) + 0.5
  
  parameters <- c(fbet1=fbet1, fbet2=fbet2, mbet1=mbet1, mbet2=mbet2,
                  sig=sig, gamm=gamm, mort=mort, p_1=p_1, p_2=p_2, b=b, 
                  m_fm=m_fm, m_mf=m_mf,
                  v=v, v_hat=v_hat, theta=theta)
  # Run resident strain until equilibrium
  init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=fI2_init, fR1=0, fR2=0, 
            fV=0, fV_E1=0, fV_I1=0, fV_E2=0, fV_I2=0,
            mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=mI2_init, mR1=0, mR2=0, 
            mV=0, mV_E1=0, mV_I1=0, mV_E2=0, mV_I2=0)
  time_eq1 <- seq(0, t_max_eq1, by = 1)
  out <- ode(y=init, times=time_eq1, eqn, parms=parameters)
  out.df <- as.data.frame(out)
  plot.out.df(out.df)
  # Introduce vaccination
  
  
  if (round(out.df$fE1[nrow(out.df)] + 
            out.df$fI1[nrow(out.df)] + 
            out.df$mE1[nrow(out.df)] +
            out.df$mI1[nrow(out.df)] +
            out.df$fV_E1[nrow(out.df)] +
            out.df$fV_I1[nrow(out.df)] + 
            out.df$mV_E1[nrow(out.df)] +
            out.df$mV_I1[nrow(out.df)]) < 1) {
    res <- 2
  } else {
    invade_init <- c(fS=out.df$fS[nrow(out.df)], fE1=out.df$fE1[nrow(out.df)], 
                     fE2=out.df$fE2[nrow(out.df)], fI1=out.df$fI1[nrow(out.df)], 
                     fI2=1, fR1=out.df$fR1[nrow(out.df)], 
                     fR2=out.df$fR2[nrow(out.df)], 
                     fV=out.df$fV[nrow(out.df)], fV_E1=out.df$fV_E1[nrow(out.df)], 
                     fV_I1=out.df$fV_I1[nrow(out.df)], fV_E2=out.df$fV_E2[nrow(out.df)], 
                     fV_I2=out.df$fV_I2[nrow(out.df)],
                     mS=out.df$mS[nrow(out.df)], mE1=out.df$mE1[nrow(out.df)], 
                     mE2=out.df$mE2[nrow(out.df)], mI1=out.df$mI1[nrow(out.df)], 
                     mI2=1, mR1=out.df$mR1[nrow(out.df)], 
                     mR2=out.df$mR2[nrow(out.df)], 
                     mV=out.df$mV[nrow(out.df)], mV_E1=out.df$mV_E1[nrow(out.df)], 
                     mV_I1=out.df$mV_I1[nrow(out.df)], mV_E2=out.df$mV_E2[nrow(out.df)], 
                     mV_I2=out.df$mV_I2[nrow(out.df)])
    time_eq2 <- seq(0, t_max_eq2, by = 1)
    out_invade <- ode(y=invade_init, times=time_eq2, eqn, parms=parameters)
    out_invade.df <- as.data.frame(out_invade)
    plot.out.df(out_invade.df)
    num_EI_res <- round(out_invade.df$fE1[nrow(out_invade.df)] + 
                          out_invade.df$fI1[nrow(out_invade.df)] + 
                          out_invade.df$fV_E1[nrow(out_invade.df)] + 
                          out_invade.df$fV_I1[nrow(out_invade.df)] + 
                          out_invade.df$mE1[nrow(out_invade.df)] +
                          out_invade.df$mI1[nrow(out_invade.df)] + 
                          out_invade.df$mV_E1[nrow(out_invade.df)] +
                          out_invade.df$mV_I1[nrow(out_invade.df)]) 
    num_EI_invader <- round(out_invade.df$fE2[nrow(out_invade.df)] + 
                            out_invade.df$fI2[nrow(out_invade.df)] + 
                            out_invade.df$fV_E2[nrow(out_invade.df)] + 
                            out_invade.df$fV_I2[nrow(out_invade.df)] + 
                            out_invade.df$mE2[nrow(out_invade.df)] +
                            out_invade.df$mI2[nrow(out_invade.df)] + 
                            out_invade.df$mV_E2[nrow(out_invade.df)] +
                            out_invade.df$mV_I2[nrow(out_invade.df)]) 
    if (num_EI_res > 0 & num_EI_invader > 0) {
      res <- 4
    } else if (num_EI_res == 0 & num_EI_invader == 0) {
      res <- 3
    } else if (num_EI_res > 0 & num_EI_invader == 0) {
      res <- 0
    } else if (num_EI_res == 0 & num_EI_invader > 0) {
      res <- 1
    }
  }
  return(c(res))
}

# Create PIP (and run in parallel) ---------------------------------------------
# vir_step combinations
combos <- list()
combos_dex <- 1
for (i in vir_steps) {
  for (j in vir_steps) {
    combos[[combos_dex]] <- c(i, j)
    combos_dex <- combos_dex + 1
  }
}

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload computer
registerDoParallel(cl)

finalMatrix <- foreach(i=combos, .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(res_vir=i[1], invade_vir=i[2])
  write.csv(tempMatrix, paste0("~/virEvol/scratch/", i[1], "_", i[2], ".csv"))
  tempMatrix
}
#stop cluster
stopCluster(cl)

# plot PIP ---------------------------------------------------------------------
pip <- matrix(finalMatrix, ncol=length(vir_steps), nrow=length(vir_steps), byrow=F)
pip <- pracma::flipud(pip) #columns stay in place, but now from bottom to top is increasing virulence
#pip <- ifelse((pip != 0 & pip != 1), NA, pip)
plot(pip)
#write.csv(pip, paste0('~/virEvol/res/', perc_sold_per_farm, '_', perc_vax, '.csv'))




