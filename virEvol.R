# ------------------------------------------------------------------------------
# Space modifies predictions of virulence evolution through differential immunization
# status and contact rate in a backyard poultry farm, live market system
#
# @author: Justin Sheen
# @description: a two-patch, deterministic model that shows the dependency of virulence 
#               evolution on space due to differential immunization status and contact rate.
#               With this model, we seek to understand whether, from a modeling standpoint,
#               a higher virulent form of the virus would be permitted, and whether this
#               higher virulent form will only be able to circulate within farms; acting
#               as a reservoir that spills into the market population and creates transitory
#               epidemics. The rate of migration may limit the amount of space the virulent 
#               virus is able to take over.
# ------------------------------------------------------------------------------

# Assumptions ------------------------------------------------------------------
# - We assume density independent transmission rates
# - We assume that there is no transmission of the vaccine, and that there are no other lower-virulent strains of the virus circulating
# - We assume that probability of death is the same between strain 1 and strain 2 after first infection
# - Vaccination completely protects chickens from mortality, and is only conducted in farms

# Load libraries ---------------------------------------------------------------
library(deSolve)

# Parameters of the model ------------------------------------------------------
# Strain 1 is a lower virulence strain
# Strain 2 is a higher virulence strain
t_max = 2e4
pop_size = 1e6
fS_init = (pop_size / 2) - 1 # initial susceptible population in farms
mS_init = (pop_size / 2) - 1 # initial susceptible population in markets
fI1_init = 1 # initial strain 1 infectious population in farms
mI1_init = 1 # initial strain 1 infectious population in markets
fI2_init = 1 # initial strain 2 infectious population in farms
mI2_init = 1 # initial strain 2 infectious population in markets
stopif(fS_init + mS_init + fI1_init + fI2_init + mI1_init + mI2_init != pop_size)
sig = 1 / 5 # transition rate of infectiousness per chicken per day
gamm = 1 / 6.25 # transition rate of recovery per chicken per day
fbet1 = 0.2 / pop_size # transmission rate of strain 1 among farms per SI contact per day
mbet1 = 0.8 / pop_size # transmission rate of strain 1 among markets per SI contact per day
strain_2_multiplication_transmission = 1.001
strain_2_multiplication_mortality = 1.5
fbet2 = fbet1 * strain_2_multiplication_transmission # transmission rate of strain 2 among farms per SI contact per day 
mbet2 = mbet1 * strain_2_multiplication_transmission # transmission rate of strain 2 among markets per SI contact per day
p_1 = 0.6 # probability of death from NDV infection of strain 1 given chicken was never infected
p_2 = p_1 * strain_2_multiplication_mortality # probability of death from NDV infection of strain 2 given chicken was never infected
p_r = 0.5 # probability of death from NDV infection of strain 1 or strain 2 given chicken was previously infected
b = ((0.75 / 30) / 15) # birth rate of new chickens in farms per susceptible chicken per day (from Table 2 of household level per month of Annapragada et al. 2019)
m_fm = ((0.0270 / 30) / 15) # migration rate of chickens from farms to markets per chicken per day (from Table 2 of household level per month of Annapragada et al. 2019)
m_mf = ((0.0176 / 30) / 15) # migration rate of chickens from markets to farms per chicken per day (from Table 2 of household level per month of Annapragada et al. 2019)
v = (1 / 90) # vaccination rate of chickens of farms per susceptible chicken of farm per day
v_hat = (1 / 60) # rate of loss of immunity due to vaccination per chicken per day
theta = (1 / 60) # rate of loss of immunity due to previous infection per chicken per day
parameters <- c(fbet1=fbet1, fbet2=fbet2, mbet1=mbet1, mbet2=mbet2,
                sig=sig, gamm=gamm, p_1, p_2, p_r, b=b, m_fm=m_fm, m_mf=m_mf,
                v=v, v_hat=v_hat, theta=theta)

# Initial conditions of the model ----------------------------------------------
init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=fI2_init, fR1=0, fR2=0, D=0, 
          fV=0, fV_E1=0, fV_I1=0, fV_E2=0, fV_I2=0,
          fV_r=0, fV_r_E1=0, fV_r_I1=0, fV_r_E2=0, fV_r_I2=0,
          fS_r=0, fE1_r=0, fE2_r=0, fI1_r=0, fI2_r=0, fR1_r=0, fR2_r=0,
          mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=mI2_init, mR1=0, mR2=0, 
          mV=0, mV_E1=0, mV_I1=0, mV_E2=0, mV_I2=0,
          mV_r=0, mV_r_E1=0, mV_r_I1=0, mV_r_E2=0, mV_r_I2=0,
          mS_r=0, mE1_r=0, mE2_r=0, mI1_r=0, mI2_r=0, mR1_r=0, mR2_r=0)

# Model equations --------------------------------------------------------------
time <- seq(0, t_max, by = t_max / (2 * length(1:t_max)))
eqn <- function(time, state, parameters){
    with(as.list(c(state, parameters)),{
      # Backyard poultry farms
      ## SEIRDV model
      dfS = b*(fS) -
            fbet1*fS*fI1 -fbet2*fS*fI2 -fbet1*fS*fI1_r -fbet2*fS*fI2_r -
            v*fS +v_hat*fV -
            m_fm*fS +m_mf*mS
      dfE1 = fbet1*fS*fI1 +fbet1*fS*fI1_r -
             sig*fE1 -
             m_fm*fE1 +m_mf*mE1
      dfE2 = fbet2*fS*fI2 +fbet2*fS*fI2_r -
             sig*fE2 -
             m_fm*fE2 +m_mf*mE2
      dfI1 = sig*fE1 -
             gamm*(1-p_1)*fI1 -gamm*p_1*fI1 -
             m_fm*fI1 +m_mf*mI1
      dfI2 = sig*fE2 -
             gamm*(1-p_2)*fI2 -gamm*p_2*fI2 -
             m_fm*fI2 + m_mf*mI2
      dfR1 = gamm*(1-p_1)*fI1 -
             theta*fR1 -
             m_fm*fR1 +m_mf*mR1
      dfR2 = gamm*(1-p_2)*fI2 -
             theta*fR2 -
             m_fm*fR2 +m_mf*mR2
      dD = gamm*p_1*fI1 +gamm*p_2*fI2 +
           gamm*p_r*fR1_r +gamm*p_r*fR2_r +
           gamm*p_1*mI1 +gamm*p_2*mI2 +
           gamm*p_r*fR2_r +gamm*p_r*mR2_r
      dfV = v*(fS) -v_hat*fV -
            fbet1*fV*fI1 -fbet2*fV*fI2 -fbet1*fV*fI1_r -fbet2*fV*fI2_r +
            gamm*fV_I1 +gamm*fV_I2 -
            m_fm*fV + m_mf*mV 
      dfV_E1 = fbet1*fV*fI1 +fbet1*fV*fI1_r -
               sig*fV_E1 -
               m_fm*fV_E1 +m_mf*mV_E1
      dfV_I1 = sig*fV_E1 -
               gamm*fV_I1 -
               m_fm*fV_I1 +m_mf*mV_I1
      dfV_E2 = fbet2*fV*fI2 +fbet2*fV*fI2_r -
               sig*fV_E2 -
               m_fm*fV_E2 +m_mf*mV_E2
      dfV_I2 = sig*fV_E2 -
               gamm*fV_I2 -
               m_fm*fV_I2 +m_mf*mV_I2
      ## SEIRDV model for recovered
      dfS_r = theta*fR1 +theta*fR2 -
              fbet1*fS_r*fI1 -fbet2*fS_r*fI2 -fbet1*fS_r*fI1_r -fbet2*fS_r*fI2_r -
              v*fS_r +v_hat*fV_r -
              m_fm*fS_r +m_mf*mS_r
      dfE1_r = fbet1*fS_r*fI1 + fbet1*fS_r*fI1_r -
               sig*fE1_r -
               m_fm*fE1_r +m_mf*mE1_r
      dfE2_r = fbet2*fS_r*fI2 + fbet2*fS_r*fI2_r -
               sig*fE2_r -
               m_fm*fE2_r +m_mf*mE2_r
      dfI1_r = sig*fE1_r -
               gamm*(1-p_r)*fI1_r -gamm*p_r*fI1_r -
               m_fm*fI1_r +m_mf*mI1_r
      dfI2_r = sig*fE2_r -
               gamm*(1-p_r)*fI2_r -gamm*p_r*fI2_r -
               m_fm*fI2_r +m_mf*mI2_r
      dfR1_r = gamm*(1-p_r)*fI1_r -
               theta*fR1_r -
               m_fm*fR1_r +m_mf*mR1_r
      dfR2_r = gamm*(1-p_r)*fI2_r -
               theta*fR2_r -
               m_fm*fR2_r +m_mf*mR2_r
      dfV_r = v*(fS_r) -v_hat*fV_r -
              fbet1*fV_r*fI1 -fbet2*fV_r*fI2 -fbet1*fV_r*fI1_r -fbet2*fV_r*fI2_r +
              gamm*fV_r_I1 +gamm*fV_r_I2 -
              m_fm*fV_r +m_mf*mV_r
      dfV_r_E1 = fbet1*fV_r*fI1 +fbet1*fV_r*fI1_r -
                 sig*fV_r_E1 -
                 m_fm*fV_r_E1 +m_mf*mV_r_E1
      dfV_r_I1 = sig*fV_r_E1 -
                 gamm*fV_r_I1 -
                 m_fm*fV_r_I1 +m_mf*mV_r_I1
      dfV_r_E2 = fbet2*fV_r*fI2 +fbet2*fV_r*fI2_r -
                 sig*fV_r_E2 -
                 m_fm*fV_r_E2 +m_mf*mV_r_E2
      dfV_r_I2 = sig*fV_r_E2 -
                 gamm*fV_r_I2 -
                 m_fm*fV_r_I2 +m_mf*mV_r_I2

      # Live bird markets
      ## SEIRDV model
      dmS = -mbet1*mS*mI1 -mbet2*mS*mI2 -mbet1*mS*mI1_r -mbet2*mS*mI2_r -
             m_mf*mS +m_fm*fS
      dmE1 = mbet1*mS*mI1 +mbet1*mS*mI1_r -
             sig*mE1 -
             m_mf*mE1 +m_fm*fE1
      dmE2 = mbet2*mS*mI2 +mbet2*mS*mI2_r -
             sig*mE2 -
             m_mf*mE2 +m_fm*fE2 
      dmI1 = sig*mE1 -
             gamm*(1-p_1)*mI1 -gamm*p_1*mI1 -
             m_mf*mI1 +m_fm*fI1
      dmI2 = sig*mE2 -
             gamm*(1-p_2)*mI2 -gamm*p_2*mI2 -
             m_mf*mI2 +m_fm*fI2
      dmR1 = gamm*(1-p_1)*mI1 -
             theta*mR1 -
             m_mf*mR1 +m_fm*fR1
      dmR2 = gamm*(1-p_2)*mI2 -
             theta*mR2 -
             m_mf*mR2 +m_fm*fR2
      dmV = -v_hat*mV -
            mbet1*mV*mI1 -mbet1*mV*mI1_r -mbet2*mV*mI2 -mbet2*mV*mI2_r +
            gamm*mV_I1 + gamm*mV_I2 +
            m_fm*fV -m_mf*mV
      dmV_E1 = mbet1*mV*mI1 +mbet1*mV*mI1_r -
               sig*mV_E1 +
               m_fm*fV_E1 -m_mf*mV_E1
      dmV_I1 = sig*mV_E1 -
               gamm*mV_I1 +
               m_fm*fV_I1 -m_mf*mV_I1
      dmV_E2 = mbet2*mV*fI2 +mbet2*mV*fI2_r -
               sig*mV_E2 +
               m_fm*fV_E2 -m_mf*mV_E2
      dmV_I2 = sig*mV_E2 -
               gamm*mV_I2 +
               m_fm*fV_I2 -m_mf*mV_I2
      ## SEIRDV model for recovered
      dmS_r = theta*mR1 +theta*mR2 -
              mbet1*mS_r*mI1 -mbet2*mS_r*mI2 -mbet1*mS_r*mI1_r -mbet2*mS_r*mI2_r +
              v_hat*mV_r +
              m_fm*fS_r -m_mf*mS_r
      dmE1_r = mbet1*mS_r*mI1 + mbet1*mS_r*mI1_r -
               sig*mE1_r -
               m_mf*mE1_r +m_fm*fE1_r
      dmE2_r = mbet2*mS_r*mI2 + mbet2*mS_r*mI2_r -
               sig*mE2_r -
               m_mf*mE2_r +m_fm*fE2_r
      dmI1_r = sig*mE1_r -
               gamm*(1-p_r)*mI1_r -gamm*p_r*mI1_r -
               m_mf*mI1_r +m_fm*fI1_r
      dmI2_r = sig*mE2_r -
               gamm*(1-p_r)*mI2_r -gamm*p_r*mI2_r -
               m_mf*mI2_r +m_fm*fI2_r
      dmR1_r = gamm*(1-p_r)*mI1_r -
               theta*mR1_r -
               m_mf*mR1_r +m_fm*fR1_r
      dmR2_r = gamm*(1-p_r)*mI2_r -
               theta*mR2_r -
               m_mf*mR2_r +m_fm*fR2_r
      dmV_r = -v_hat*mV_r -
              mbet1*mV_r*mI1 -mbet1*mV_r*mI1_r -mbet2*mV_r*mI2 -mbet2*mV_r*mI2_r +
              gamm*mV_I1 +gamm*mV_I2 +
              m_fm*fV_r -m_mf*mV_r 
      dmV_r_E1 = mbet1*mV_r*mI1 +mbet1*mV_r*mI1_r -
                 sig*mV_r_E1 +
                 m_fm*fV_r_E1 -m_mf*mV_r_E1
      dmV_r_I1 = sig*mV_r_E1 -
                 gamm*mV_r_I1 +
                 m_fm*fV_r_I1 -m_mf*mV_r_I1
      dmV_r_E2 = mbet2*mV_r*fI2 +mbet2*mV_r*fI2_r -
                 sig*mV_r_E2 +
                 m_fm*fV_r_E2 -m_mf*mV_r_E2
      dmV_r_I2 = sig*mV_r_E2 -
                 gamm*mV_r_I2 +
                 m_fm*fV_r_I2 -m_mf*mV_r_I2
      return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR1, dfR2, dD, 
                    dfV, dfV_E1, dfV_I1, dfV_E2, dfV_I2, 
                    dfV_r, dfV_r_E1, dfV_r_I1, dfV_r_E2, dfV_r_I2, 
                    dfS_r, dfE1_r, dfE2_r, dfI1_r, dfI2_r, dfR1_r, dfR2_r,
                    dmS, dmE1, dmE2, dmI1, dmI2, dmR1, dmR2, 
                    dmV, dmV_E1, dmV_I1, dmV_E2, dmV_I2, 
                    dmV_r, dmV_r_E1, dmV_r_I1, dmV_r_E2, dmV_r_I2, 
                    dmS_r, dmE1_r, dmE2_r, dmI1_r, dmI2_r, dmR1_r, dmR2_r)))})}
out <- ode(y=init, times=time, eqn, parms=parameters)
out.df <- as.data.frame(out)

# Plotting ---------------------------------------------------------------------
# Population size
plot(out.df$time, rowSums(out.df[-c(1)]), main='Population size', ylab='N', xlab='days', type='l')

# Susceptible pool divided into markets and farms
plot(out.df$time, out.df$fS, col='black', type='l', ylim=c(0, max(out.df$fS, out.df$mS)), main='Susceptible divided into markets and farms', ylab='S', xlab='days')
lines(out.df$time, out.df$mS, col='red')
legend('topright', legend=c("farms", "markets"),
       col=c("black", "red"), lty=1, cex=0.8)

# Infectious pool + Infectious and previously recovered pool divided into markets and farms
plot(out.df$time, out.df$fI1, type='l', ylim=c(0, max(out.df$fI1, out.df$mI1, out.df$fI2, out.df$mI2)),
     main='Incubation and Infectious divided into markets and farms', ylab='I + I_r', xlab='days', col='orange', lty=1)
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

# Vaccination divided into markets and farms
plot(out.df$time, out.df$mV + out.df$mV_r, col='red', type='l', ylim=c(0, max(out.df$mV + out.df$mV_r, out.df$fV + out.df$fV_r)), main='Vaccination divided into markets and farms', ylab='V', xlab='days')
lines(out.df$time, out.df$fV + out.df$fV_r, col='black')
legend('topright', legend=c("farms", "markets"),
       col=c("black", "red"), lty=1, cex=0.8)

