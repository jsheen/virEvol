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
t_max = 4e4 # time of simulation (days)
pop_size = 1e5 # population size
fS_init = (pop_size * 1/2) - 2 # initial susceptible population in farms
mS_init = (pop_size * 1/2) - 2 # initial susceptible population in markets
fI1_init = 1 # initial strain 1 infectious population in farms
mI1_init = 1 # initial strain 1 infectious population in markets
fI2_init = 1 # initial strain 2 infectious population in farms
mI2_init = 1 # initial strain 2 infectious population in markets
sig = 1 / 5 # transition rate of infectiousness per chicken per day
gamm = 1 / 12#4.5 # transition rate of recovery per chicken per day
mort = 1 / 4 # disease mortality rate per chicken per day
nat_mort = 1 / 730 # natural mortality rate per chicken per day

res_vir = 30.01
invade_vir = 80.01
fbet1 = ((((0.005 * res_vir)^0.05) / 100) + 0.5) / pop_size
mbet1 = fbet1 * mfbet_ratio
fbet2 = ((((0.005 * invade_vir)^0.05) / 100) + 0.5) / pop_size
mbet2 = fbet2 * mfbet_ratio
p_1 = ((res_vir * 0.5) / 100) + 0.5
p_2 = ((invade_vir * 0.5) / 100) + 0.5

# fbet1 = (15 / 30) / pop_size # transmission rate of strain 1 among farms per SI contact per day
# mbet1 = (40 / 7) / pop_size # transmission rate of strain 1 among markets per SI contact per day
# strain_2_multiplication_transmission = 1.0001
# strain_2_multiplication_mortality = 1.2
# fbet2 = fbet1 * strain_2_multiplication_transmission # transmission rate of strain 2 among farms per SI contact per day
# mbet2 = mbet1 * strain_2_multiplication_transmission # transmission rate of strain 2 among markets per SI contact per day
# p_1 = 0.85 # probability of death from NDV infection of strain 1 given chicken was never infected
# p_2 = p_1 * strain_2_multiplication_mortality # probability of death from NDV infection of strain 2 given chicken was never infected

b = ((0.75 / 30) / 15) * 5 # birth rate of new chickens in farms per susceptible chicken per day (from Table 2 of household level per month of Annapragada et al. 2019)
perc_sold_per_farm = 0.2 # percent sold in interval
inter_sell_time_per_farm = 365 # days between successive sales of chickens of a farm
m_fm = 0#perc_sold_per_farm / inter_sell_time_per_farm # migration rate of chickens from farms to markets per chicken per day
m_mf = (1 / 7) # migration rate of chickens from markets to farms per chicken per day
perc_vax = 0.1 # percent vaccinated at each campaign
inter_vax_time = 365 # time that perc_vax is vaccinated
v = 0#perc_vax / inter_vax_time # vaccination rate of chickens of farms per susceptible chicken of farm per day
v_hat = (1 / 126) # rate of loss of immunity due to vaccination per chicken per day
theta = (1 / 126) # rate of loss of immunity due to previous infection per chicken per day
parameters <- c(fbet1=fbet1, fbet2=fbet2, mbet1=mbet1, mbet2=mbet2,
                sig=sig, gamm=gamm, mort=mort, p_1=p_1, p_2=p_2, b=b, 
                m_fm=m_fm, m_mf=m_mf,
                v=v, v_hat=v_hat, theta=theta)

# virulences <- seq(0.01, 100, 0.1)
# morts <- ((virulences * 0.4) / 100) + 0.6
# betas <- ((((0.005 * virulences)^0.05) / 100) + 0.5)
# plot(morts, betas, type='l')

# Initial conditions of the model ----------------------------------------------
init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=fI2_init, fR1=0, fR2=0, 
          fV=0, fV_E1=0, fV_I1=0, fV_E2=0, fV_I2=0,
          mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=mI2_init, mR1=0, mR2=0, 
          mV=0, mV_E1=0, mV_I1=0, mV_E2=0, mV_I2=0)

# Model equations --------------------------------------------------------------
time <- seq(0, t_max, by = t_max / (2 * length(1:t_max)))
eqn <- function(time, state, parameters){
    with(as.list(c(state, parameters)),{
      # Backyard poultry farms
      dfS = b*(fS + fE1 + fE2 + fI1 + fI2 + fR1 + fR2) -
            fbet1*fS*fI1 -fbet2*fS*fI2 -
            fbet1*fS*fV_I1 -fbet2*fS*fV_I2 -
            v*fS +v_hat*fV -
            m_fm*fS +m_mf*mS -
            nat_mort*fS +
            theta*fR1 + theta*fR2
      dfE1 = fbet1*fS*fI1 +fbet1*fS*fV_I1 -
             sig*fE1 -
             m_fm*fE1 +m_mf*mE1 -
             nat_mort*fE1
      dfE2 = fbet2*fS*fI2 +fbet2*fS*fV_I2 -
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
            fbet1*fV*fI1 -fbet2*fV*fI2 -
            fbet1*fV*fV_I1 -fbet2*fV*fV_I2 +
            gamm*fV_I1 +gamm*fV_I2 -
            m_fm*fV + m_mf*mV -
            nat_mort*fV
      dfV_E1 = fbet1*fV*fI1 +fbet1*fV*fV_I1 -
               sig*fV_E1 -
               m_fm*fV_E1 +m_mf*mV_E1 -
               nat_mort*fV_E1
      dfV_I1 = sig*fV_E1 -
               gamm*fV_I1 -
               m_fm*fV_I1 +m_mf*mV_I1 -
               nat_mort*fV_I1
      dfV_E2 = fbet2*fV*fI2 +fbet2*fV*fV_I2 -
               sig*fV_E2 -
               m_fm*fV_E2 +m_mf*mV_E2 -
               nat_mort*fV_E2
      dfV_I2 = sig*fV_E2 -
               gamm*fV_I2 -
               m_fm*fV_I2 +m_mf*mV_I2 -
               nat_mort*fV_I2

      # Live bird markets
      dmS =  b*(mS + mE1 + mE2 + mI1 + mI2 + mR1 + mR2)-
             mbet1*mS*mI1 -mbet2*mS*mI2 -
             m_mf*mS +m_fm*fS -
             nat_mort*mS +
             theta*mR1 + theta*mR2
      dmE1 = mbet1*mS*mI1 -
             sig*mE1 -
             m_mf*mE1 +m_fm*fE1 -
             nat_mort*mE1
      dmE2 = mbet2*mS*mI2 -
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
            mbet1*mV*mI1 -mbet2*mV*mI2 +
            gamm*mV_I1 + gamm*mV_I2 +
            m_fm*fV -m_mf*mV -
            nat_mort*mV
      dmV_E1 = mbet1*mV*mI1 -
               sig*mV_E1 +
               m_fm*fV_E1 -m_mf*mV_E1 -
               nat_mort*mV_E1
      dmV_I1 = sig*mV_E1 -
               gamm*mV_I1 +
               m_fm*fV_I1 -m_mf*mV_I1 -
               nat_mort*mV_I1
      dmV_E2 = mbet2*mV*fI2 -
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

# Vaccination divided into markets and farms
plot(out.df$time, out.df$mV, col='red', type='l', ylim=c(0, max(out.df$mV, out.df$fV)), main='Vaccination divided into markets and farms', ylab='V', xlab='days')
lines(out.df$time, out.df$fV, col='black')
legend('topright', legend=c("farms", "markets"),
       col=c("black", "red"), lty=1, cex=0.8)




