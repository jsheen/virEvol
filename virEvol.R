# ------------------------------------------------------------------------------
# Space modifies predictions of virulence evolution through differential immunization
# status and contact rate in a backyard poultry farm, live market system
#
# @author: Justin Sheen
# @date: a two-patch, deterministic model that shows the dependency of virulence 
#        evolution on space due to differential immunization status and contact rate
# ------------------------------------------------------------------------------

# Assumptions ------------------------------------------------------------------
# - We assume density independent transmission rates
# - We assume that there is no transmission of the vaccine

# Load libraries ---------------------------------------------------------------
library(deSolve)

# Parameters of the model ------------------------------------------------------
# Strain 1 is a lower virulence strain
# Strain 2 is a higher virulence strain (evolved)
t = 20000 # time
fS_init = 0.5 - 1e-6 # initial susceptible population in farms
mS_init = 0.5 - 1e-6 # initial susceptible population in markets
fI1_init = 1e-6 # initial strain 1 infectious population in farms
fI2_init = 0 # initial strain 2 infectious population in farms
mI1_init = 1e-6 # initial strain 1 infectious population in markets
mI2_init = 0 # initial strain 2 infectious population in markets
if (fS_init + mS_init + fI1_init + fI2_init + mI1_init + mI2_init != 1) {
  print("Error in initial conditions")
}
fbet1 = 0 # transmission rate of strain 1 among farms (instantaneous rate)
fbet2 = 0.2 # transmission rate of strain 2 among farms (instantaneous rate)
mbet1 = 0 # transmission rate of strain 1 among markets (instantaneous rate)
mbet2 = 0.05 # transmission rate of strain 2 among markets (instantaneous rate)
sig = 1 / 5 # transition rate of infectiousness per chicken (instantaneous rate)
gamm = 1 / 12 # transition rate of recovery per chicken (instantaneous rate)
p = 0.8 # probability of death from NDV infection of strain 2 given chicken was never infected
p_r = 0.6 # probability of death from NDV infection of strain 1 or strain 2 given chicken was previously infected
b = 1 / ((0.75 / 15) / 30) # birth rate of new chickens in farms per susceptible chicken per day (from Table 2 of household level per month of Annapragada et al. 2019)
m_fm = 1 / ((0.0270 / 15) / 30) # migration rate of chickens from farms to markets per day (from Table 2 of household level per month of Annapragada et al. 2019)
m_mf = 1 / ((0.0176 / 15) / 30) # migration rate of chickens from markets to farms per day (from Table 2 of household level per month of Annapragada et al. 2019)
v =  # vaccination rate of chickens of farms per susceptible chicken of farm
v_hat = 0.006 # rate of loss of immunity due to vaccination
theta = 1 / 20 # rate of loss of immunity due to previous infection
parameters <- c(fbet1=fbet1, fbet2=fbet2, mbet1=mbet1, mbet2=mbet2,
                sig=sig, gamm=gamm, p, p_r, b=b, m_fm=m_fm, m_mf=m_mf,
                v=v, v_hat=v_hat, theta=theta)

# Initial conditions of the model ----------------------------------------------
init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=fI2_init, fR1=0, fR2=0, D=0, fV=0, fV_r=0,
          fS_r=0, fE1_r=0, fE2_r=0, fI1_r=0, fI2_r=0, fR1_r=0, fR2_r=0,
          mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=mI2_init, mR1=0, mR2=0, mV=0, mV_r=0,
          mS_r=0, mE1_r=0, mE2_r=0, mI1_r=0, mI2_r=0, mR1_r=0, mR2_r=0)

# Model equations --------------------------------------------------------------
time <- seq(0, t, by = t / (2 * length(1:t)))
eqn <- function(time,state,parameters){
    with(as.list(c(state,parameters)),{
      # Backyard poultry farms
      dfS = -fbet1*fS*fI1 -fbet2*fS*fI2 +b*(fS) -v*fS +v_hat*fV -m_fm*fS +m_mf*mS
      dfE1 = fbet1*fS*fI1 -sig*fE1 -m_fm*fE1 +m_mf*mE1
      dfE2 = fbet2*fS*fI2 -sig*fE2 -m_fm*fE2 +m_mf*mE2
      dfI1 = sig*fE1 -gamm*fI1
      dfI2 = sig*fE2 -gamm*(1-p)*fI2 - gamm*p*fI2
      dfR1 = gamm*fI1 -theta*fR1 -m_fm*fR1 +m_mf*mR1
      dfR2 = gamm*(1-p)*fI2 -theta*fR2 -m_fm*fR2 +m_mf*mR2
      dD = gamm*p*fI2 +gamm*p_r*fR2_r +gamm*p*mI2 +gamm*p_r*mR2_r
      dfV = v*(fS) -v_hat*fV -m_fm*fV + m_mf*mV
      dfS_r = theta*fR1 +theta*fR2 -fbet1*fS_r*fI1 -fbet2*fS_r*fI2 -
              fbet1*fS_r*fI1_r -fbet2*fS_r*fI2_r -
              v*fS_r +v_hat*fV_r -
              m_fm*fS_r +m_mf*mS_r
      dfE1_r = fbet1*fS_r*fI1 + fbet1*fS_r*fI1_r -sig*fE1_r -m_fm*fE1_r +m_mf*mE1_r
      dfE2_r = fbet2*fS_r*fI2 + fbet2*fS_r*fI2_r -sig*fE2_r -m_fm*fE2_r +m_mf*mE2_r
      dfI1_r = sig*fE1_r -gamm*fI1_r
      dfI2_r = sig*fE2_r -gamm*(1-p_r)*fI2_r -gamm*p_r*fI2_r
      dfR1_r = gamm*fI1_r -theta*fR1_r -m_fm*fR1_r +m_mf*mR1_r
      dfR2_r = gamm*(1-p_r)*fI2_r -theta*fR2_r -m_fm*fR2_r +m_mf*mR2_r
      dfV_r = v*(fS_r) -v_hat*fV_r -m_fm*fV_r +m_mf*mV_r
      # Live bird markets
      dmS = -mbet1*mS*mI1 -mbet2*mS*mI2 -m_mf*mS +m_fm*fS # No vaccination within live markets
      dmE1 = mbet1*mS*mI1 -sig*mE1 -m_mf*mE1 +m_fm*fE1
      dmE2 = mbet2*mS*mI2 -sig*mE2 -m_mf*mE2 +m_fm*fE2 
      dmI1 = sig*mE1 -gamm*mI1
      dmI2 = sig*mE2 -gamm*(1-p)*mI2 - gamm*p*mI2
      dmR1 = gamm*mI1 -theta*mR1 -m_mf*mR1 +m_fm*fR1
      dmR2 = gamm*(1-p)*mI2 -theta*mR2 -m_mf*mR2 +m_fm*fR2
      dmV = v*(mS) -v_hat*mV +m_fm*fV -m_mf*mV
      dmS_r = theta*mR1 +theta*mR2 -mbet1*mS_r*mI1 -mbet2*mS_r*mI2 -
        mbet1*mS_r*mI1_r -mbet2*mS_r*mI2_r -
        v*mS_r +v_hat*mV_r +
        m_fm*fS_r -m_mf*mS_r
      dmE1_r = mbet1*mS_r*mI1 + mbet1*mS_r*mI1_r -sig*mE1_r -m_mf*mE1_r +m_fm*fE1_r
      dmE2_r = mbet2*mS_r*mI2 + mbet2*mS_r*mI2_r -sig*mE2_r -m_mf*mE2_r +m_fm*fE2_r
      dmI1_r = sig*mE1_r -gamm*mI1_r
      dmI2_r = sig*mE2_r -gamm*(1-p_r)*mI2_r -gamm*p_r*mI2_r
      dmR1_r = gamm*mI1_r -theta*mR1_r -m_mf*mR1_r +m_fm*fR1_r
      dmR2_r = gamm*(1-p_r)*mI2_r -theta*mR2_r -m_mf*mR2_r +m_fm*fR2_r
      dmV_r = v*(mS_r) -v_hat*mV_r -m_mf*mV_r +m_fm*fV_r
      return(list(c(dfS, dfE1, dfE2, dfI1, dfI2, dfR1, dfR2, dD, dfV, dfV_r,
                    dfS_r, dfE1_r, dfE2_r, dfI1_r, dfI2_r, dfR1_r, dfR2_r,
                    dmS, dmE1, dmE2, dmI1, dmI2, dmR1, dmR2, dmV, dmV_r,
                    dmS_r, dmE1_r, dmE2_r, dmI1_r, dmI2_r, dmR1_r, dmR2_r)))})}
out<-ode(y=init,times=time,eqn,parms=parameters)
out.df<-as.data.frame(out)
# Susceptible pool
plot(out.df$time, out.df$fS, type='l', ylim=c(0,1))
lines(out.df$time, out.df$mS, col='red', ylim=c(0,1))
# Infectious pool
plot(out.df$time, out.df$mI1 + out.df$mI2, col='red', type='l', ylim=c(0,0.02))
lines(out.df$time, out.df$fI1 + out.df$fI2)
# Infectious pool + incubation pool
plot(out.df$time, out.df$mI1 + out.df$mI2 + out.df$mE1 + out.df$mE2, col='red', type='l', ylim=c(0,0.02))
lines(out.df$time, out.df$fI1 + out.df$fI2 + out.df$fE1 + out.df$fE2)
# Comparing strains 1 and 2
plot(out.df$time, out.df$mI1 + out.df$fI1, col='red', type='l', ylim=c(0,0.02))
lines(out.df$time, out.df$mI2 + out.df$fI2)
# Vaccination
plot(out.df$time, out.df$mV, col='red', type='l', ylim=c(0,1))
lines(out.df$time, out.df$fV)

