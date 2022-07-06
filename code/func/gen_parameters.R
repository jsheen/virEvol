# Libraries
library(deSolve)
library(foreach)
library(doParallel)
library(pracma)
library(plot.matrix)

# Time of simulation to first and second equilibria (days)
t_max_eq1 = t_max_eq2 = 2e3

# Population size
pop_size = 1e6

# Transition rate of infectiousness per chicken per day
sig = 1 / 5 

# Transition rate of recovery per chicken per day
gamm = 1 / 10

# Natural birth rate per chicken per day
b = 1 / 120
# b = 0.07 # Ian's parameter

# Natural mortality rate per chicken per day
nat_mort = 1 / 365
# nat_mort = 0.05 # Ian's parameter

# Mortality rate due to disease per chicken per day
mort = 1 / 4

# Threshold value for extinction
threshold_extinction = 2

# Parameters to define transmission rate from virulence level
c1 = 0.1
c2 = 0.3

# Virulence steps virus is allowed to take
vir_steps = seq(2.01, 100.01, 5)




# Specific parameters to use in later models
# fS_init = (pop_size * 1/2) - 1 # initial susceptible population in farms
# mS_init = (pop_size * 1/2) - 1 # initial susceptible population in markets
# fI1_init = 1 # initial strain 1 infectious population in farms
# mI1_init = 1 # initial strain 1 infectious population in markets
# perc_sold_per_farm = as.numeric(args[1]) / 100 # percent sold in interval
# inter_sell_time_per_farm = 120 # days between successive sales of chickens of a farm
# m_fm = perc_sold_per_farm / inter_sell_time_per_farm # migration rate of chickens from farms to markets per chicken per day, if unvaccinated
# m_fm_vax = (perc_sold_per_farm * 0.1) / inter_sell_time_per_farm # migration rate of chickens from farms to markets per chicken, if vaccinated
# m_mf = 1 / 7 # migration rate of chickens from markets to farms per chicken per day
# perc_vax = as.numeric(args[2]) / 100 # percent vaccinated at each campaign
# inter_vax_time = 120 # time that perc_vax is vaccinated
# v = perc_vax / inter_vax_time # vaccination rate of chickens of farms per susceptible chicken of farm per day
# v_hat = (1 / 126) # rate of loss of immunity due to vaccination per chicken per day
# # Ratio of transmission rate within markets to farms
# mfbet_ratio = 10

# Parameters
# parameters <- c(fbet1=fbet1, fbet2=fbet2, mbet1=mbet1, mbet2=mbet2,
#                 sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, 
#                 m_fm=m_fm, m_fm_vax=m_fm_vax, m_mf=m_mf, 
#                 v=v, v_hat=v_hat, theta=theta)
# # Initial conditions
# init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR1=0, fR2=0, 
#           fV=0, fV_E1=0, fV_I1=0, fV_E2=0, fV_I2=0,
#           mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=0, mR1=0, mR2=0, 
#           mV=0, mV_E1=0, mV_I1=0, mV_E2=0, mV_I2=0)

# 
# # Read arguments from command line ---------------------------------------------
# args <- commandArgs(trailingOnly=TRUE)
# print(args)
