# Libraries
library(deSolve)
library(foreach)
library(doParallel)
library(pracma)
library(plot.matrix)

# Time of simulation to first equilibrium (days)
t_max_eq1 = 4e3

# Time of simulation to second equilibrium (days)
t_max_eq2 = 4e3 

# Population size
pop_size = 1e6

# Transition rate of infectiousness per chicken per day
sig = 1 / 5 

# Transition rate of recovery per chicken per day
gamm = 1 / 10

# Natural birth rate per chicken per day
b = 1 / 120

# Natural mortality rate per chicken per day
nat_mort = 1 / 730

# Mortality rate due to disease per chicken per day
mort = 1 / 4

# Threshold value for extinction
threshold_extinction = 2

# Parameters to define transmission rate from virulence level
c1 = 0.1
c2 = 0.3

# Virulence steps virus is allowed to take
vir_steps = seq(2.01, 100.01, 5)





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



