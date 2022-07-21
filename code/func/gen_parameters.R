# Libraries
library(deSolve)
library(foreach)
library(doParallel)
library(pracma)
library(plot.matrix)

# Time of simulation to first and second equilibria (days)
t_max_eq1 = t_max_eq2 = 4e3

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
nat_mort = 1 / 730
# nat_mort = 0.05 # Ian's parameter

# Mortality rate due to disease per chicken per day
mort = 1 / 4

# Parameters to define transmission rate from virulence level
c1 = 0.1
c2 = 0.3

# Virulence steps virus is allowed to take
vir_steps = seq(2, 100, 10)

# Plot tradeoff curve
# virulences <- seq(0, 100, 0.01)
# transmissions <- (c1 * virulences)^c2
# mortalities <- virulences / 100
# plot(mortalities, transmissions, type='l', main='Transmission-Mortality Tradeoff Curve',
#      xlab=expression(mu), ylab=expression(beta), lwd=5, col='red', cex.lab=1.2)



