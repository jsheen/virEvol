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

# Natural mortality rate per chicken per day
nat_mort = 1 / 365

# Mortality rate due to disease per chicken per day
mort = 1 / 4

# Virulence steps virus is allowed to take
vir_steps = seq(2, 100, 1)