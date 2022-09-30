# ------------------------------------------------------------------------------
# Model 4: SEIR with vaccination and migration and slaughter
# ------------------------------------------------------------------------------
# Source functional scripts ----------------------------------------------------
source('~/virEvol/code/plot/plot_functions.R')
source('~/virEvol/code/func/gen_parameters.R')
source('~/virEvol/code/func/model_eqns.R')
source('~/virEvol/code/func/test_invade.R')

# Set model 4 specific parameters and functions --------------------------------
c1 = 1
c2 = 0.4
# Initial susceptible population in farms
fS_init = (pop_size * 1/2) - 1
# Initial susceptible population in markets
mS_init = (pop_size * 1/2) - 1
# Initial strain 1 infectious population in farms
fI1_init = 1
# Initial strain 1 infectious population in markets
mI1_init = 1
# Percent of susceptible chickens vaccinated in each time period
perc_vax = 0.33 
# Time that perc_vax is vaccinated
inter_vax_time = 120 
# Vaccination rate of chickens of farms per susceptible chicken of farm per day
v = perc_vax / inter_vax_time 
# Rate of loss of immunity due to vaccination per chicken per day
v_hat = 1 / 120
# Percent sold in interval
perc_sold_per_farm = 0.33
# Days between successive sales of chickens of a farm
inter_sell_time_per_farm = 120 
# Migration rate of chickens from farms to markets per chicken per day, if unvaccinated
m_fm = perc_sold_per_farm / inter_sell_time_per_farm
# Migration rate of chickens from markets to farms per chicken per day
m_mf = 1 / 7
# Ratio of contact rate in markets vs. farms
bet_mf_ratio = 10
# Threshold value for extinction
threshold_extinction = 2.2
# Percentage of market chickens that are to be immediately slaughtered
p_s = 0.8

# Assign model 4 specific equation and test_invade -----------------------------
eqn <- eqn_mod4
test_invade <- test_invade_mod4

# Create and save PIP with no differential migration due to vaccination --------
m_fm_vax = m_fm
combos <- list()
combos_dex <- 1
for (i in vir_steps) {
  for (j in vir_steps) {
    combos[[combos_dex]] <- c(i, j)
    combos_dex <- combos_dex + 1
  }
}
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix <- foreach(i=combos, .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  tempMatrix = test_invade(res_vir=i[1], invade_vir=i[2])
  # write.csv(tempMatrix, paste0("~/virEvol/scratch/", i[1], "_", i[2], ".csv")) # for debugging purposes
  tempMatrix
}
stopCluster(cl)
pip <- matrix(finalMatrix, ncol=length(vir_steps), nrow=length(vir_steps), byrow=F)
pip <- pracma::flipud(pip)
write.csv(pip, paste0('~/virEvol/code_output/pips/main/mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens4.csv'))
