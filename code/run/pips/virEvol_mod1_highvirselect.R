# ------------------------------------------------------------------------------
# Model 1: SEIR
# ------------------------------------------------------------------------------
# Source functional scripts ----------------------------------------------------
source('~/virEvol/code/plot/plot_functions.R')
source('~/virEvol/code/func/gen_parameters.R')
source('~/virEvol/code/func/model_eqns.R')
source('~/virEvol/code/func/test_invade.R')

# Set model 1 specific parameters and functions --------------------------------
# Initial susceptible population
fS_init = pop_size - 1
# Initial strain 1 infectious population in farms
fI1_init = 1
# Threshold value for extinction
threshold_extinction = 1.1

# Assign model 1 specific equation and test_invade -----------------------------
eqn <- eqn_mod1
test_invade <- test_invade_mod1

# Create and save PIP ----------------------------------------------------------
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
  write.csv(tempMatrix, paste0("~/virEvol/scratch/", i[1], "_", i[2], ".csv")) # for debugging purposes
  tempMatrix
}
stopCluster(cl)
pip <- matrix(finalMatrix, ncol=length(vir_steps), nrow=length(vir_steps), byrow=F)
pip <- pracma::flipud(pip)
write.csv(pip, paste0('~/virEvol/code_output/pips/main/mod1_highvirselect.csv'))

