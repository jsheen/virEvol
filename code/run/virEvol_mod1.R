# ------------------------------------------------------------------------------
# Model 1: SEIR model
# ------------------------------------------------------------------------------

# Read arguments from command line ---------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
print(args)

# Source functional scripts ----------------------------------------------------
source('~/virEvol/code/func/gen_parameters.R')
source('~/virEvol/code/func/model_eqns.R')
source('~/virEvol/code/func/test_invade.R')

# Set model 1 specific parameters and functions --------------------------------
# Initial susceptible population
fS_init = pop_size - 1
# Initial strain 1 infectious population in farms
fI1_init = 1
# Parameters
parameters <- c(fbet1=fbet1, fbet2=fbet2, sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b)
# Initial conditions
init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR=0)

# Assign model 1 to equation ---------------------------------------------------
eqn <- eqn_mod1

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
  tempMatrix = test_invade(res_vir=i[1], invade_vir=i[2], markets=FALSE)
  # write.csv(tempMatrix, paste0("~/virEvol/scratch/", i[1], "_", i[2], ".csv")) # for debugging purposes
  tempMatrix
}
stopCluster(cl)
pip <- matrix(finalMatrix, ncol=length(vir_steps), nrow=length(vir_steps), byrow=F)
pip <- pracma::flipud(pip)
write.csv(pip, paste0('~/virEvol/code_output/pips/', perc_sold_per_farm, '_', perc_vax, '_mod1.csv'))




