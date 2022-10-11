# A script in order to visualize the equilibriums for all tests done

# Source functional scripts ----------------------------------------------------
source('~/virEvol/code/plot/plot_functions.R')
source('~/virEvol/code/func/gen_parameters.R')
source('~/virEvol/code/func/model_eqns.R')
source('~/virEvol/code/func/test_invade.R')

# Function for mod2
plot_equil <- function(res_vir, invade_vir, name) {
  # Strain specific virulence parameters
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  fbet2 <- (c1 * (invade_vir)^c2) / pop_size
  p_1 <- ((res_vir) / 100)
  p_2 <- ((invade_vir) / 100)
  
  # Parameters
  parameters <- c(fbet1=fbet1, fbet2=fbet2, sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, nat_mort=nat_mort, 
                  v=v, v_hat=v_hat)
  
  # Run resident strain until equilibrium (no vaccination)
  init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR=0, fV=0, fVE1=0, fVE2=0, fVI1=0, fVI2=0)
  time_eq1 <- seq(0, t_max_eq1, by = 1)
  out_eq1 <- ode(y=init, times=time_eq1, func=eqn, parms=parameters)
  out_eq1.df <- as.data.frame(out_eq1)
  out.df <- out_eq1.df
  # Find how much things are changing (max)
  max_change <- 0
  start_check <- which(out.df$time == (out.df$time[nrow(out.df)] - 100))
  max_tot_change <- -1
  for (col_dex in 2:ncol(out.df)) {
    #tot_change_100_days <- sum(abs(out.df[start_check:(nrow(out.df) - 1),col_dex] - out.df[(start_check + 1):nrow(out.df),col_dex]))
    tot_change_100_days <- max(out.df[start_check:nrow(out.df),col_dex]) - min(out.df[start_check:nrow(out.df),col_dex] )
    if (tot_change_100_days > max_tot_change) {
      max_tot_change <- tot_change_100_days
    }
  }
  jpeg(paste0('~/virEvol/code_output/testing/equil/', name, '.jpg'))
  plot(out.df$time, out.df$fI1, type='l', ylim=c(0, max(out.df$fI1)),
       main='', ylab='I', xlab='days', col='orange', lty=1)
  dev.off()
  return(max_tot_change)
}

# Baseline
# Set model 2 specific parameters and functions --------------------------------
# For low transmission tradeoff
c1 = 0.1
c2 = 0.3
# Initial susceptible population
fS_init = pop_size - 1
# Initial strain 1 infectious population in farms
fI1_init = 1
# Percent of susceptible chickens vaccinated in each time period
perc_vax = 0.26
# Time that perc_vax is vaccinated
inter_vax_time = 120 
# Vaccination rate of chickens of farms per susceptible chicken of farm per day
v = perc_vax / inter_vax_time 
# Rate of loss of immunity due to vaccination per chicken per day
v_hat = (1 / 120)
# Threshold value for extinction
threshold_extinction = 1
# Assign model 2 specific equation and test_invade -----------------------------
eqn <- eqn_mod2
test_invade <- test_invade_mod2
# Create and save PIP with no differential migration due to vaccination --------
m_fm_vax = m_fm
combos <- list()
combos_dex <- 1
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix <- foreach(i=vir_steps, .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  max_tot_change <- plot_equil(res_vir=i, invade_vir=i, name=paste0('baseline_mod2_tradeoff1 ', i))
  max_tot_change
}
jpeg('~/virEvol/code_output/testing/equil_plot/baseline_mod2_tradeoff1.jpg')
plot(vir_steps, finalMatrix, type='l', main='baseline_mod2_tradeoff1', ylab='tot change in past 100 days')
dev.off()