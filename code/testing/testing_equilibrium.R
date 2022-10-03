# A script in order to visualize the equilibriums for all tests done

# Source functional scripts ----------------------------------------------------
source('~/virEvol/code/plot/plot_functions.R')
source('~/virEvol/code/func/gen_parameters.R')
source('~/virEvol/code/func/model_eqns.R')
source('~/virEvol/code/func/test_invade.R')

# Function for mod4
plot_equil <- function(res_vir, invade_vir, name) {
  # Strain specific virulence parameters
  fbet1 <- (c1 * (res_vir)^c2) / pop_size
  mbet1 <- fbet1 * bet_mf_ratio
  fbet2 <- (c1 * (invade_vir)^c2) / pop_size
  mbet2 <- fbet2 * bet_mf_ratio
  p_1 <- ((res_vir) / 100)
  p_2 <- ((invade_vir) / 100)
  # Parameters
  parameters <- c(fbet1=fbet1, fbet2=fbet2, sig=sig, gamm=gamm, p_1=p_1, p_2=p_2, mort=mort, b=b, nat_mort=nat_mort, 
                  v=v, v_hat=v_hat, m_fm=m_fm, m_fm_vax=m_fm_vax, m_mf=m_mf, mbet1=mbet1, mbet2=mbet2, p_s=p_s)
  # Run resident strain until equilibrium (no vaccination)
  init <- c(fS=fS_init, fE1=0, fE2=0, fI1=fI1_init, fI2=0, fR=0, fV=0, fVE1=0, fVE2=0, fVI1=0, fVI2=0,
            mS=mS_init, mE1=0, mE2=0, mI1=mI1_init, mI2=0, mR=0, mV=0, mVE1=0, mVE2=0, mVI1=0, mVI2=0,
            M_FM=0, M_MF=0)
  time_eq1 <- seq(0, t_max_eq1, by = 1)
  out_eq1 <- ode(y=init, times=time_eq1, eqn, parms=parameters)
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
  plot(out.df$time, out.df$fI1, type='l', ylim=c(0, max(out.df$fI1, out.df$mI1, out.df$fI2, out.df$mI2)),
       main='Incubation and Infectious divided into markets and farms', ylab='I', xlab='days', col='orange', lty=1)
  lines(out.df$time, out.df$mI1, col='orange', lty=2)
  lines(out.df$time, out.df$fI2, col='blue', lty=1)
  lines(out.df$time, out.df$mI2, col='blue', lty=2)
  legend('topright', legend=c("farms_strain1", "markets_strain1", "farms_strain2", "markets_strain2", paste0('change: ', max_tot_change)),
         col=c("orange", "orange", "blue", "blue", 'black'), lty=c(1,2,1,2,1), cex=0.8)
  dev.off()
  return(max_tot_change)
}

# Baseline (high virulence)
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
bet_mf_ratio = 5
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
cores <- detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)
finalMatrix <- foreach(i=vir_steps, .combine=cbind) %dopar% {
  library(deSolve)
  library(foreach)
  library(doParallel)
  max_tot_change <- plot_equil(res_vir=i, invade_vir=i, name=paste0('baseline ', i))
  max_tot_change
}
jpeg('~/virEvol/code_output/testing/equil_plot/baseline.jpg')
plot(vir_steps, finalMatrix, type='l', main='baseline', ylab='tot change in past 100 days')
dev.off()