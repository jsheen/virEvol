# Plotting function ------------------------------------------------------------
plot.out.df.migrate <- function(out.df) {
  # Global parameters
  unvax_farm_patch <- out.df$fS + out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2 + out.df$fR1 + out.df$fR2
  vax_farm_patch <- out.df$fV + out.df$fV_E1 + out.df$fV_E2 + out.df$fV_I1 + out.df$fV_I2
  unvax_market_patch <- out.df$mS + out.df$mE1 + out.df$mI1 + out.df$mE2 + out.df$mI2 + out.df$mR1 + out.df$mR2
  vax_market_patch <- out.df$mV + out.df$mV_E1 + out.df$mV_E2 + out.df$mV_I1 + out.df$mV_I2
  total_infectious <- out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2 + out.df$mE1 + out.df$mI1 + out.df$mE2 + out.df$mI2 + 
    out.df$fV_E1 + out.df$fV_E2 + out.df$fV_I1 + out.df$fV_I2 + out.df$mV_E1 + out.df$mV_E2 + out.df$mV_I1 + out.df$mV_I2
  
  # Population size
  plot(out.df$time, rowSums(out.df[-c(1)]), main='Population size', ylab='N', xlab='days', type='l', lwd=2, ylim=c(0, max(rowSums(out.df[-c(1)]))))
  lines(out.df$time, unvax_farm_patch, type='l', col='black')
  lines(out.df$time, vax_farm_patch, type='l', col='black', lty='dashed')
  lines(out.df$time, unvax_market_patch, type='l', col='red')
  lines(out.df$time, vax_market_patch, type='l', col='red', lty='dashed')
  legend('topright', legend=c("farms", "markets", 'farms (vax)', 'markets (vax)'),
         col=c("black", "red", 'black', 'red'), lty=c(1, 1, 2, 2), cex=0.8)
  
  print('At the last time step (pseudo-equilibrium), what % of total population is in...:')
  print(paste0('...unvaccinated, farm patch?: ', round(100 * (unvax_farm_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print(paste0('...vaccinated, farm patch?: ', round(100 * (vax_farm_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print(paste0('...unvaccinated, market patch?: ', round(100 * (unvax_market_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print(paste0('...vaccinated, market patch?: ', round(100 * (vax_market_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]),digits=2)))
  print('--------------------------------------------------------------')
  print('Of N in each patch, how infectious is each patch? (including incubation):')
  print(paste0('Farm patch: ', round(100 * ((out.df$fE1 + out.df$fI1 + out.df$fV_E1 + out.df$fV_I1 + out.df$fE2 + out.df$fI2 + out.df$fV_E2 + out.df$fV_I2)[nrow(out.df)] / (unvax_farm_patch[nrow(out.df)] + vax_farm_patch[nrow(out.df)])), digits=2)))
  print(paste0('Market patch: ', round(100 * ((out.df$mE1 + out.df$mI1 + out.df$mV_E1 + out.df$mV_I1 + out.df$mE2 + out.df$mI2 + out.df$mV_E2 + out.df$mV_I2)[nrow(out.df)] / (unvax_market_patch[nrow(out.df)] + vax_market_patch[nrow(out.df)])), digits=2)))
  print(paste0('Unvaccinated, farm patch: Resident strain = ', round(100 * ((out.df$fE1 + out.df$fI1)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$fE2 + out.df$fI2)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, farm patch: Resident strain = ', round(100 * ((out.df$fV_E1 + out.df$fV_I1)[nrow(out.df)] / vax_farm_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$fV_E2 + out.df$fV_I2)[nrow(out.df)] / vax_farm_patch[nrow(out.df)]), digits=2)))
  print(paste0('Unvaccinated, market patch: Resident strain = ', round(100 * ((out.df$mE1 + out.df$mI1)[nrow(out.df)] / unvax_market_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$mE2 + out.df$mI2)[nrow(out.df)] / unvax_market_patch[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, market patch: Resident strain = ', round(100 * ((out.df$mV_E1 + out.df$mV_I1)[nrow(out.df)] / vax_market_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$mV_E2 + out.df$mV_I2)[nrow(out.df)] / vax_market_patch[nrow(out.df)]), digits=2)))
  print('--------------------------------------------------------------')
  print('Of total N of infectious chickens, in which patch are they located? (including incubation):')
  print(paste0('Unvaccinated, farm patch: ', round(100 * ((out.df$fE1[nrow(out.df)] + out.df$fE2[nrow(out.df)] + out.df$fI1[nrow(out.df)] + out.df$fI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, farm patch: ', round(100 * ((out.df$fV_E1[nrow(out.df)] + out.df$fV_E2[nrow(out.df)] + out.df$fV_I1[nrow(out.df)] + out.df$fV_I2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  print(paste0('Unvaccinated, market patch: ', round(100 * ((out.df$mE1[nrow(out.df)] + out.df$mE2[nrow(out.df)] + out.df$mI1[nrow(out.df)] + out.df$mI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, market patch: ', round(100 * ((out.df$mV_E1[nrow(out.df)] + out.df$mV_E2[nrow(out.df)] + out.df$mV_I1[nrow(out.df)] + out.df$mV_I2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  
  # Susceptible pool divided into markets and farms
  plot(out.df$time, out.df$fS, col='black', type='l', ylim=c(0, max(out.df$fS, out.df$mS)), main='Susceptible divided into markets and farms', ylab='S', xlab='days')
  lines(out.df$time, out.df$mS, col='red')
  legend('topright', legend=c("farms", "markets"),
         col=c("black", "red"), lty=1, cex=0.8)
  
  # Infectious pool divided into markets and farms (excluding incubation)
  plot(out.df$time, out.df$fI1, type='l', ylim=c(0, max(out.df$fI1, out.df$mI1, out.df$fI2, out.df$mI2)),
       main='Incubation and Infectious divided into markets and farms', ylab='I', xlab='days', col='orange', lty=1)
  lines(out.df$time, out.df$mI1, col='orange', lty=2)
  lines(out.df$time, out.df$fI2, col='blue', lty=1)
  lines(out.df$time, out.df$mI2, col='blue', lty=2)
  legend('topright', legend=c("farms_strain1", "markets_strain1", "farms_strain2", "markets_strain2"),
         col=c("orange", "orange", "blue", "blue"), lty=c(1,2,1,2), cex=0.8)
  
  # Comparing strains 1 and 2 (excluding incubation)
  plot(out.df$time, out.df$mI1 + out.df$fI1, col='orange', type='l', ylim=c(0, max((out.df$fI1 + out.df$mI1), (out.df$fI2 + out.df$mI2))),
       main='Strain 1 vs. Strain 2', ylab='I', xlab='days')
  lines(out.df$time, out.df$mI2 + out.df$fI2, col='blue')
  legend('topright', legend=c("strain 1", "strain 2"),
         col=c("orange", "blue"), lty=1, cex=0.8)
  
  # Plot recovered pool
  plot(out.df$time, out.df$mR1 + out.df$fR1, col='orange', type='l', ylim=c(0, max((out.df$fR1 + out.df$mR1), (out.df$fR2 + out.df$mR2))),
       main='Strain 1 vs. Strain 2', ylab='R', xlab='days')
  lines(out.df$time, out.df$mR2 + out.df$fR2, col='blue')
  legend('topright', legend=c("strain 1", "strain 2"),
         col=c("orange", "blue"), lty=1, cex=0.8)
  
  # Vaccination divided into markets and farms
  plot(out.df$time, out.df$mV, col='red', type='l', ylim=c(0, max(out.df$mV, out.df$fV,
                                                                  out.df$mV_I1 + out.df$mV_I2,
                                                                  out.df$fV_I1 + out.df$fV_I2)), main='Vaccination divided into markets and farms', ylab='V', xlab='days')
  lines(out.df$time, out.df$fV, col='black')
  lines(out.df$time, out.df$mV_I1 + out.df$mV_I2, col='red', lty=2)
  lines(out.df$time, out.df$fV_I1 + out.df$fV_I2, col='black', lty=2)
  legend('topright', legend=c("farms", "markets", 'farms (vaccinated, infectious)', 'markets (vaccinated, infectious)'),
         col=c("black", "red", 'black', 'red'), lty=c(1, 1, 2, 2), cex=0.8)
  
  # Plot market migration flow
  # plot(migration_df$time, migration_df$market_to_farm, ylim=c(0, max(migration_df$market_to_farm, migration_df$farm_to_market)), 
  #      main='Market migration', ylab='N', xlab='days', col='black', type='l')
  # lines(migration_df$time, migration_df$farm_to_market, col='red')
  # legend('topright', legend=c('Market to farm', 'Farm to market'), col=c('black', 'red'), lty=c(1,1), cex=0.8)
}

# Migration flow in equation
# Counter for migration flows
# FM = c(m_fm*state[1:7], m_fm_vax*state[8:12])
# MF = m_mf*state[13:24]
# new_row <- data.frame(matrix(ncol=3, nrow=1))
# colnames(new_row) <- c('time', 'market_to_farm', 'farm_to_market')
# new_row$time[1] <- time
# new_row$market_to_farm[1] <- sum(MF)
# new_row$farm_to_market[1] <- sum(FM)
# migration_ls[[migration_ls_dex]] <<- new_row
# migration_ls_dex <<- migration_ls_dex + 1

# Put the following before and after out_eq1 if want to see migration flow
migration_df <- do.call(rbind, migration_ls)
migration_ls <- list()
migration_ls_dex <- 1
