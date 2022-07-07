# Model 1: SEIR ----------------------------------------------------------------
plot.out.df.mod1 <- function(out.df) {
  # Global parameters
  unvax_farm_patch <- out.df$fS + out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2 + out.df$fR
  total_infectious <- out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2
  
  # Population size
  plot(out.df$time, rowSums(out.df[-c(1)]), main='Population size', ylab='N', xlab='days', type='l', lwd=2, ylim=c(0, max(rowSums(out.df[-c(1)]))))
  lines(out.df$time, unvax_farm_patch, type='l', col='black')
  
  print('At the last time step, what % of total population is in...:')
  print(paste0('...unvaccinated, farm patch?: ', round(100 * (unvax_farm_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print('--------------------------------------------------------------')
  print('Of N in each patch, how infectious is each patch? (including incubation):')
  print(paste0('Unvaccinated, farm patch: Res strain = ', round(100 * ((out.df$fE1 + out.df$fI1)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$fE2 + out.df$fI2)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2)))
  print('--------------------------------------------------------------')
  print('Of total N of infectious chickens, in which patch are they located? (including incubation):')
  print(paste0('Unvaccinated, farm patch: ', round(100 * ((out.df$fE1[nrow(out.df)] + out.df$fE2[nrow(out.df)] + out.df$fI1[nrow(out.df)] + out.df$fI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))

  # Susceptible pool
  plot(out.df$time, out.df$fS, col='black', type='l', ylim=c(0, max(out.df$fS)), main='Susceptible', ylab='S', xlab='days')

  # Comparing strains 1 and 2 (excluding incubation)
  plot(out.df$time, out.df$fI1, col='orange', type='l', ylim=c(0, max((out.df$fI1), (out.df$fI2))),
       main='Strain 1 vs. Strain 2', ylab='I', xlab='days')
  lines(out.df$time, out.df$fI2, col='blue')
  legend('topright', legend=c("strain 1", "strain 2"),
         col=c("orange", "blue"), lty=1, cex=0.8)
  
  # Plot recovered pool
  plot(out.df$time, out.df$fR, col='orange', type='l', ylim=c(0, max(out.df$fR)), ylab='R', xlab='days')
}

# Model 2: SEIR with vaccination -----------------------------------------------
plot.out.df.mod2 <- function(out.df) {
  # Global parameters
  unvax_farm_patch <- out.df$fS + out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2 + out.df$fR
  vax_farm_patch <- out.df$fV + out.df$fVE1 + out.df$fVE2 + out.df$fVI1 + out.df$fVI2
  total_infectious <- out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2 +out.df$fVE1 + out.df$fVE2 + out.df$fVI1 + out.df$fVI2
  
  # Population size
  plot(out.df$time, rowSums(out.df[-c(1)]), main='Population size', ylab='N', xlab='days', type='l', lwd=2, ylim=c(0, max(rowSums(out.df[-c(1)]))))
  lines(out.df$time, unvax_farm_patch, type='l', col='black')
  lines(out.df$time, vax_farm_patch, type='l', col='black', lty='dashed')
  legend('topright', legend=c("farms (unvax)",'farms (vax)'),
         col=c("black", 'black'), lty=c(1, 2), cex=0.8)
  
  print('At the last time step, what % of total population is in...:')
  print(paste0('...unvaccinated, farm patch?: ', round(100 * (unvax_farm_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print(paste0('...vaccinated, farm patch?: ', round(100 * (vax_farm_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print('--------------------------------------------------------------')
  print('Of N in each patch, how infectious is each patch? (including incubation):')
  print(paste0('Unvaccinated, farm patch: Res strain = ', round(100 * ((out.df$fE1 + out.df$fI1)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$fE2 + out.df$fI2)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, farm patch: Res strain = ', round(100 * ((out.df$fVE1 + out.df$fVI1)[nrow(out.df)] / vax_farm_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$fVE2 + out.df$fVI2)[nrow(out.df)] / vax_farm_patch[nrow(out.df)]), digits=2)))
  print('--------------------------------------------------------------')
  print('Of total N of infectious chickens, in which patch are they located? (including incubation):')
  print(paste0('Unvaccinated, farm patch: ', round(100 * ((out.df$fE1[nrow(out.df)] + out.df$fE2[nrow(out.df)] + out.df$fI1[nrow(out.df)] + out.df$fI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, farm patch: ', round(100 * ((out.df$fVE1[nrow(out.df)] + out.df$fVE2[nrow(out.df)] + out.df$fVI1[nrow(out.df)] + out.df$fVI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))

  # Susceptible pool divided into markets and farms
  plot(out.df$time, out.df$fS, col='black', type='l', ylim=c(0, max(out.df$fS, out.df$mS)), main='Susceptible divided into markets and farms', ylab='S', xlab='days')
  legend('topright', legend=c("farms"),
         col=c("black"), lty=1, cex=0.8)
  
  # Comparing strains 1 and 2 (excluding incubation)
  plot(out.df$time, out.df$fI1, col='orange', type='l', ylim=c(0, max((out.df$fI1), (out.df$fI2))),
       main='Strain 1 vs. Strain 2', ylab='I', xlab='days')
  lines(out.df$time, out.df$fI2, col='blue')
  legend('topright', legend=c("strain 1", "strain 2"),
         col=c("orange", "blue"), lty=1, cex=0.8)
  
  # Plot recovered pool
  plot(out.df$time, out.df$fR, col='orange', type='l', ylim=c(0, max(out.df$fR)),
       main='', ylab='R', xlab='days')
  
  # Vaccination
  plot(out.df$time, out.df$fV, col='red', type='l', ylim=c(0, max(out.df$fV, out.df$fVI1 + out.df$fVI2)), main='Vaccination divided into markets and farms', ylab='V', xlab='days')
  lines(out.df$time, out.df$fVI1 + out.df$fVI2, col='black', lty=2)
  legend('topright', legend=c("farms", 'farms (infectious)'),
         col=c("red", "red"), lty=c(1, 1), cex=0.8)
}

# Model 3: SEIR with vaccination and migration ---------------------------------
plot.out.df.mod3 <- function(out.df) {
  # Global parameters
  unvax_farm_patch <- out.df$fS + out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2 + out.df$fR
  vax_farm_patch <- out.df$fV + out.df$fVE1 + out.df$fVE2 + out.df$fVI1 + out.df$fVI2
  unvax_market_patch <- out.df$mS + out.df$mE1 + out.df$mI1 + out.df$mE2 + out.df$mI2 + out.df$mR
  vax_market_patch <- out.df$mV + out.df$mVE1 + out.df$mVE2 + out.df$mVI1 + out.df$mVI2
  total_infectious <- out.df$fE1 + out.df$fI1 + out.df$fE2 + out.df$fI2 + out.df$mE1 + out.df$mI1 + out.df$mE2 + out.df$mI2 + 
    out.df$fVE1 + out.df$fVE2 + out.df$fVI1 + out.df$fVI2 + out.df$mVE1 + out.df$mVE2 + out.df$mVI1 + out.df$mVI2
  
  # Population size
  plot(out.df$time, rowSums(out.df[-c(1)]), main='Population size', ylab='N', xlab='days', type='l', lwd=2, ylim=c(0, max(rowSums(out.df[-c(1)]))))
  lines(out.df$time, unvax_farm_patch, type='l', col='black')
  lines(out.df$time, vax_farm_patch, type='l', col='black', lty='dashed')
  lines(out.df$time, unvax_market_patch, type='l', col='red')
  lines(out.df$time, vax_market_patch, type='l', col='red', lty='dashed')
  legend('topright', legend=c("farms", "markets", 'farms (vax)', 'markets (vax)'),
         col=c("black", "red", 'black', 'red'), lty=c(1, 1, 2, 2), cex=0.8)
  
  print('At the last time step, what % of total population is in...:')
  print(paste0('...unvaccinated, farm patch?: ', round(100 * (unvax_farm_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print(paste0('...vaccinated, farm patch?: ', round(100 * (vax_farm_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print(paste0('...unvaccinated, market patch?: ', round(100 * (unvax_market_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]), digits=2)))
  print(paste0('...vaccinated, market patch?: ', round(100 * (vax_market_patch[nrow(out.df)] / rowSums(out.df[-c(1)])[nrow(out.df)]),digits=2)))
  print('--------------------------------------------------------------')
  print('Of N in each patch, how infectious is each patch? (including incubation):')
  print(paste0('Unvaccinated, farm patch: Res strain = ', round(100 * ((out.df$fE1 + out.df$fI1)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$fE2 + out.df$fI2)[nrow(out.df)] / unvax_farm_patch[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, farm patch: Res strain = ', round(100 * ((out.df$fVE1 + out.df$fVI1)[nrow(out.df)] / vax_farm_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$fVE2 + out.df$fVI2)[nrow(out.df)] / vax_farm_patch[nrow(out.df)]), digits=2)))
  print(paste0('Unvaccinated, market patch: Res strain = ', round(100 * ((out.df$mE1 + out.df$mI1)[nrow(out.df)] / unvax_market_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$mE2 + out.df$mI2)[nrow(out.df)] / unvax_market_patch[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, market patch: Res strain = ', round(100 * ((out.df$mVE1 + out.df$mVI1)[nrow(out.df)] / vax_market_patch[nrow(out.df)]), digits=2), 
               '; Invader strain = ',  round(100 * ((out.df$mVE2 + out.df$mVI2)[nrow(out.df)] / vax_market_patch[nrow(out.df)]), digits=2)))
  print('--------------------------------------------------------------')
  print('Of total N of infectious chickens, in which patch are they located? (including incubation):')
  print(paste0('Unvaccinated, farm patch: ', round(100 * ((out.df$fE1[nrow(out.df)] + out.df$fE2[nrow(out.df)] + out.df$fI1[nrow(out.df)] + out.df$fI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, farm patch: ', round(100 * ((out.df$fVE1[nrow(out.df)] + out.df$fVE2[nrow(out.df)] + out.df$fVI1[nrow(out.df)] + out.df$fVI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  print(paste0('Unvaccinated, market patch: ', round(100 * ((out.df$mE1[nrow(out.df)] + out.df$mE2[nrow(out.df)] + out.df$mI1[nrow(out.df)] + out.df$mI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  print(paste0('Vaccinated, market patch: ', round(100 * ((out.df$mVE1[nrow(out.df)] + out.df$mVE2[nrow(out.df)] + out.df$mVI1[nrow(out.df)] + out.df$mVI2[nrow(out.df)]) / total_infectious[nrow(out.df)]), digits=2)))
  
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
  plot(out.df$time, out.df$mR + out.df$fR, col='orange', type='l', ylim=c(0, max(out.df$fR + out.df$mR)),
       main='', ylab='R', xlab='days')
  
  # Vaccination divided into markets and farms
  plot(out.df$time, out.df$mV, col='red', type='l', ylim=c(0, max(out.df$mV, out.df$fV,
                                                                  out.df$mVI1 + out.df$mVI2,
                                                                  out.df$fVI1 + out.df$fVI2)), main='Vaccination divided into markets and farms', ylab='V', xlab='days')
  lines(out.df$time, out.df$fV, col='black')
  lines(out.df$time, out.df$mVI1 + out.df$mVI2, col='red', lty=2)
  lines(out.df$time, out.df$fVI1 + out.df$fVI2, col='black', lty=2)
  legend('topright', legend=c("farms", "markets", 'farms (infectious)', 'markets (infectious)'),
         col=c("black", "red", 'black', 'red'), lty=c(1, 1, 2, 2), cex=0.8)
}
