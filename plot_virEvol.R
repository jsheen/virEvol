# ------------------------------------------------------------------------------
# @author: Justin Sheen
# @description: a script to plot pairwise invasibility plots
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)

pdf('~/virEvol/output.pdf')
par(mfrow = c(5, 5))
for (perc_sold_per_farm in seq(0, 0.8, by=0.2)) {
  for (perc_vax in seq(0, 0.8, by=0.2)) {
    pip <- read.csv(paste0('~/virEvol/res/', perc_sold_per_farm, '_', perc_vax, '.csv'))
    pip <- pip[,c(2:ncol(pip))]
    pip <- data.matrix(pip)
    pip <- ifelse((pip != 0 & pip != 1), NA, pip)
    colnames(pip) <- seq(0.01, 100, 5)
    rownames(pip) <- seq(0.01, 100, 5)
    melted <- melt(pip)
    ggplot(melted, aes(x = Var2, y = Var1, fill = factor(value))) + geom_tile() +
      #scale_y_reverse() +
      scale_fill_manual(values = c("black", "white"), na.value='grey')
  }
}
dev.off()

pip <- matrix(finalMatrix, ncol=length(vir_steps), nrow=length(vir_steps), byrow=F)
pip <- pracma::flipud(pip) #columns stay in place, but now from bottom to top is increasing virulence
#pip <- ifelse((pip != 0 & pip != 1), NA, pip)
#plot(pip)
write.csv(pip, paste0('~/virEvol/res/', perc_sold_per_farm, '_', perc_vax, '.csv'))
