# ------------------------------------------------------------------------------
# @author: Justin Sheen
# @description: a script to plot pairwise invasibility plots
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(gridExtra)

l <- list()
l_dex <- 1
for (perc_sold_per_farm in seq(0, 0.33, by=0.33)) {
  for (perc_vax in seq(0, 0.33, by=0.33)) {
    pip <- read.csv(paste0('~/virEvol/res/', perc_sold_per_farm, '_', perc_vax, '_diff.csv'))
    pip <- pip[,c(2:ncol(pip))]
    pip <- data.matrix(pip)
    colnames(pip) <- seq(0.01, 99.01, 10)
    rownames(pip) <- rev(seq(0.01, 99.01, 10))
    pip_toPlot <- ifelse((pip == 3), 1, pip)
    pip_toPlot <- ifelse((pip_toPlot == 5), NA, pip_toPlot)
    pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot)
    if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.',
                                                                             ' %sold:', perc_sold_per_farm,
                                                                             ' %vax:', perc_vax)) }
    melted <- melt(pip_toPlot)
    colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      #scale_y_reverse() +
      scale_fill_manual(values = c("black", "white"), na.value='grey') + 
      geom_abline(intercept = 0, slope = 1, col='red') +
      theme(legend.position="none")
    l[[l_dex]] <- temp_plot
    l_dex <- l_dex + 1
  }
}
ggsave(filename="~/virEvol/output_diff.tiff", marrangeGrob(grobs = l, nrow=3, ncol=3, top='Percent sold per farm in 3 months -->',
                                                     left='<-- Percent vaccinated in 3 months'),
       width=6, height=6, units='in', dpi=600)


