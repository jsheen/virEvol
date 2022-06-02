# ------------------------------------------------------------------------------
# @author: Justin Sheen
# @description: a script to plot pairwise invasibility plots
# ------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(gridExtra)

l <- list()
l_dex <- 1
for (perc_sold_per_farm in seq(0, 0.8, by=0.2)) {
  for (perc_vax in seq(0, 0.8, by=0.2)) {
    pip <- read.csv(paste0('~/virEvol/res/', perc_sold_per_farm, '_', perc_vax, '_sensitivity.csv'))
    pip <- pip[,c(2:ncol(pip))]
    pip <- data.matrix(pip)
    pip <- ifelse((pip != 0 & pip != 1), NA, pip)
    colnames(pip) <- seq(0.01, 100, 5)
    rownames(pip) <- seq(95.01, 0.01, -5)
    melted <- melt(pip)
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
ggsave(filename="~/virEvol/output.tiff", marrangeGrob(grobs = l, nrow=5, ncol=5, top='Percent sold per farm in 3 months -->',
                                                     left='<-- Percent vaccinated in 3 months'),
       width=8, height=8, units='in', dpi=600)
ggsave(filename="~/virEvol/output_sensitivity.tiff", marrangeGrob(grobs = l, nrow=5, ncol=5, top='Percent sold per farm in 3 months -->',
                                                      left='<-- Percent vaccinated in 3 months'),
       width=8, height=8, units='in', dpi=600)


