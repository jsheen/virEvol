# ------------------------------------------------------------------------------
# Plotting functions
# ------------------------------------------------------------------------------
# Load libraries ---------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

# Load global variables --------------------------------------------------------
vir_steps <- seq(2, 100, 1)

# Function to find singular strategy -------------------------------------------
find_singular_strat <- function(name, vir_steps) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
  singular_strats <- c()
  for (col_dex in 1:ncol(pip)) {
    if (length(which(pip[,col_dex] == 0)) == nrow(pip)) {
      singular_strats <- c(singular_strats, col_dex)
    }
  }
  col_singular_strat <- singular_strats[floor(length(singular_strats) / 2)]
  vir_singular_strat <- vir_steps[col_singular_strat]
  print(paste0(name, ": ", vir_singular_strat))
}
find_singular_strat('mod1', vir_steps=vir_steps)
find_singular_strat('mod2', vir_steps=vir_steps)
find_singular_strat('mod3_nodiff', vir_steps=vir_steps)
find_singular_strat('mod3_diff', vir_steps=vir_steps)
find_singular_strat('mod3_nodiff_mpatch', vir_steps=vir_steps)
find_singular_strat('mod3_diff_mpatch', vir_steps=vir_steps)

# Function to plot each of the results -----------------------------------------
plot_pip <- function(name, vir_steps, title) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray"), na.value='white') +geom_abline(intercept = 0, slope = 1, col='red') +theme(legend.position="none") +
    xlab('Resident virulence') + ylab('Invader virulence') + ggtitle(title)
  l <- list()
  l[[1]] <- temp_plot
  ggsave(filename=paste0("~/virEvol/code_output/plots/virEvol_", name, ".png"), marrangeGrob(grobs = l, nrow=1, ncol=1, top=NULL), width=3, height=3, units='in', dpi=600)
}

# Plot each model result -------------------------------------------------------
plot_pip('mod1', vir_steps=vir_steps, 'Without vaccination\n(Model 1)')
plot_pip('mod2', vir_steps=vir_steps, 'With vaccination\n(Model 2)')
plot_pip('mod3_nodiff', vir_steps=vir_steps, 'With vaccination and\nnon-differential migration\n(Model 3)')
plot_pip('mod3_diff', vir_steps=vir_steps, 'With vaccination and\ndifferential migration\n(Model 3)')
plot_pip('mod3_nodiff_mpatch', vir_steps=vir_steps, 'With vaccination and\nnon-differential migration\n(Model 3 mpatch)')
plot_pip('mod3_diff_mpatch', vir_steps=vir_steps, 'With vaccination and\ndifferential migration\n(Model 3 mpatch)')

# Plot figure 1 of paper -------------------------------------------------------
l <- list()
l_dex <- 1
titles <- c("(A) Without vaccination\n", "(B) With vaccination\n", 
            "(C) With vaccination &\nnon-differential migration", 
            "(D) With vaccination &\ndifferential migration")
for (name in c('mod1', 'mod2', 'mod3_nodiff_mpatch', 'mod3_diff_mpatch')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  if (l_dex == 1) {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('Resident virulence') + ylab('Invader virulence') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + 
      theme(plot.title = element_text(size=18)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('Resident virulence') + ylab('') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + 
      theme(plot.title = element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/Fig1.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=4, top=NULL, common.legend = TRUE, legend="bottom"), width=15, height=4, units='in', dpi=600)




