# ------------------------------------------------------------------------------
# Plotting functions
# ------------------------------------------------------------------------------
# Load libraries ---------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(gridExtra)

# Model 1: SEIR ---------------------------------------------------------------- 
pip <- read.csv(paste0('~/virEvol/code_output/pips/mod1.csv'))
pip <- pip[,c(2:ncol(pip))]
pip <- data.matrix(pip)
colnames(pip) <- vir_steps
rownames(pip) <- rev(vir_steps)
pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
melted <- melt(pip_toPlot)
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value='grey') +geom_abline(intercept = 0, slope = 1, col='red') +theme(legend.position="none")
l <- list()
l[[1]] <- temp_plot
ggsave(filename="~/virEvol/code_output/plots/virEvol_mod1.png", marrangeGrob(grobs = l, nrow=1, ncol=1), width=3, height=3, units='in', dpi=600)

# Model 2: SEIR with vaccination -----------------------------------------------
pip <- read.csv(paste0('~/virEvol/code_output/pips/mod2.csv'))
pip <- pip[,c(2:ncol(pip))]
pip <- data.matrix(pip)
colnames(pip) <- vir_steps
rownames(pip) <- rev(vir_steps)
pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
melted <- melt(pip_toPlot)
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value='grey') +geom_abline(intercept = 0, slope = 1, col='red') +theme(legend.position="none")
l <- list()
l[[1]] <- temp_plot
ggsave(filename="~/virEvol/code_output/plots/virEvol_mod2.png", marrangeGrob(grobs = l, nrow=1, ncol=1), width=3, height=3, units='in', dpi=600)

# Model 3: SEIR with vaccination and migration, differential migration ---------
pip <- read.csv(paste0('~/virEvol/code_output/pips/mod3_diff.csv'))
pip <- pip[,c(2:ncol(pip))]
pip <- data.matrix(pip)
colnames(pip) <- vir_steps
rownames(pip) <- rev(vir_steps)
pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
melted <- melt(pip_toPlot)
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value='grey') +geom_abline(intercept = 0, slope = 1, col='red') +theme(legend.position="none")
l <- list()
l[[1]] <- temp_plot
ggsave(filename="~/virEvol/code_output/plots/virEvol_mod3_diff.png", marrangeGrob(grobs = l, nrow=1, ncol=1), width=3, height=3, units='in', dpi=600)

# Model 3: SEIR with vaccination and migration, non-differential migration -----
pip <- read.csv(paste0('~/virEvol/code_output/pips/mod3_nodiff.csv'))
pip <- pip[,c(2:ncol(pip))]
pip <- data.matrix(pip)
colnames(pip) <- vir_steps
rownames(pip) <- rev(vir_steps)
pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
melted <- melt(pip_toPlot)
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value='grey') +geom_abline(intercept = 0, slope = 1, col='red') +theme(legend.position="none")
l <- list()
l[[1]] <- temp_plot
ggsave(filename="~/virEvol/code_output/plots/virEvol_mod3_nodiff.png", marrangeGrob(grobs = l, nrow=1, ncol=1), width=3, height=3, units='in', dpi=600)

# Model 3 (mpatch): SEIR with vaccination and migration, differential migration ---------
pip <- read.csv(paste0('~/virEvol/code_output/pips/mod3_diff_mpatch.csv'))
pip <- pip[,c(2:ncol(pip))]
pip <- data.matrix(pip)
colnames(pip) <- vir_steps
rownames(pip) <- rev(vir_steps)
pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
melted <- melt(pip_toPlot)
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value='grey') +geom_abline(intercept = 0, slope = 1, col='red') +theme(legend.position="none")
l <- list()
l[[1]] <- temp_plot
ggsave(filename="~/virEvol/code_output/plots/virEvol_mod3_diff_mpatch.png", marrangeGrob(grobs = l, nrow=1, ncol=1), width=3, height=3, units='in', dpi=600)

# Model 3 (mpatch): SEIR with vaccination and migration, non-differential migration -----
pip <- read.csv(paste0('~/virEvol/code_output/pips/mod3_nodiff_mpatch.csv'))
pip <- pip[,c(2:ncol(pip))]
pip <- data.matrix(pip)
colnames(pip) <- vir_steps
rownames(pip) <- rev(vir_steps)
pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
melted <- melt(pip_toPlot)
colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
  scale_fill_manual(values = c("black", "white"), na.value='grey') +geom_abline(intercept = 0, slope = 1, col='red') +theme(legend.position="none")
l <- list()
l[[1]] <- temp_plot
ggsave(filename="~/virEvol/code_output/plots/virEvol_mod3_nodiff_mpatch.png", marrangeGrob(grobs = l, nrow=1, ncol=1), width=3, height=3, units='in', dpi=600)

