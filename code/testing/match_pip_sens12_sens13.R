# Load libraries ---------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(treemapify)

# Do testing -------------------------------------------------------------------
vir_steps <- seq(2, 100, 1)
dir <- '~/virEvol/code_output/pips/'
names <- c('mod1', 'mod1_highvirselect', 'mod2', 'mod2_highvirselect', 'mod4', 'mod4_highvirselect')
#names <- c('mod1', 'mod1_highvirselect')
for (name in names) {
  print(name)
  main_pip <- read.csv(paste0(dir, name, '.csv'))
  main_pip <- main_pip[,2:ncol(main_pip)]
  sens12_pip <- read.csv(paste0(dir, name, '_sens12.csv'))
  sens12_pip <- sens12_pip[,2:ncol(sens12_pip)]
  sens13_pip <- read.csv(paste0(dir, name, '_sens13.csv'))
  sens13_pip <- sens13_pip[,2:ncol(sens13_pip)]
  if (nrow(main_pip) != ncol(main_pip) | nrow(sens12_pip) != ncol(sens12_pip) | nrow(sens13_pip) != ncol(sens13_pip)) {
    stop('Error in number of rows or number of columns.')
  }
  if (nrow(main_pip) != length(vir_steps) | nrow(sens12_pip) != length(vir_steps) | nrow(sens13_pip) != length(vir_steps)) {
    stop('Error in number of vir_steps.')
  }
  if ((nrow(main_pip) != nrow(sens12_pip)) | (ncol(main_pip) != ncol(sens12_pip))) {
    stop('Error in sens12')
  }
  if ((nrow(main_pip) != nrow(sens13_pip)) | (ncol(main_pip) != ncol(sens13_pip))) {
    stop('Error in sens12')
  }
  main_pip_sens12 <- main_pip # 100 will mean there is a difference
  main_pip_sens13 <- main_pip
  print('+180 days res:')
  for (i in 1:nrow(main_pip)) {
    for (j in 1:ncol(main_pip)) {
      if (main_pip[i,j] != sens12_pip[i,j]) {
        print(paste0(vir_steps[i], '; ', vir_steps[j]))
        main_pip_sens12[i, j] <- 100
      }
    }
  }
  print('-180 days res:')
  for (i in 1:nrow(main_pip)) {
    for (j in 1:ncol(main_pip)) {
      if (main_pip[i,j] != sens13_pip[i,j]) {
        print(paste0(vir_steps[i], '; ', vir_steps[j]))
        main_pip_sens13[i, j] <- 100
      }
    }
  }
  # Plot: main_pip with red where there is a difference from sens12
  l <- list()
  l_dex <- 1
  pip <- main_pip_sens12
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0 & pip_toPlot != 100)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  l[[l_dex]] <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray", 'red'), na.value='white') +theme(legend.position="none") +
    xlab('') + ylab('') + ggtitle(paste0(name, '_sens12')) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  l_dex <- l_dex + 1
  # Plot: main_pip with red where there is a difference from sens13
  pip <- main_pip_sens13
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0 & pip_toPlot != 100)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  l[[l_dex]] <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray", 'red'), na.value='white') +theme(legend.position="none") +
    xlab('') + ylab('') + ggtitle(paste0(name, '_sens13')) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  # Save jpeg
  ggsave(filename=paste0("~/virEvol/code_output/testing/sens12_sens13/", name, ".jpg"), marrangeGrob(grobs = l, nrow=2, ncol=1, top=NULL, common.legend = FALSE), width=4, height=8, units='in', dpi=600)
}







