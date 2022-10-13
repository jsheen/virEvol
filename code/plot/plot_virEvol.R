# ------------------------------------------------------------------------------
# Plotting functions
# ------------------------------------------------------------------------------
# Load libraries ---------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(treemapify)

# Load global variables --------------------------------------------------------
vir_steps <- seq(2, 100, 1)

# Function to plot each of the results -----------------------------------------
plot_pip <- function(name, vir_steps, title) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 4), 0, pip_toPlot) # Both extinct
  pip_toPlot <- ifelse((pip_toPlot == 2), 0, pip_toPlot) # Res extinct
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

# Function to find singular strategy -------------------------------------------
find_singular_strat <- function(name, vir_steps) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before invader introduced
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both resident and invader
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
  singular_strats <- c()
  for (col_dex in 1:ncol(pip)) {
    if (length(which(pip[,col_dex] == 0)) == nrow(pip)) {
      singular_strats <- c(singular_strats, col_dex)
    }
  }
  if (length(singular_strats) == 1) {
    print(paste0(name, ": ", vir_steps[singular_strats]))
  } else {
    vir_singular_strat <- mean(vir_steps[singular_strats])
    #print(paste0(name, ": ", vir_singular_strat))
    return(vir_singular_strat)
  }
}
# Tradeoff curve 1 sensitivity tests
find_singular_strat('mod4', vir_steps=seq(2, 100, 1)) # Ref
find_singular_strat('mod4_sens1', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens2', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens3', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens5', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens6', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens7', vir_steps=seq(4, 100, 1))
find_singular_strat('mod4_sens8', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens9', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens10', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens11', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod3_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens13', vir_steps=seq(2, 100, 1))

# Tradeoff curve 2 sensitivity tests
find_singular_strat('mod4_highvirselect', vir_steps=seq(2, 100, 1)) # Ref
find_singular_strat('mod4_highvirselect_sens1', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens2', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens3', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens5', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens6', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens7', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens8', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens9', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens10', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens11', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_highvirselect_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_highvirselect_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_highvirselect_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_highvirselect_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens13', vir_steps=seq(2, 100, 1))

# Figure 3
find_singular_strat('mod1', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_highvirselect', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_highvirselect', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect', vir_steps=seq(2, 100, 1))

# Figure S2
find_singular_strat('mod4_sens7', vir_steps=seq(4, 100, 1))
find_singular_strat('mod4_sens2', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens8', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens7', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens2', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens8', vir_steps=seq(2, 100, 1))

# Figure S3
find_singular_strat('mod4_sens3', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens3', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens4', vir_steps=seq(2, 100, 1))

# Figure S4
find_singular_strat('mod4_sens9', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens5', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens9', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens5', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect', vir_steps=seq(2, 100, 1))

# Figure S5
find_singular_strat('mod4_sens10', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens1', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens10', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens1', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect', vir_steps=seq(2, 100, 1))

# Figure S6
find_singular_strat('mod4_sens6', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens11', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens6', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens11', vir_steps=seq(2, 100, 1))

# Figure S7
find_singular_strat('mod1_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_highvirselect_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_highvirselect_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_highvirselect_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_highvirselect_sens13', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens12', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_highvirselect_sens13', vir_steps=seq(2, 100, 1))

# Figure 1: Conceptual diagram -------------------------------------------------
# Done in Google slides

# Figure 2: Model digratm ------------------------------------------------------
# Done in Latex

# Figure 3: Main comparison ----------------------------------------------------
l <- list()
l_dex <- 1
titles <- c("(A) Without vaccination\n  (SEIR)", "",  
            "(B) With vaccination\n  (SEIRV)", "", 
            "(C) With vaccination\n      & markets (SEIRV-m)", "")
for (name in c('mod1','mod1_highvirselect',
               'mod2', 'mod2_highvirselect',
               'mod4','mod4_highvirselect')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  l[[l_dex]] <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
    xlab('') + ylab('') + ggtitle(titles[l_dex]) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps=seq(2, 100, 1)), y=find_singular_strat(name, vir_steps=seq(2, 100, 1))), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/main/Fig3.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=3, top=NULL, common.legend = FALSE), width=12, height=8, units='in', dpi=600)

# Figure S1: Tradeoff curves ---------------------------------------------------
c1 = 1
c2 = 0.4
virulences <- seq(0, 100, 0.01)
transmissions <- (c1 * virulences)^c2
mortalities <- virulences / 100
plot(mortalities, transmissions, type='l', main='Transmission-Mortality Tradeoff Curve',
     xlab=expression(mu), ylab=expression(beta), lwd=5, col='red', cex.lab=1.2)
c1 = 0.1
c2 = 0.3
virulences <- seq(0, 100, 0.01)
transmissions <- (c1 * virulences)^c2
mortalities <- virulences / 100
lines(mortalities, transmissions, lwd=5, col='blue')
legend("topleft", 
       legend = c("Tradeoff curve 1", "Tradeoff curve 2"), 
       col = c('blue', 'red'),
       bty = "n", 
       pch = c(16, 16), 
       pt.cex = 1.2,
       cex = 1, 
       text.col = "black", 
       horiz = F)

# Figure S2: Migration rate will increase virulence (sensitivity tests 7, 2, main ref, and 8)
l <- list()
l_dex <- 1
titles <- c('0% / 4 mo.', "",  
            '5% / 4 mo.', "", 
            '10% / 4 mo.', "",
            '20% / 4 mo.', "")
for (name in c('mod4_sens7','mod4_highvirselect_sens7',
               'mod4_sens2', 'mod4_highvirselect_sens2',
               'mod4','mod4_highvirselect',
               'mod4_sens8','mod4_highvirselect_sens8')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  if (l_dex %in% c(1)) {
    vir_steps <- seq(4, 100, 1)
  }
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  if (l_dex %in% c(1, 3, 5, 7)) {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle(bquote(δ[f]==.(titles[l_dex]))) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle('') +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS2.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=4, top=NULL, common.legend = FALSE), width=16, height=8, units='in', dpi=600)

# Figure S3: Rho will not affect results (sensitivity tests 3, main ref, and 4)
l <- list()
l_dex <- 1
titles <- c('ρ = 1', "",  
            'ρ = 5', "", 
            'ρ = 10', "")
for (name in c('mod4_sens3','mod4_highvirselect_sens3',
               'mod4','mod4_highvirselect',
               'mod4_sens4','mod4_highvirselect_sens4')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  if (l_dex %in% c(1, 3, 5)) {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18))  + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle('') +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS3.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=3, top=NULL, common.legend = FALSE), width=12, height=8, units='in', dpi=600)

# Figure S4: Increasing probability of slaughter will increase ESS  (sensitivity tests 9, 5, and main)
l <- list()
l_dex <- 1
titles <- c('0', "",  
            '0.5', "", 
            '0.8', "")
for (name in c('mod4_sens9','mod4_highvirselect_sens9',
               'mod4_sens5','mod4_highvirselect_sens5',
               'mod4','mod4_highvirselect')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  if (l_dex %in% c(1, 3 ,5)) {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle(bquote(p[s]==.(titles[l_dex]))) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle('') +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS4.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=3, top=NULL, common.legend = FALSE), width=12, height=8, units='in', dpi=600)

# Figure S5: Increasing the rate of vaccination will increase the ESS (sensitivity tests 10, 1, and main)
l <- list()
l_dex <- 1
titles <- c('v = 0% / 4 mo.', "",  
            'v = 13% / 4 mo.', "", 
            'v = 26% / 4 mo.', "")
for (name in c('mod4_sens10','mod4_highvirselect_sens10',
               'mod4_sens1','mod4_highvirselect_sens1',
               'mod4','mod4_highvirselect')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  if (l_dex %in% c(1, 5)) {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS5.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=3, top=NULL, common.legend = FALSE), width=12, height=8, units='in', dpi=600)

# Figure S6: Miscellaneous tests (sensitivity tests 6 and 11)
l <- list()
l_dex <- 1
titles <- c('No migration of infectious', "",  
            '1 / 60 days', "")
for (name in c('mod4_sens6','mod4_highvirselect_sens6',
               'mod4_sens11','mod4_highvirselect_sens11')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  vir_steps <- seq(2, 100, 1)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  if (l_dex %in% c(1)) {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else if (l_dex %in% c(3)){
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle(bquote(δ[m]==.(titles[l_dex]))) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('') + ylab('') + ggtitle("") +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + geom_point(aes_q(x=find_singular_strat(name, vir_steps), y=find_singular_strat(name, vir_steps)), colour="red", size=3) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS6.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=2, top=NULL, common.legend = FALSE), width=8, height=8, units='in', dpi=600)

# Figure 7: robustness to introduction of invader at other time steps ----------
vir_steps <- seq(2, 100, 1)
dir <- '~/virEvol/code_output/pips/'
names <- c('mod1', 'mod1_highvirselect', 'mod2', 'mod2_highvirselect', 'mod4', 'mod4_highvirselect')
titles <- c("(A) Without vaccination\n  (SEIR)", "",  
            "(B) With vaccination\n  (SEIRV)", "", 
            "(C) With vaccination\n      & markets (SEIRV-m)", "")
paper_ls <- list()
paper_ls_dex <- 1
for (name in names) {
  print(name)
  main_pip <- read.csv(paste0(dir, name, '.csv'))
  main_pip <- main_pip[,2:ncol(main_pip)]
  sens12_pip <- read.csv(paste0(dir, name, '_sens12.csv'))
  sens12_pip <- sens12_pip[,2:ncol(sens12_pip)]
  sens13_pip <- read.csv(paste0(dir, name, '_sens13.csv'))
  sens13_pip <- sens13_pip[,2:ncol(sens13_pip)]
  sens14_pip <- read.csv(paste0(dir, name, '_sens14.csv'))
  sens14_pip <- sens14_pip[,2:ncol(sens14_pip)]
  sens15_pip <- read.csv(paste0(dir, name, '_sens15.csv'))
  sens15_pip <- sens15_pip[,2:ncol(sens15_pip)]
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
  main_pip_sens14 <- main_pip
  main_pip_sens15 <- main_pip
  paper_pip <- main_pip
  print('+180 days res:')
  for (i in 1:nrow(main_pip)) {
    for (j in 1:ncol(main_pip)) {
      if (main_pip[i,j] != sens12_pip[i,j]) {
        print(paste0(vir_steps[i], '; ', vir_steps[j]))
        main_pip_sens12[i, j] <- 100
        paper_pip[i, j] <- 100
      }
    }
  }
  print('-180 days res:')
  for (i in 1:nrow(main_pip)) {
    for (j in 1:ncol(main_pip)) {
      if (main_pip[i,j] != sens13_pip[i,j]) {
        print(paste0(vir_steps[i], '; ', vir_steps[j]))
        main_pip_sens13[i, j] <- 100
        paper_pip[i, j] <- 100
      }
    }
  }
  print('+365 days res:')
  for (i in 1:nrow(main_pip)) {
    for (j in 1:ncol(main_pip)) {
      if (main_pip[i,j] != sens14_pip[i,j]) {
        print(paste0(vir_steps[i], '; ', vir_steps[j]))
        main_pip_sens14[i, j] <- 100
        paper_pip[i, j] <- 100
      }
    }
  }
  print('-365 days res:')
  for (i in 1:nrow(main_pip)) {
    for (j in 1:ncol(main_pip)) {
      if (main_pip[i,j] != sens15_pip[i,j]) {
        print(paste0(vir_steps[i], '; ', vir_steps[j]))
        main_pip_sens15[i, j] <- 100
        paper_pip[i, j] <- 100
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
  # Plot: pip for paper (sens12 through sens15 differences)
  pip <- paper_pip
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
  paper_ls[[paper_ls_dex]] <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray", 'red'), na.value='white') +theme(legend.position="none") +
    xlab('') + ylab('') + ggtitle(titles[paper_ls_dex]) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  paper_ls_dex <- paper_ls_dex + 1
  # Save jpeg
  ggsave(filename=paste0("~/virEvol/code_output/testing/sens12_sens15/", name, ".jpg"), marrangeGrob(grobs = l, nrow=1, ncol=2, top=NULL, common.legend = FALSE), width=4, height=8, units='in', dpi=600)
}
# Save jpeg for paper
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS7.jpg"), marrangeGrob(grobs = paper_ls, nrow=2, ncol=3, top=NULL, common.legend = FALSE), width=12, height=8, units='in', dpi=600)

