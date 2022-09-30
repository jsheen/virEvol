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
vir_steps <- seq(2, 100, 5)

# Function to plot each of the results -----------------------------------------
plot_pip <- function(name, vir_steps, title) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/main/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip == 4), 0, pip_toPlot) # Both extinct
  pip_toPlot <- ifelse((pip == 2), 0, pip_toPlot) # Res extinct
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
    col_singular_strat <- singular_strats[floor(length(singular_strats) / 2)]
    vir_singular_strat <- vir_steps[col_singular_strat]
    print(paste0(name, ": ", vir_singular_strat))
  }
}
# Figure 3
find_singular_strat('mod1', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff', vir_steps=seq(2, 100, 1))
find_singular_strat('mod1_highvirselect', vir_steps=seq(2, 100, 1))
find_singular_strat('mod2_highvirselect', vir_steps=seq(2, 100, 1))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect', vir_steps=seq(2, 100, 5))

# Sensitivity test 1
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff', vir_steps=seq(2, 100, 1)) # Ref
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_sens1', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_sens2', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_sens3', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_sens4', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_sens5', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_sens6', vir_steps=seq(2, 100, 5))

# Sensitivity test 2
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect', vir_steps=seq(2, 100, 5)) # Ref
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens1', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens2', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens3', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens4', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens5', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens6', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect_sens7', vir_steps=seq(2, 100, 5))

# Figure 1: Conceptual diagram -------------------------------------------------
# Done in Google slides

# Figure 2: Model digratm ------------------------------------------------------
# Done in Latex

# Figure 3: Main comparison ----------------------------------------------------
l <- list()
l_dex <- 1
titles <- c("(A) Without vaccination\n  (SEIR)", "(D) Without vaccination\n  (SEIR)",  
            "(B) With vaccination\n  (SEIRV)", "(E) With vaccination\n  (SEIRV)", 
            "(C) With vaccination\n      & markets (SEIRV-m)", "(F) With vaccination\n      & markets (SEIRV-m)")
for (name in c('mod1','mod1_highvirselect',
               'mod2', 'mod2_highvirselect',
               'mod4_v33_mfm33_mmf7_c5_nodiff','mod4_v33_mfm33_mmf7_c5_nodiff_highvirselect')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  if (l_dex %in% c(6)) {
    vir_steps <- seq(2, 100, 5)
  } else {
    vir_steps <- seq(2, 100, 1)
  }
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Extinction of resident before introduction of invader
  pip_toPlot <- ifelse((pip_toPlot == 4), NA, pip_toPlot) # Extinction of both strains by the end of ten years
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip: ', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
    xlab('') + ylab('') + ggtitle(titles[l_dex]) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/main/Fig3.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=3, top=NULL, common.legend = TRUE, legend="bottom"), width=12, height=8, units='in', dpi=600)

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

# Figure S2: Vary vaccination and migration (sensitivity tests 4 and 5) --------
vir_steps <- seq(2, 100, 5)
l <- list()
l_dex <- 1
titles <- c("(A) Vax = 33%\nMigrate = 33%",  "(C) Vax = 66%\nMigrate = 33%", 
            "(B) Vax = 33%\nMigrate = 66%", "(D) Vax = 66%\nMigrate = 66%")
for (name in c('mod4_v33_mfm33_mmf7_c5_nodiff', 'mod4_v66_mfm33_mmf7_c5_nodiff', 
               'mod4_v33_mfm66_mmf7_c5_nodiff', 'mod4_v66_mfm66_mmf7_c5_nodiff')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/main/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Coexistence counted as successful invasion
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
    xlab('Resident virulence') + ylab('Invader virulence') + ggtitle(titles[l_dex]) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS2.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=2, top=NULL, common.legend = TRUE, legend="bottom"), width=8, height=8, units='in', dpi=600)

# Figure S3: sensitivity tests 1 through 3 -------------------------------------
vir_steps <- seq(2, 100, 5)
l <- list()
l_dex <- 1
titles <- c("(A) Reference",  "(C) No symptomatic migration", 
            "(B) Differential migration", "(D) No slaughtered poultry")
for (name in c('mod4_v33_mfm33_mmf7_c5_nodiff', 'mod5_v33_mfm33_mmf7_c5_nodiff', 
               'mod3_v33_mfm33_mmf7_c10_diff', 'mod3_v33_mfm33_mmf7_c10_nodiff')) {
  if (l_dex == 1 | l_dex == 2 | l_dex == 4) {
    pip <- read.csv(paste0('~/virEvol/code_output/pips/main/', name, '.csv'))
  } else {
    pip <- read.csv(paste0('~/virEvol/code_output/pips/supplementary/', name, '.csv'))
  }
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  if (l_dex == 3 | l_dex == 4) {
    vir_steps <- seq(2, 100, 1)
  } else {
    vir_steps <- seq(2, 100, 5)
  }
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Coexistence counted as successful invasion
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.', name)) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
    xlab('Resident virulence') + ylab('Invader virulence') + ggtitle(titles[l_dex]) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS3.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=4, top=NULL, common.legend = TRUE, legend="bottom"), width=18, height=4, units='in', dpi=600)

# Figure 2: turnover rate from markets and contact rate ------------------------
vary_removal <- data.frame(matrix(c(7, 47,
                                    30, 42,
                                    60, 37), byrow=T, nrow=3, ncol=2))
colnames(vary_removal) <- c('Duration', 'ESS')
vary_contact <- data.frame(matrix(c(1, 44,
                                    2, 44,
                                    5, 38,
                                    10, 37), byrow=T, nrow=4, ncol=2))
colnames(vary_contact) <- c('Ratio', 'ESS')
l <- list()
l_dex <- 1
titles <- c(paste("(A) ESS vs. duration in\nmarkets"," (ρ = 10)"), paste("(B) ESS vs. ", "ρ", "(duration in\nmarkets = 60 days)"))
for (plot.df in list(vary_removal, vary_contact)) {
  if (l_dex == 1) {
    temp_plot <- ggplot(plot.df, aes(y = ESS, x = Duration)) + geom_line(size=1) + 
      theme(legend.position="none") + geom_hline(yintercept=47, linetype='dashed', color='red', size=2) + geom_point(size=5) +
      xlab('Duration in markets (days)') + ylab(expression(paste("α"^"*"))) + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18))
  } else if (l_dex == 2) {
    temp_plot <- ggplot(plot.df, aes(y = ESS, x = Ratio)) + geom_line(size=1) +
      theme(legend.position="none") + geom_hline(yintercept=47, linetype='dashed', color='red', size=2) + geom_point(size=5) +
      xlab(expression(paste("ρ"))) + ylab(expression(paste("α"^"*"))) + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/main/Fig2.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=2, top=NULL, common.legend = TRUE, legend="bottom"), width=8, height=4, units='in', dpi=600)

# Figure S4: Faster rate of migration from market patch results and lower
#            transmission rate in markets --------------------------------------
l <- list()
l_dex <- 1
titles <- c("(A) Duration in Markets\n(days) = 7 (ρ = 10)", 
            "(B) Duration in Markets\n(days) = 30 (ρ = 10)",
            "(C) Duration in Markets\n(days) = 60 (ρ = 10)",
            "(D) ρ = 1 (duration in\nmarkets = 60 days)",
            "(E) ρ = 2 (duration in\nmarkets = 60 days)",
            "(F) ρ = 5 (duration in\nmarkets = 60 days)")
for (name in c('mod3_v33_mfm33_mmf7_c10_nodiff', 
               'mod3_v33_mfm33_mmf30_c10_nodiff',
               'mod3_v33_mfm33_mmf60_c10_nodiff',
               'mod3_v33_mfm33_mmf60_c1_nodiff', 
               'mod3_v33_mfm33_mmf60_c2_nodiff',
               'mod3_v33_mfm33_mmf60_c5_nodiff')) {
  if (l_dex == 1) {
    pip <- read.csv(paste0('~/virEvol/code_output/pips/main/', name, '.csv'))
  } else {
    pip <- read.csv(paste0('~/virEvol/code_output/pips/supplementary/', name, '.csv'))
  }
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  if (l_dex == 2) {
    vir_steps <- seq(2, 100, 5)
  } else {
    vir_steps <- seq(2, 100, 1)
  }
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
    scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
    xlab('Resident virulence') + ylab('Invader virulence') + ggtitle(titles[l_dex]) +
    theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS4.jpg"), marrangeGrob(grobs = l, ncol=3, nrow=2, top=NULL, layout_matrix=matrix(1:6, 2, 3, TRUE), common.legend = TRUE, legend="bottom"), width=12, height=8, units='in', dpi=600)

# Plot Figure 3 (virulence = 50)  ----------------------------------------------
l <- list()
vir_50_mod2 <- data.frame(matrix(c(76.32, 23.68,
                                          'Unvax, Farm', 'Vax, Farm'), ncol=2, nrow=2))
colnames(vir_50_mod2) <- c('infectious', 'patch')
vir_50_mod2$infectious <- as.numeric(vir_50_mod2$infectious)
l[[1]] <- ggplot(vir_50_mod2, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(A) No Migration From Markets\nNo Market Transmission') + theme(plot.title = element_text(size=19))
vir_50_mod3_nodiff <- data.frame(matrix(c(75.63, 23.46, 0.88, 0.03,
                                          'Unvax, Farm', 'Vax, Farm',
                                          'Unvax, Market', 'Vax, Market'), ncol=2, nrow=4))
colnames(vir_50_mod3_nodiff) <- c('infectious', 'patch')
vir_50_mod3_nodiff$infectious <- as.numeric(vir_50_mod3_nodiff$infectious)
l[[2]] <- ggplot(vir_50_mod3_nodiff, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(B) Fast Migration From Markets\nHigh Market Transmission') + theme(plot.title = element_text(size=19))
vir_50_mod3_nodiff_slow <- data.frame(matrix(c(46.21, 12.84, 40.5, 0.46,
                                          'Unvax, Farm', 'Vax, Farm',
                                          'Unvax, Market', 'Vax, Market'), ncol=2, nrow=4))
colnames(vir_50_mod3_nodiff_slow) <- c('infectious', 'patch')
vir_50_mod3_nodiff_slow$infectious <- as.numeric(vir_50_mod3_nodiff_slow$infectious)
l[[3]] <- ggplot(vir_50_mod3_nodiff_slow, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(C) Slow Migration From Markets\nHigh Market Transmission') + theme(plot.title = element_text(size=19))
vir_50_mod3_nodiff_slow_lowc <- data.frame(matrix(c(76.16, 21.65, 2.12, 0.06,
                                               'Unvax, Farm', 'Vax, Farm',
                                               'Unvax, Market', 'Vax, Market'), ncol=2, nrow=4))
colnames(vir_50_mod3_nodiff_slow_lowc) <- c('infectious', 'patch')
vir_50_mod3_nodiff_slow_lowc$infectious <- as.numeric(vir_50_mod3_nodiff_slow_lowc$infectious)
l[[4]] <- ggplot(vir_50_mod3_nodiff_slow_lowc, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(D) Slow Migration From Markets\nLow Market Transmission') + theme(plot.title = element_text(size=19))
ggsave(filename=paste0("~/virEvol/code_output/plots/main/Fig3.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=4, top=NULL, legend="bottom"), width=16, height=4.5, units='in', dpi=600)


