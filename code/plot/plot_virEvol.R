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
find_singular_strat <- function(name, vir_steps, supplementary=FALSE) {
  if (supplementary) {
    pip <- read.csv(paste0('~/virEvol/code_output/pips/supplementary/', name, '.csv'))
  } else {
    pip <- read.csv(paste0('~/virEvol/code_output/pips/main/', name, '.csv'))
  }
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
# Figure 1 singular strategies
find_singular_strat('mod1', vir_steps=vir_steps)
find_singular_strat('mod2', vir_steps=vir_steps)
find_singular_strat('mod4_v33_mfm33_mmf7_c5_nodiff', vir_steps=seq(2, 100, 5))

# Figure S2 singular strategies
find_singular_strat('mod4_v66_mfm33_mmf7_c5_nodiff', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v33_mfm66_mmf7_c5_nodiff', vir_steps=seq(2, 100, 5))
find_singular_strat('mod4_v66_mfm66_mmf7_c5_nodiff', vir_steps=seq(2, 100, 5))

# Figure S3 singular strategies
find_singular_strat('mod3_v33_mfm33_mmf7_c10_diff', vir_steps=seq(2, 100, 1), supplementary=TRUE)
find_singular_strat('mod5_v33_mfm33_mmf7_c5_nodiff', vir_steps=seq(2, 100, 5))
find_singular_strat('mod3_v33_mfm33_mmf7_c10_nodiff', vir_steps=seq(2, 100, 1))


find_singular_strat('mod3_v33_mfm33_mmf7_c10_nodiff', vir_steps=vir_steps)
find_singular_strat('mod3_v33_mfm33_mmf60_c5_nodiff', vir_steps=vir_steps, supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf60_c10_nodiff', vir_steps=vir_steps, supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf30_c10_nodiff', vir_steps=seq(2, 100, 5), supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf45_c10_nodiff', vir_steps=seq(2, 100, 5), supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf60_c10_diff', vir_steps=vir_steps, supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf7_c10_diff', vir_steps=vir_steps, supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf60_c1_nodiff', vir_steps=vir_steps, supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf60_c1_diff', vir_steps=vir_steps, supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf60_c2_nodiff', vir_steps=vir_steps, supplementary=TRUE)
find_singular_strat('mod3_v33_mfm33_mmf60_c2_diff', vir_steps=vir_steps, supplementary=TRUE)


# Figure 1: Main comparison ----------------------------------------------------
l <- list()
l_dex <- 1
titles <- c("(A) Without vaccination\n  (SEIR)", "(B) With vaccination\n  (SEIRV)", 
            "(C) With vaccination\n      & markets (SEIRV-m)")
for (name in c('mod1', 'mod2', 'mod4_v33_mfm33_mmf7_c5_nodiff')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/main/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  if (l_dex == 3) {
    vir_steps <- seq(2, 100, 5)
  } else {
    vir_steps <- seq(2, 100, 1)
  }
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), NA, pip_toPlot) # Should be grayed
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.')) }
  melted <- melt(pip_toPlot)
  colnames(melted) <- c('invader_virulence', 'resident_virulence', 'value')
  if (l_dex == 1) {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('Resident virulence') + ylab('Invader virulence') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else {
    temp_plot <- ggplot(melted, aes(y = invader_virulence, x = resident_virulence, fill = factor(value))) + geom_tile() +
      scale_fill_manual(values = c("black", "gray"), na.value='white') +theme(legend.position="none") +
      xlab('Resident virulence') + ylab('') + ggtitle(titles[l_dex]) +
      theme(text = element_text(size = 24)) + theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=18)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  l[[l_dex]] <- temp_plot
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/main/Fig1.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=3, top=NULL, common.legend = TRUE, legend="bottom"), width=12, height=4, units='in', dpi=600)

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



# Plot Figure 2 (virulence = 50)  ----------------------------------------------
l <- list()
vir_50_mod2 <- data.frame(matrix(c(76.32, 23.68,
                                          'Unvax, Farm', 'Vax, Farm'), ncol=2, nrow=2))
colnames(vir_50_mod2) <- c('infectious', 'patch')
vir_50_mod2$infectious <- as.numeric(vir_50_mod2$infectious)
l[[1]] <- ggplot(vir_50_mod2, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(A) High Contact\nNo Migration Market to Farm')
vir_50_mod3_nodiff <- data.frame(matrix(c(40.26, 10.24, 45.03, 4.46,
                                          'Unvax, Farm', 'Vax, Farm',
                                          'Unvax, Market', 'Vax, Market'), ncol=2, nrow=4))
colnames(vir_50_mod3_nodiff) <- c('infectious', 'patch')
vir_50_mod3_nodiff$infectious <- as.numeric(vir_50_mod3_nodiff$infectious)
l[[2]] <- ggplot(vir_50_mod3_nodiff, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(B) High Contact,\nLow Migration Rate Market to Farm')
vir_50_mod3_nodiff_fast <- data.frame(matrix(c(75.79, 22.99, 0.92, 0.31,
                                          'Unvax, Farm', 'Vax, Farm',
                                          'Unvax, Market', 'Vax, Market'), ncol=2, nrow=4))
colnames(vir_50_mod3_nodiff_fast) <- c('infectious', 'patch')
vir_50_mod3_nodiff_fast$infectious <- as.numeric(vir_50_mod3_nodiff_fast$infectious)
l[[3]] <- ggplot(vir_50_mod3_nodiff_fast, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(C) High Contact,\nHigh Migration Rate Market to Farm')
vir_50_mod3_nodiff_lowc <- data.frame(matrix(c(76.9, 20.02, 2.47, 0.61,
                                               'Unvax, Farm', 'Vax, Farm',
                                               'Unvax, Market', 'Vax, Market'), ncol=2, nrow=4))
colnames(vir_50_mod3_nodiff_lowc) <- c('infectious', 'patch')
vir_50_mod3_nodiff_lowc$infectious <- as.numeric(vir_50_mod3_nodiff_lowc$infectious)
l[[4]] <- ggplot(vir_50_mod3_nodiff_lowc, aes(area = infectious, fill = infectious, label = patch)) + geom_treemap() +
  geom_treemap_text(fontface = "italic", colour = "white", place = "centre",
                    grow = TRUE) + scale_fill_continuous(name = "% Infectious") + theme(legend.position="none") +
  ggtitle('(D) Low Contact,\nLow Migration Rate Market to Farm')
ggsave(filename=paste0("~/virEvol/code_output/plots/main/Fig2.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=4, top=NULL, legend="bottom"), width=12, height=4, units='in', dpi=600)



# Figure S1: Lower contact rate in markets -------------------------------------
l <- list()
l_dex <- 1
titles <- c("(A) Contact rate ratio in\nMarket vs. Farms = 1", 
            "(B) Contact rate ratio in\nMarket vs. Farms = 2",
            "(C) Contact rate ratio in\nMarket vs. Farms = 10")
for (name in c('mod3_v33_mfm33_mmf60_c1_nodiff', 
               'mod3_v33_mfm33_mmf60_c2_nodiff',
               'mod3_v33_mfm33_mmf60_c10_nodiff')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/supplementary/', name, '.csv'))
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
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS1.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=3, top=NULL, common.legend = TRUE, legend="bottom"), width=12, height=4, units='in', dpi=600)

# Figure S2: Faster rate of migration from market patch results ----------------
l <- list()
l_dex <- 1
titles <- c("(A) Migration from market\n to farm = 7", 
            "(B) Migration from market\n to farm = 30")
for (name in c('mod3_v33_mfm33_mmf7_c10_nodiff', 
               'mod3_v33_mfm33_mmf30_c10_nodiff')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/supplementary/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  if (l_dex == 1) {
    vir_steps = vir_steps = vir_steps = seq(2, 100, 1)
  } else {
    vir_steps = vir_steps = vir_steps = seq(2, 100, 5)
  }
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
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS2.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=2, top=NULL, common.legend = TRUE, legend="bottom"), width=8, height=4, units='in', dpi=600)

# Figure S3 ---------------------------------------------------------------------
l <- list()
l_dex <- 1
titles <- c("(A) Vax = 33%;\nMigrate = 33%", "(B) Vax = 33%;\nMigrate = 66%", 
            "(C) Vax = 66%;\nMigrate = 33%", "(D) Vax = 66%;\nMigrate = 66%")
for (name in c('mod3_v33_mfm33_mmf7_c10_nodiff', 'mod3_v33_mfm66_mmf7_c5_nodiff',
               'mod3_v66_mfm33_mmf7_c5_nodiff', 'mod3_v66_mfm66_mmf7_c5_nodiff')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/supplementary/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  if (l_dex == 1) {
    vir_steps <- seq(2, 100, 1)
  } else {
    vir_steps <- seq(2, 100, 5)
  }
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
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS3.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=2, top=NULL, common.legend = TRUE, legend="bottom"), width=8, height=8, units='in', dpi=600)

# Figure S5 --------------------------------------------------------------------
vir_steps <- seq(2, 100, 5)
l <- list()
l_dex <- 1
titles <- c("(A) Vax = 33%;\nMigrate = 33%", "(B) Vax = 33%;\nMigrate = 66%", 
            "(C) Vax = 66%;\nMigrate = 33%", "(D) Vax = 66%;\nMigrate = 66%")
for (name in c('mod5_v33_mfm33_mmf7_c5_nodiff', 'mod5_v33_mfm66_mmf7_c5_nodiff',
               'mod5_v66_mfm33_mmf7_c5_nodiff', 'mod5_v66_mfm66_mmf7_c5_nodiff')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/main/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), 0, pip_toPlot) # TO REPLACE
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.', name)) }
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
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS5.jpg"), marrangeGrob(grobs = l, nrow=2, ncol=2, top=NULL, common.legend = TRUE, legend="bottom"), width=8, height=8, units='in', dpi=600)

# Figure S6 --------------------------------------------------------------------
vir_steps <- seq(2, 100, 3)
l <- list()
l_dex <- 1
titles <- c("Lowered Transmission")
for (name in c('mod1_c0.5')) {
  pip <- read.csv(paste0('~/virEvol/code_output/pips/supplementary/', name, '.csv'))
  pip <- pip[,c(2:ncol(pip))]
  pip <- data.matrix(pip)
  colnames(pip) <- vir_steps
  rownames(pip) <- rev(vir_steps)
  pip_toPlot <- ifelse((pip == 3), 1, pip) # Coexistence counted as successful invasion
  pip_toPlot <- ifelse((pip_toPlot == 2), 0, pip_toPlot) # TO REPLACE
  if (length(which(pip_toPlot != 1 & pip_toPlot != 0)) > 0) { print(paste0('Error in pip.', name)) }
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
ggsave(filename=paste0("~/virEvol/code_output/plots/supplementary/FigS6.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=1, top=NULL, common.legend = TRUE, legend="bottom"), width=8, height=8, units='in', dpi=600)

# Should put a supplementary plot with differential migration of vaccinated individuals


