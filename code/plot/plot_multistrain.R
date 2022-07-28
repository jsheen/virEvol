# Load libraries ---------------------------------------------------------------
library(ggplot2)

# Plotting function for low vs. high virulence strains -------------------------
plot_prevs <- function(forhist, threshold) {
  forhist <- as.numeric(forhist)
  low_prev <- c()
  high_prev <- c()
  for (year in 1:101) {
    forhist <- c()
    for (row_dex in seq(year, 10100, 101)) {
      forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
    }
    forhist <- as.numeric(forhist)
    low_prev <- c(low_prev, length(which(forhist < threshold)) / length(forhist))
    high_prev <- c(high_prev, length(which(forhist >= threshold)) / length(forhist))
  }
  plot(1:101, low_prev, ylim=c(0, 1), type='l', xlab='year', ylab='prevalence', main='low vs. high virulent strain prev.')
  lines(1:101, high_prev, ylim=c(0, 1), type='l', col='red')
  legend('topright', legend=c("low vir.", "high vir."),
         col=c("black", "red"), lty=c(1, 1), cex=0.8)
}

# Model 1: SEIR ----------------------------------------------------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod1_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 11000, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod1_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 101000, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100, xlim=c(0,100))

# Model 2: SEIR with vaccination -----------------------------------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod2_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 1100, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod2_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 101000, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100, xlim=c(0, 100))

# Model 3: SEIR with vaccination and differential migration --------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 1100, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 101000, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

# Model 3: SEIR with vaccination and differential migration (large market patch)
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_mpatch_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 101000, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

# Figure 2 ---------------------------------------------------------------------
year_horizon <- 5
titles <- c('(A) Without vaccination\n', '(B) With vaccination\n', '(C) With vaccination &\ndifferential migration')
l <- list()
l_dex <- 1
for (name in c('mod1_interyear1_maxyear100', 'mod2_interyear1_maxyear100',
               'mod3_mpatch_interyear1_maxyear100')) {
  finalMatrix <- read.csv(paste0('~/virEvol/code_output/multistrain_res/', name, '.csv'))
  finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
  forhist <- c()
  for (row_dex in seq(101, 101000, 101)) {
    new_strains <- c()
    base_row <- finalMatrix[(row_dex - (year_horizon + 1)),]
    base_row <- base_row[!is.na(base_row)]
    for (i in year_horizon:1) {
      new_row <- finalMatrix[(row_dex - i),]
      new_row <- new_row[!is.na(new_row)]
      new_strain <- setdiff(new_row, base_row)
      if (length(new_strain) != 1) {
        print(row_dex)
        print(new_strain)
        stop('Error in loop.')
      }
      new_strains <- c(new_strains, new_strain)
      base_row <- new_row
    }
    last_row <- finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))]
    forhist <- c(forhist, setdiff(last_row, new_strains)[1,1])
  }
  forhist.df <- data.frame(forhist)
  if (l_dex == 1) {
    l[[l_dex]] <- ggplot(forhist.df, aes(x = forhist)) + 
      geom_histogram(aes(y = ..density..),
                     colour = 'white', fill = "cornflowerblue", bins=30) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      xlab('Final virulences') + ylab('Density') + scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
      ggtitle(titles[l_dex]) + theme(text = element_text(size = 24)) + ylim(0, 0.1) +
      theme(plot.title = element_text(size=18))
  } else {
    l[[l_dex]] <- ggplot(forhist.df, aes(x = forhist)) + 
      geom_histogram(aes(y = ..density..),
                     colour = 'white', fill = "cornflowerblue", bins=30) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      xlab('Final virulences') + ylab('') + scale_x_continuous(limits = c(0, 100), oob = scales::oob_keep) +
      ggtitle(titles[l_dex]) + theme(text = element_text(size = 24)) + ylim(0, 0.1) +
      theme(plot.title = element_text(size=18))
  }
  print(paste0(titles[l_dex], ": ", mean(forhist.df$forhist), " [", var(forhist.df$forhist), "]"))
  l_dex <- l_dex + 1
}
ggsave(filename=paste0("~/virEvol/code_output/plots/Fig2.jpg"), marrangeGrob(grobs = l, nrow=1, ncol=3, top=NULL), width=13, height=4, units='in', dpi=600)

