# Model 1: SEIR ----------------------------------------------------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod1_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 1100, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist))

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod1_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 10100, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

# Model 2: SEIR with vaccination -----------------------------------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod2_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 1100, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', breaks=100)

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod2_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 10100, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

# Model 3: SEIR with vaccination and migration ---------------------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 1100, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), breaks=100)

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 10100, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

# Model 3: SEIR with vaccination and migration (larger market patch) -----------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_mpatch_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 10100, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist), xlab='Distribution of strain virulences', main='', breaks=100)

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


