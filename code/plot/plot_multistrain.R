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
hist(as.numeric(forhist))

# Model 2: SEIR with vaccination -----------------------------------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod2_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 1100, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist))

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod2_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 10100, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist))

# Model 3: SEIR with vaccination and migration ---------------------------------
finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_interyear10_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(11, 1100, 11)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist))

finalMatrix <- read.csv('~/virEvol/code_output/multistrain_res/mod3_interyear1_maxyear100.csv')
finalMatrix <- finalMatrix[,2:ncol(finalMatrix)]
forhist <- c()
for (row_dex in seq(101, 10100, 101)) {
  forhist <- c(forhist, finalMatrix[row_dex,][which(!is.na(finalMatrix[row_dex,]))])
}
hist(as.numeric(forhist))