morts <- seq(0, 100, 0.1)
trans <- (0.3 * morts) ^ 0.4
jpeg('~/virEvol/code_output/plots/main/conceptual.jpeg', width=5, height=5, unit='in', res=600)
plot(morts, trans, type='l', xlab='increasing mortality rate', 
     ylab='increasing transmission rate', xaxt='n', yaxt='n', line=1, col='red', lwd=2)
dev.off()
