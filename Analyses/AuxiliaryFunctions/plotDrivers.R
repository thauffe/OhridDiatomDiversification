plotDrivers <- function(Bf, Klim, GammaLab = expression(gamma [lambda]), Mtext = "", SubPlot = ""){
  par(las = 1, mar = c(4.5, 11, 0.2, 0.1), cex.axis = 0.7, cex.lab = 0.7)
  layout(matrix(1:2, nrow = 1, ncol = 2), widths = c(0.64, 0.36))
  N <- nrow(Bf)
  plot(1, 1, type = "n", 
       xlim = c(Klim, 0), ylim = c(1, N),
       xlab = expression(paste(Delta, "K")), ylab = "", yaxt = "n")
  axis(side = 2, at = 1:N, labels = rownames(Bf)) 
  mtext(Mtext, side = 2, line = 8, las = 0, cex = 0.7)
  mtext(SubPlot, side = 2, line = 8, las = 1, at = N + 0.5)
  for(i in 1:N){
    BfTmp <- Bf[i, -c((ncol(Bf)-3):ncol(Bf))]
    points(BfTmp, rep(i, length(BfTmp)), pch = 19, cex = 0.6,
           col = adjustcolor(Bf[i, 8], 1/length(BfTmp)))
  }
  par(las = 1, mar = c(4.5, 0.1, 0.2, 0.1))
  plot(1, 1, type = "n", 
       xlim = c(-1.5, 1.5), ylim = c(1, N),
       xlab = GammaLab, ylab = "", yaxt = "n")
  abline(v = 0, lty = 2, col = "grey")
  for(i in 1:N){
    lines(Bf[i, 6:7], rep(i, 2), col = Bf[i, 8])
  }
}
