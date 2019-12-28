plotCovariate <- function(Time, Data, LR04, ylab = "", 
                          col = "black", ylim = NULL, Raw = NULL){
  par(las = 1, mar = c(5.1, 4.5, 0.1, 0.8))
  plot(Time, Data, 
       xlim = c(1.4, 0), ylim = ylim, xaxs = "i",
       ylab = ylab, 
       xlab = "Age (Myr ago)", type = "n")
  # MIS
  for(i in 1:nrow(LR04)){
    if( LR04[i,2] < par("usr")[1] & (LR04[i,1] %% 2 == 1) ){
      rect(xleft = LR04[i, 2], ybottom = par("usr")[3], 
           xright = LR04[i-1, 2], ytop = par("usr")[4], 
           border = NA, col = adjustcolor("grey", alpha.f = 0.1))
    }
  }
  if(!is.null(Raw)){
    lines(Raw[, 1] / 1000, Raw[, 2], col = "grey")
  }
  lines(Time / 1000, Data, col = col)
  box()
}
