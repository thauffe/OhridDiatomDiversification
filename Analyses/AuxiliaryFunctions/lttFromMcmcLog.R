lttFromMcmcLog <- function(x, ColTs, ColTe){
  TimeSeq <- seq(0, 13.7, 0.05)
  LttDf <- data.frame(Age = c(unlist(x[ColTs]), x[ColTe][x[ColTe] != 0], 0), # Add 0 for the present!
                      S = c(rep(1, length(x[ColTs])), 
                            rep(-1, length(x[ColTe][x[ColTe] != 0])), 0))
  LttDf <- LttDf[order(LttDf$Age, decreasing = TRUE), ]
  D <- approx(x = LttDf$Age, y = cumsum(LttDf$S), xout = TimeSeq, method = "constant")$y
  return(D)
}