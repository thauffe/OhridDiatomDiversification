extractAgesFromMcmcLog <- function(BdsMcmcLog, N = 10, Out){
  ColumnTs <- grep("TS", colnames(BdsMcmcLog))
  ColumnTe <- grep("TE", colnames(BdsMcmcLog))
  # Repeat for N samples
  for(n in 1:N){
    Row <- sample(1:nrow(BdsMcmcLog), 1)
    Ts <- unlist(BdsMcmcLog[Row, ColumnTs])
    Te <- unlist(BdsMcmcLog[Row, ColumnTe])  
    # If many species have exactly the same TS, PyRate does not work
    W <- which(Ts > 13.65999999)
    Ts[W] <- sample(seq(13.66, 15, length.out = 2000), length(W))
    SeEst = data.frame(clade = 0, species = 1:length(Ts), ts = Ts, te = Te)
    write.table(SeEst, paste0(Out, n, ".txt"), 
                sep = "\t", row.names = FALSE)
  }
}