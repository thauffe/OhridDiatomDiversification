getRatesHPD <- function(SpThroughTime,
                        ExThroughTime,
                        NetDivThroughTime,
                        Prob){
  NetDivThroughTimeHPD <- ExThroughTimeHPD <- SpThroughTimeHPD <- vector(mode = "list",
                                                                         length(Prob)) 
  names(SpThroughTimeHPD) <- Prob
  names(ExThroughTimeHPD) <- Prob
  names(NetDivThroughTimeHPD) <- Prob
  for (i in 1:length(Prob)) {
    SpTmpHPD <- ExTmpHPD <- NetDivTmpHPD <- matrix(NA_real_,
                                                   nrow = 2,
                                                   ncol = ncol(SpThroughTime))
    for (y in 1:ncol(SpThroughTime)) {
      SpTmpHPD[, y] <-  HPDinterval(coda:::mcmc(SpThroughTime[, y]), prob = Prob[i])
      ExTmpHPD[, y] <-  HPDinterval(coda:::mcmc(ExThroughTime[, y]), prob = Prob[i])
      NetDivTmpHPD[, y] <-  HPDinterval(coda:::mcmc(NetDivThroughTime[, y]), prob = Prob[i])
    }
    SpThroughTimeHPD[[i]] <- SpTmpHPD
    ExThroughTimeHPD[[i]] <- ExTmpHPD
    NetDivThroughTimeHPD[[i]] <- NetDivTmpHPD
  }
  Res <- list()
  Res[[1]] <- SpThroughTimeHPD
  Res[[2]] <- ExThroughTimeHPD
  Res[[3]] <- NetDivThroughTimeHPD
  names(Res) <- c("SpHPDs", "ExHPDs", "NetDivHPDs")
  return(Res)
}

divDynamic <- function(Path,
                       Prob = seq(0.95, 0.05, by = -0.05),
                       Subphylum = NULL){
  SpRj <- list.files(Path, pattern = "sp_rates")
  ExRj <- list.files(Path, pattern = "ex_rates")
  McmcLogs <- list.files(Path, pattern = "mcmc.log")
  if (!is.null(Subphylum)) {
    SpRj <- SpRj[grepl(Subphylum, SpRj)]
    ExRj <- ExRj[grepl(Subphylum, ExRj)]
    McmcLogs <- ExRj[grepl(Subphylum, McmcLogs)]
  }
  EdgeShift <- 13.65
  TimeSeq <- seq(EdgeShift, 0, by = -0.05) 
  SpTmp <- read.table(paste0(Path, SpRj[1]), 
                      sep = "\t", fill = NA, col.names = paste0("V", 1:20))
  Nrow <- nrow(SpTmp)
  ExThroughTime <- SpThroughTime <- matrix(NA_real_, 
                                           ncol = length(TimeSeq), 
                                           nrow = length(SpRj)*Nrow)
  ExShifts <- SpShifts <- matrix(NA_real_, ncol = 8, nrow = length(SpRj) * Nrow)
  Counter <- 1
  
  for(i in 1:length(SpRj)){
    # Give a high number of col.names because each line may differ in length
    SpTmp <- read.table(paste0(Path, SpRj[i]), 
                        sep = "\t", fill = NA, col.names = paste0("V", 1:20))
    ExTmp <- read.table(paste0(Path, ExRj[i]), 
                        sep = "\t", fill = NA, col.names = paste0("V", 1:20))
    LogTmp <- read.table(paste0(Path, McmcLogs[i]), sep = "\t", header = TRUE)
    for(y in 1:nrow(SpTmp)){
      # Speciation
      Sp1Tmp <- unlist( SpTmp[y, ] )
      Kbirth <- LogTmp[y, "k_birth"]
      # Rarely there is a shift inferred earlier then 13.65
      Shifts <- Sp1Tmp[(Kbirth+1):length(Sp1Tmp)]
      Shifts <- na.omit(Shifts)
      if (all(Shifts <= EdgeShift)) {
        SpThroughTime[Counter, ] <- approx(x = c(Shifts, 0), y = Sp1Tmp[1:Kbirth], 
                                           xout = TimeSeq, method = "constant")$y
        SpShifts[Counter, 1:length(Shifts)] <- Shifts
      }
      # Extinction
      Ex1Tmp <- unlist( ExTmp[y, ] )
      Kdeath <- LogTmp[y, "k_death"]
      Shifts <- Ex1Tmp[(Kdeath+1):length(Sp1Tmp)]
      Shifts <- na.omit(Shifts)
      if (all(Shifts <= EdgeShift)) {
        ExThroughTime[Counter, ] <- approx(x = c(Shifts, 0), y = Ex1Tmp[1:Kdeath], 
                                           xout = TimeSeq, method = "constant")$y
        ExShifts[Counter, 1:length(Shifts)] <- Shifts
      }
      Counter <- Counter + 1
    }
  }
  # Omit the fixed edge shift of 13.65 million years
  SpThroughTime <- SpThroughTime[, -1]
  ExThroughTime <- ExThroughTime[, -1]
  SpShifts <- SpShifts[, -1]
  ExShifts <- ExShifts[, -1]
  # Omit all rows with occasional NAs 
  Keep <- 1:(Counter-1) # Easier
  SpThroughTime <- SpThroughTime[Keep, ]
  ExThroughTime <- ExThroughTime[Keep, ]
  SpShifts <- SpShifts[Keep, ]
  ExShifts <- ExShifts[Keep, ]
  # Rescale rates from 13.6 million years to 1.36
  SpThroughTime <- SpThroughTime * 10
  ExThroughTime <- ExThroughTime * 10
  NetDivThroughTime <- SpThroughTime - ExThroughTime
  SpShifts <- SpShifts / 10
  ExShifts <- ExShifts / 10
  Res <- list()
  # Summarize shifts
  ###########################################################################
  # Fequency of shifts found across all birth-death analyses
  PerSpShift <- apply(SpShifts, 2, function(x) sum(!is.na(x)) / length(x))
  PerExShift <- apply(ExShifts, 2, function(x) sum(!is.na(x)) / length(x))
  # 
  MedSpShift <- apply(SpShifts, 2, median, na.rm = TRUE) 
  MedExShift <- apply(ExShifts, 2, median, na.rm = TRUE) 
  
  ExShiftsHpd <- SpShiftsHpd <- matrix(NA_real_, nrow = 2, ncol = ncol(SpShifts))
  for(i in 1:ncol(SpShifts)){
    TmpShift <- na.omit(SpShifts[, i])
    if(length(TmpShift) > 1){
      SpShiftsHpd[, i] <- HPDinterval(coda:::mcmc(TmpShift))
    }
    TmpShift <- na.omit(ExShifts[, i])
    if(length(TmpShift) > 1){
      ExShiftsHpd[, i] <- HPDinterval(coda:::mcmc(TmpShift))
    }
  }
  
  SummarySpShift <- rbind(PerSpShift, MedSpShift, SpShiftsHpd) 
  SummarySpShift <- SummarySpShift[, SummarySpShift[1, ] != 0] 
  rownames(SummarySpShift) <- c("Frequency", "Age", "LwrHPD", "UprHPD")
  colnames(SummarySpShift) <- paste("Shift", 1:ncol(SummarySpShift))
  Res[[1]] <- SummarySpShift
  
  SummaryExShift <- rbind(PerExShift, MedExShift, ExShiftsHpd) 
  SummaryExShift <- SummaryExShift[, SummaryExShift[1, ] != 0] 
  rownames(SummaryExShift) <- c("Frequency", "Age", "LwrHPD", "UprHPD")
  colnames(SummaryExShift) <- paste("Shift", 1:ncol(SummaryExShift))
  Res[[2]] <- SummaryExShift
  Res[[3]] <- SpShifts
  Res[[4]] <- ExShifts
  #######################################################################
  RateHPDs <- getRatesHPD(SpThroughTime,
                          ExThroughTime,
                          NetDivThroughTime,
                          Prob)
  Res[[5]] <- RateHPDs
  Res[[6]] <- SpThroughTime
  Res[[7]] <- SpThroughTime
  Res[[8]] <- NetDivThroughTime
  #######################################################################
  names(Res) <- c("SummarySpShifts", "SummaryExShifts",
                  "SpShifts", "ExShifts",
                  "RatesHPDs",
                  "SpThroughTime", "ExThroughTime", "NetDivThroughTime")
  return(Res)
}