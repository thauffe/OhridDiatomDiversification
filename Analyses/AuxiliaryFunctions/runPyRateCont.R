runPyRateCont <- function(i, Comb, Path){
  ArgsConstDD <- ""
  # Covariate
  ArgsCovar <- paste0("-c Results/Covariates/", Comb[i, "Covariate"], ".txt") 
  if (Comb[i, "Covariate"] == "DD"){
    # Get non-endemic diversity
    Div <- read.table("Results/Covariates/NonEnd.txt", 
                      sep = "\t", header = TRUE)
    Div <- Div[nrow(Div):1, ] # Present to past
    # Get endemic diversity
    FixTmp <- read.table(paste0("Results/BD/", Comb[i, "Fix"], ".txt"),
                         sep = "\t", header = TRUE)
    FirstApp <- sort(FixTmp$ts, decreasing = TRUE)
    Increase <- rep(1, length(FirstApp))
    LastApp <- sort(FixTmp$te, decreasing = TRUE)
    Decrease <- rep(-1, length(LastApp))
    Decrease[LastApp == 0] <- 0
    Att <- data.frame(Time = c(FirstApp, LastApp), Change = c(Increase, Decrease), SR = NA_integer_)
    Att <- Att[order(Att$Time, decreasing = TRUE), ]
    Att$SR <- cumsum(Att$Change)
    Att <- Att[min(which(Att$Time <= 13.66)):which(Att$Time == 0)[1], ]
    End <- approx(Att[, c(1, 3)], xout = Div$time, method ="constant", rule = 2)$y
    Div$Div <- Div$NonEnd + End
    AddTime <- seq(13.67, 50.00, length.out = 100)
    WriteDiv <- data.frame(time = c(Div$time, AddTime),
                           SR = c(scale(Div$Div), rep(0, length(AddTime))))
    ArgsConstDD <- "-constDD"
    colnames(WriteDiv)[2] <- "DD"
    DdName <- paste0("DD_", Comb[Args[1], "Fix"],".txt")
    write.table(WriteDiv, 
                file = paste0("Results/Covariates/", DdName), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    ArgsCovar <- paste0("-c Results/Covariates/", DdName) 
  }
  system( paste("python2 PyRate/PyRateContinuous2.py",
                paste0("-d Results/", Path,"/", Comb[i, "Fix"], ".txt"),
                ArgsCovar, # Covariate
                "-m", Comb[i, "Model"],
                "-n 10", 
                "-s 1", 
                "-p 10", 
                "-stimesL", Comb[i, "stimesL"],
                "-stimesM", Comb[i, "stimesM"],
                ArgsConstDD, 
                "-r 1",
                "-A 1")) # TI
}