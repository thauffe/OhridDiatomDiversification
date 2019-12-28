finitDeriv <- function(Gam, NewData, H = 0.001, Scale = TRUE){
  # Gam: gam for the time series
  # NewData: Time for which the derivative should be calculated
  # H: Difference, here change per year. But this is for numerical precision only. 
  #                                      The output will be in ka.
  PlusH <- predict(Gam, newdata = NewData + H, type = "response", se.fit = FALSE)
  MinusH <- predict(Gam, newdata = NewData - H, type = "response", se.fit = FALSE)
  FirstDeriv <- (PlusH - MinusH) / 2*H
  FirstDeriv <- FirstDeriv*1000 # Change per 1000 years
  if(Scale){
    FirstDeriv <- scale(FirstDeriv)
  }
  Res <- data.frame(NewData$Age, FirstDeriv)
  return(Res)
}