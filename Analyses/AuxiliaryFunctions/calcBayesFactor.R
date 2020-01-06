calcBayesFactor <- function(Lik){
  MaxLl <- max(Lik, na.rm = TRUE)
  Dif <- Lik - MaxLl
  Res <- 2*-Dif
  return(Res)
}