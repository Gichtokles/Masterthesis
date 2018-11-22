

mlag <- function(X, lag)
{
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(0, Traw, p * N)
  for (ii in 1:p) {
    Xlag[(ii+1):Traw, (N * (ii - 1) + 1):(N * ii)] = X[2:(Traw-ii+1), (1:N)]
  }
  return(Xlag)
  
}

X  <- data
lag <- 10

for(ii in 1:p){
  print(length(2:(Traw-ii+1)))
}