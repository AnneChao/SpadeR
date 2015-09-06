VarEstFun <-
function(x1, x2, diffFun, FunName) {
  n1 <- sum(x1)
  n2 <- sum(x2)
  D12 <- sum(x1 > 0 & x2 > 0)
  f11 <- sum(x1 == 1 & x2 == 1)
  f22 <- sum(x1 == 2 & x2 == 2)
  f12 <- sum(x1 == 1 & x2 == 2)
  f21 <- sum(x1 == 2 & x2 == 1)
  f1p <- sum(x1 == 1 & x2 >= 1)
  fp1 <- sum(x1 >= 1 & x2 == 1)
  f2p <- sum(x1 == 2 & x2 >= 1)
  fp2 <- sum(x1 >= 1 & x2 == 2)
  K1 <- (n1 - 1) / n1
  K2 <- (n2 - 1) / n2
  ksi <- c(D12, f11, fp1, f1p, fp2, f2p, f22, K1, K2)
  
  S12 <- FunName(x1, x2)
  CovFun <- function(i, j) {
    if (i == j) {
      cov <- ksi[i] - ksi[i]^2 / S12
    } else if ((i == 1 & j == 2) | (j == 1 & i == 2)) {
      cov <- f11 - ksi[i] * ksi[j] / S12
    } else if ((i == 1 & j == 3) | (j == 1 & i == 3)) {
      cov <- fp1 - ksi[i] * ksi[j] / S12
    } else if ((i == 1 & j == 4) | (j == 1 & i == 4)) {
      cov <- f1p - ksi[i] * ksi[j] / S12
    } else if ((i == 1 & j == 5) | (j == 1 & i == 5)) {
      cov <- fp2 - ksi[i] * ksi[j] / S12
    } else if ((i == 1 & j == 6) | (j == 1 & i == 6)) {
      cov <- f2p - ksi[i] * ksi[j] / S12
    } else if ((i == 2 & j == 3) | (j == 2 & i == 3)) {
      cov <- f11 - ksi[i] * ksi[j] / S12
    } else if ((i == 2 & j == 4) | (j == 2 & i == 4)) {
      cov <- f11 - ksi[i] * ksi[j] / S12
    } else if ((i == 3 & j == 4) | (j == 3 & i == 4)) {
      cov <- f11 - ksi[i] * ksi[j] / S12
    } else if ((i == 3 & j == 6) | (j == 3 & i == 6)) {
      cov <- f21 - ksi[i] * ksi[j] / S12
    } else if ((i == 4 & j == 5) | (j == 4 & i == 5)) {
      cov <- f12 - ksi[i] * ksi[j] / S12
    } else if ((i == 5 & j == 6) | (j == 5 & i == 6)) {
      cov <- f22 - ksi[i] * ksi[j] / S12
    } else if ((i == 1 & j == 7) | (j == 1 & i == 7)) {  #  Pan has 7 parameter
      cov <- f22 - ksi[i] * ksi[j] / S12
    } else if ((i == 5 & j == 7) | (j == 5 & i == 7)) {  #  Pan has 7 parameter
      cov <- f22 - ksi[i] * ksi[j] / S12
    } else if ((i == 6 & j == 7) | (j == 6 & i == 7)) {  #  Pan has 7 parameter
      cov <- f22 - ksi[i] * ksi[j] / S12
    } else {
      cov <- 0 - ksi[i] * ksi[j] / S12
    }
    return(cov)
  }
  
  diff <- diffFun
  i <- rep(c(1:7), 7)
  j <- rep(c(1:7), each=7)
  var <- sum(mapply(function(i, j) diff(ksi, i) * diff(ksi, j) * CovFun(i, j), i, j))
  se <- sqrt(var)
  return(se)
}
