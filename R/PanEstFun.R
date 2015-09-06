PanEstFun <-
function(x1, x2) {
  n1 <- sum(x1)
  n2 <- sum(x2)
  D12 <- sum(x1 > 0 & x2 > 0)
  f11 <- sum(x1 == 1 & x2 == 1)
  f22 <- sum(x1 == 2 & x2 == 2)
  f1p <- sum(x1 == 1 & x2 >= 1)
  fp1 <- sum(x1 >= 1 & x2 == 1)
  f2p <- sum(x1 == 2 & x2 >= 1)
  fp2 <- sum(x1 >= 1 & x2 == 2)
  K1 <- (n1 - 1) / n1
  K2 <- (n2 - 1) / n2
  if (f2p == 0 || fp2 == 0 || f22 == 0) {
    est <- D12 + K1 * f1p * (f1p - 1) / 2 / (f2p + 1) + 
      K2 * fp1 * (fp1 - 1) / 2 / (fp2 + 1) + 
      K1 * K2 * f11 * (f11 - 1) / 4 / (f22 + 1)
  } else {
    est <- D12 + K1 * f1p^2 / 2 / f2p + K2 * fp1^2 / 2 / fp2 + 
      K1 * K2 * f11^2 / 4 / f22
  }
  return(est)
}
