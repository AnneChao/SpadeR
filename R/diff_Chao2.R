diff_Chao2 <-
function(ksi, s) {
  if (sum(ksi == 0) != 0) {
    ksi[which(ksi == 0)] <- 1
  }
  D12 <- ksi[1]
  f11 <- ksi[2]
  fp1 <- ksi[3]
  f1p <- ksi[4]
  fp2 <- ksi[5]
  f2p <- ksi[6]
  K1 <- ksi[8]
  K2 <- ksi[9]
  if (s == 1){
    d <- 1
  } else if (s == 2) {
    d <- K1 * K2 * f1p * fp1 / (4 * f2p * fp2)
  } else if (s == 3) {
    d <- K1 * K2 * f11 * f1p / (4 * f2p * fp2) + K2 * fp1 / fp2
  } else if (s == 4) {
    d <- K1 * K2 * f11 * fp1 / (4 * f2p * fp2) + K1 * f1p / f2p
  } else if (s == 5) {
    d <- - K1 * K2 * f11 * f1p * fp1 / (4 * f2p) / fp2^2 - K2 * (fp1 / fp2)^2 / 2
  } else if (s == 6) {
    d <- - K1 * K2 * f11 * f1p * fp1 / (4 * fp2) / f2p^2 - K1 * (f1p / f2p)^2 / 2
  } else {
    d <- 0
  }
  return(d)
}
