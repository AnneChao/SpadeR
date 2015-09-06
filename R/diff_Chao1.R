diff_Chao1 <-
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
  if (s == 1){
    d <- 1
  } else if (s == 2) {
    d <- f1p * fp1 / (4 * f2p * fp2)
  } else if (s == 3) {
    d <- f11 * f1p / (4 * f2p * fp2) + fp1 / fp2
  } else if (s == 4) {
    d <- f11 * fp1 / (4 * f2p * fp2) + f1p / f2p
  } else if (s == 5) {
    d <- - f11 * f1p * fp1 / (4 * f2p) / fp2^2 - (fp1 / fp2)^2 / 2
  } else if (s == 6) {
    d <- - f11 * f1p * fp1 / (4 * fp2) / f2p^2 - (f1p / f2p)^2 / 2
  } else {
    d <- 0
  }
  return(d)
}
