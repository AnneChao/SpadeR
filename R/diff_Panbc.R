diff_Panbc <-
function(ksi, s) {
  D12 <- ksi[1]
  f11 <- ksi[2]
  fp1 <- ksi[3]
  f1p <- ksi[4]
  fp2 <- ksi[5]
  f2p <- ksi[6]
  f22 <- ksi[7]
  K1 <- ksi[8]
  K2 <- ksi[9]
  if (s == 1){
    d <- 1
  } else if (s == 2) {
    d <- K1 * K2 * (2 * f11 - 1) / (4 * (f22 + 1))
  } else if (s == 3) {
    d <- K2 * (2 * fp1 -1) / (2 * (fp2 + 1))
  } else if (s == 4) {
    d <- K1 * (2 * f1p - 1) / (2 * (f2p + 1))
  } else if (s == 5) {
    d <- - K2 * fp1 * (fp1 - 1) / (fp2 + 1)^2 / 2
  } else if (s == 6) {
    d <- - K1 * f1p * (f1p - 1) / (f2p + 1)^2 / 2
  } else {
    d <- - K1 * K2 * f11 * (f11 - 1) / 4 / (f22 + 1)^2
  }
  return(d)
}
