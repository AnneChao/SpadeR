diff_Chao1bc <-
function(ksi, s) {
  D12 <- ksi[1]
  f11 <- ksi[2]
  fp1 <- ksi[3]
  f1p <- ksi[4]
  fp2 <- ksi[5]
  f2p <- ksi[6]
  if (s == 1){
    d <- 1
  } else if (s == 2) {
    d <- f1p * fp1 / (4 * (f2p + 1) * (fp2 + 1))
  } else if (s == 3) {
    d <- f11 * f1p / (4 * (f2p + 1) * (fp2 + 1)) + (2 * fp1 - 1)/ (2 * (fp2 + 1))
  } else if (s == 4) {
    d <- f11 * fp1 / (4 * (f2p + 1) * (fp2 + 1)) + (2 * f1p - 1)/ (2 * (f2p + 1))
  } else if (s == 5) {
    d <- - f11 * f1p * fp1 / (4 * (f2p + 1)) / (fp2 + 1)^2 - fp1 * (fp1 - 1) / 2 / (fp2 + 1)^2
  } else if (s == 6) {
    d <- - f11 * f1p * fp1 / (4 * (fp2 + 1)) / (f2p + 1)^2 - f1p * (f1p - 1) / 2 / (f2p + 1)^2
  } else {
    d <- 0
  }
  return(d)
}
