Chao1_bcEstFun <-
function(x1, x2) {
  D12 <- sum(x1 > 0 & x2 > 0)
  x1_share <- x1[which(x1 > 0 & x2 > 0)]
  x2_share <- x2[which(x1 > 0 & x2 > 0)]
  f11 <- sum(x1_share == 1 & x2_share == 1)
  f1.plus <- sum(x1_share == 1 & x2_share >= 1)
  fplus.1 <- sum(x2_share == 1 & x1_share >= 1)
  f2.plus <- sum(x1_share == 2 & x2_share >= 1)
  fplus.2 <- sum(x2_share == 2 & x1_share >= 1)
  est <- D12 + f11 * f1.plus * fplus.1 / (4 * (f2.plus + 1) * (fplus.2 + 1)) + 
    f1.plus * (f1.plus - 1) / (2 * (f2.plus + 1)) + 
    fplus.1 * (fplus.1 - 1) / (2 * (fplus.2 + 1))
  return(est)
}
