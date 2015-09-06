Chao2_sharedEstFun <-
function(y1, y2) {
  t1 <- y1[1]
  t2 <- y2[1]
  c1 <- (t1 - 1) / t1
  c2 <- (t2 - 1) / t2
  x1 <- y1[-1]
  x2 <- y2[-1]
  D12 <- sum(x1 > 0 & x2 > 0)  
  Q11 <- sum(x1 == 1 & x2 == 1)
  Q1.plus <- sum(x1 == 1 & x2 >= 1)
  Qplus.1 <- sum(x2 == 1 & x1 >= 1)
  Q2.plus <- sum(x1 == 2 & x2 >= 1)
  Qplus.2 <- sum(x2 == 2 & x1 >= 1)  
  if (Q2.plus == 0 || Qplus.2 == 0) {
    est <- D12 + Q11 * c1 * c2 * Q1.plus * Qplus.1 / (4 * (Q2.plus + 1) * (Qplus.2 + 1)) + 
      c1 * Q1.plus * (Q1.plus - 1) / (2 * (Q2.plus + 1)) + 
      c2 * Qplus.1 * (Qplus.1 - 1) / (2 * (Qplus.2 + 1))
  } else {
    est <- D12 + Q11 * c1 * c2 * Q1.plus * Qplus.1 / (4 * Q2.plus * Qplus.2) + 
      c1 * Q1.plus^2 / (2 * Q2.plus) + c2 * Qplus.1^2 / (2 * Qplus.2)
  }
  return(est)
}
