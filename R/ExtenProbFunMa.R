ExtenProbFunMa <-
function(z1, z2) {
  x1 <- z1; x2 <- z2  # Sorted data
  n1 <- sum(x1); n2 <- sum(x2)
  D1 <- sum(x1 > 0); D2 <- sum(x2 > 0)
  D12 <- sum(x1 > 0 & x2 > 0)
  f11 <- sum(x1 == 1 & x2 == 1)
  f22 <- sum(x1 == 2 & x2 == 2)
  f1p <- sum(x1 == 1 & x2 >= 1)
  fp1 <- sum(x1 >= 1 & x2 == 1)
  f2p <- sum(x1 == 2 & x2 >= 1)
  fp2 <- sum(x1 >= 1 & x2 == 2)
  K1 <- (n1 - 1) / n1; K2 <- (n2 - 1) / n2
  f0p <- ifelse(f2p == 0,
                f0p <- ceiling(K1 * f1p * (f1p - 1) / 2 / (f2p + 1)), 
                f0p <- ceiling(K1 * f1p^2 / 2 / f2p))
  fp0 <- ifelse(fp2 == 0, 
                fp0 <- ceiling(K2 * fp1 * (fp1 - 1) / 2 / (fp2 + 1)),
                fp0 <- ceiling(K2 * fp1^2 / 2 / fp2 ))
  f00 <- ifelse(f22 == 0,
                f00 <- ceiling(K1 * K2 * f11 * (f11 - 1) / 4 / (f22 + 1)),
                f00 <- ceiling(K1 * K2 * f11^2 / 4 / f22))
  
  ga1 <- D1 - D12; ga2 <- D2 - D12
  #   d1 <- max(f0p, ga2)  # maximum of f0+ & gamma2
  #   d2 <- max(fp0, ga1)  # maximum of f+0 & gamma1
  
  
  ## Community 1
  tmp1 <- Cf0Fun(x1)
  Chat1 <- tmp1[1] ; f0_1 <- tmp1[2]
  more1 <- max(f0_1, (max(ga1, fp0) - ga1 + f0p + f00))  # unseen species in commuity 1
  add1 <- max(f0_1 - max(ga1, fp0) + ga1 - f0p - f00, 0)  # unseen endemic in commuity 1
  lambda1 <- (1 - Chat1) / sum(x1 / n1 * (1 - x1 / n1)^n1)
  pi1 <- x1 / n1 * (1 - lambda1 * (1 - x1 /n1)^n1)
  p0_1 <- (1 - Chat1) / more1
  
  ## Community 2
  tmp2 <- Cf0Fun(x2)
  Chat2 <- tmp2[1] ; f0_2 <- tmp2[2]
  more2 <- max(f0_2, (max(ga2, f0p) - ga2 + fp0 + f00))  # unseen species in commuity 2
  add2 <- max(f0_2 - max(ga2, f0p) + ga2 - fp0 - f00, 0)  # unseen endemic in commuity 2
  lambda2 <- (1 - Chat2) / sum(x2 / n2 * (1 - x2 / n2)^n2)
  pi2 <- x2 / n2 * (1 - lambda2 * (1 - x2 /n2)^n2)
  p0_2 <- (1 - Chat2) / more2
  
  ## Extension probility for Community 1
  prob1 <- c(pi1[pi1 > 0], rep(p0_1, max(fp0-ga1, 0)), rep(p0_1, f0p),
             rep(0, max(ga2-f0p, 0)), rep(p0_1, f00), rep(p0_1, add1), 
             rep(0, add2))
  
  ## Extension probility for Community 2
  prob2 <- c(pi2[1:D12], rep(p0_2, fp0), rep(0, max(ga1-fp0, 0)), 
             pi2[(D12 + ga1 + 1):(D12 + ga1 + ga2)], rep(p0_2, max(f0p-ga2, 0)),
             rep(p0_2, f00), rep(0, add1), rep(p0_2, add2))
  out <- list(prob1=prob1, prob2=prob2)
  return(out)
}
