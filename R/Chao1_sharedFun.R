Chao1_sharedFun <-
function(x1, x2, conf=0.95) {
  f2p <- sum(x1 == 2 & x2 >= 1)
  fp2 <- sum(x1 >= 1 & x2 == 2)
  
  est <- Chao1_sharedEstFun(x1, x2)
  if (f2p == 0 || fp2 == 0) { 
    se <- VarEstFun(x1, x2, diffFun=diff_Chao1bc, FunName=Chao1_bcEstFun)
  } else {
    se <- VarEstFun(x1, x2, diff_Chao1, FunName=Chao1_sharedEstFun)
  }
  
  CI <- logCI(x1, x2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Chao1(shared)")
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
