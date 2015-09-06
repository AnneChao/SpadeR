Chao2_sharedFun <-
function(y1, y2, conf=0.95) {
  x1 <- y1[-1]
  x2 <- y2[-1]
  f2p <- sum(x1 == 2 & x2 >= 1)
  fp2 <- sum(x1 >= 1 & x2 == 2)
  
  est <- Chao2_sharedEstFun(y1, y2)
  if (f2p == 0 || fp2 == 0) {
    se <- VarEstFun.Sam(y1, y2, diffFun=diff_Chao2bc, FunName=Chao2_bcEstFun)
  } else {
    se <- VarEstFun.Sam(y1, y2, diffFun=diff_Chao2, FunName=Chao2_sharedEstFun)
  }
  CI <- logCI(y1, y2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Chao2-shared")
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
