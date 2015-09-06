PanFun <-
function(x1, x2, conf=0.95) {
  f22 <- sum(x1 == 2 & x2 == 2)
  f2p <- sum(x1 == 2 & x2 >= 1)
  fp2 <- sum(x1 >= 1 & x2 == 2)
  est <- PanEstFun(x1, x2)
  if (f2p == 0 || fp2 == 0 || f22 == 0) {
    se <- VarEstFun(x1, x2, diffFun=diff_Panbc, FunName=PanbcEstFun)
  } else {
    se <- VarEstFun(x1, x2, diffFun=diff_Pan, FunName=PanEstFun)
  }
  CI <- logCI(x1, x2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Lower-bound")
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
