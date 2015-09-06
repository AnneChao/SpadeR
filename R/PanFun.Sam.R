PanFun.Sam <-
function(y1, y2, conf=0.95) {
  x1 <- y1[-1]
  x2 <- y2[-1]
  f22 <- sum(x1 == 2 & x2 == 2)
  f2p <- sum(x1 == 2 & x2 >= 1)
  fp2 <- sum(x1 >= 1 & x2 == 2)
  est <- PanEstFun.Sam(y1, y2)
  if (f2p == 0 || fp2 == 0 || f22 == 0) {
    se <- VarEstFun.Sam(y1, y2, diffFun=diff_Panbc, FunName=PanbcEstFun.Sam)
  } else {
    se <- VarEstFun.Sam(y1, y2, diffFun=diff_Pan, FunName=PanEstFun.Sam)
  }
  CI <- logCI(y1, y2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Lower-bound")
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
