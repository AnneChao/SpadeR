PanbcFun <-
function(x1, x2, conf=0.95) {
  est <- PanbcEstFun(x1, x2)
  se <- VarEstFun(x1, x2, diffFun=diff_Panbc, FunName=PanbcEstFun)
  CI <- logCI(x1, x2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Lower-bound-bc")
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
