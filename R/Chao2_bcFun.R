Chao2_bcFun <-
function(y1, y2, conf=0.95) {
  est <- Chao2_bcEstFun(y1, y2)
  se <- VarEstFun.Sam(y1, y2, diffFun=diff_Chao2bc, FunName=Chao2_bcEstFun)
  CI <- logCI(y1, y2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Chao2-shared-bc")
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
