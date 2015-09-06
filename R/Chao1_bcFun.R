Chao1_bcFun <-
function(x1, x2, conf=0.95) {
  est <- Chao1_bcEstFun(x1, x2)
  se <- VarEstFun(x1, x2, diffFun=diff_Chao1bc, FunName=Chao1_bcEstFun)
  CI <- logCI(x1, x2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Chao1-shared-bc")
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
