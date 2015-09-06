HeteroFun <-
function(x1, x2, B, conf=0.95) {
  est <- HeteroEstFun(x1, x2)
  se <- BootstrapFunMa(x1, x2, B, FunName=HeteroEstFun)
  CI <- logCI(x1, x2, est, se, conf)
  out <- matrix(c(est, se, CI), nrow = 1)
  out <- data.frame(out)
  rownames(out) <- c("Heterogeneous(ACE-shared)")
  
  colnames(out) <- c("Estimator", "Est_s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}
