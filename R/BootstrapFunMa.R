BootstrapFunMa <-
function(x1, x2, B, FunName) {
  n1 <- sum(x1); n2 <- sum(x2)
  z <- SortDataFun(x1, x2)
  z1 <- z[, 1]
  z2 <- z[, 2]
  newprob <- ExtenProbFunMa(z1, z2)
  p1 <- newprob$prob1
  p2 <- newprob$prob2
  set.seed(123)
  X1 <- rmultinom(B, n1, p1)
  set.seed(123)
  X2 <- rmultinom(B, n2, p2)
  X <- rbind(X1, X2)
  
  se <- sd(apply(X, 2, function(x) {
    y <- matrix(x, ncol=2)
    y1 <- y[, 1]
    y2 <- y[, 2]
#     y1 <- x[1 : length(p1)]
#     y2 <- x[(length(p1) + 1) : (2 * length(p1))]
    FunName(y1, y2)
  }), na.rm=T)
  
  return(se)
}
