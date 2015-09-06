logCI <-
function(y1, y2, est, se, conf) {
  x1 <- y1[-1]
  x2 <- y2[-1]
  D12 <- sum(x1 > 0 & x2 > 0)
  t <- est - D12
  z <- qnorm((1-conf)/2, lower.tail=F)
  K <- exp(z * sqrt(log(1 + se^2 / t^2)))
  CI <- c(D12 + t / K, D12 + t * K)
  return(CI)
}
