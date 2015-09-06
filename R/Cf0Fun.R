Cf0Fun <-
function(x) {
  n <- sum(x)
  f1 <- sum(x == 1); f2 <- sum(x == 2)
  if (f2 > 0) {
    C <- 1 - f1 / n * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))
    f0 <- (n - 1) / n * f1^2 / (2 * f2)
  } else if (f2 == 0 & f1 != 0) {
    C <- 1 - f1 / n * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
  } else {
    C <- 1
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
  }
  f0 <- ceiling(f0)
  return(c(C, f0))
}
