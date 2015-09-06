SortDataFun <-
function(x1, x2) {
  D1 <- sum(x1 > 0); D2 <- sum(x2 > 0)
  D12 <- sum(x1 > 0 & x2 > 0)
  sit <- order(x1 > 0 & x2 > 0, decreasing=T)
  #   set.seed(123)
  #   sit <- sample(sit[1:D12], replace=F)  # random sample
  common <- cbind(x1[sit[1:D12]], x2[sit[1:D12]])
  sit1 <- order(x1 > 0 & x2 == 0, decreasing=T)
  set.seed(123)
  sit1 <- sample(sit1[1:(D1-D12)], replace=F)  # random sample
  special1 <- cbind(x1[sit1[1:(D1-D12)]], x2[sit1[1:(D1-D12)]])
  sit2 <- order(x1 == 0 & x2 > 0, decreasing=T)
  set.seed(123)
  sit2 <- sample(sit2[1:(D2-D12)], replace=F)  # random sample
  special2 <-  cbind(x1[sit2[1:(D2-D12)]], x2[sit2[1:(D2-D12)]])
  z <- rbind(common, special1, special2)
  return(z)
}
