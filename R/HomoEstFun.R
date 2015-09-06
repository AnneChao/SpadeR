HomoEstFun <-
function(x1, x2) {
  D1 <- sum(x1 > 0)
  D2 <- sum(x2 > 0)
  D12 <- sum(x1 > 0 & x2 > 0)
  x1_share <- x1[which(x1 > 0 & x2 > 0)]
  x2_share <- x2[which(x1 > 0 & x2 > 0)]
  f11 <- sum(x1_share == 1 & x2_share == 1)
  
  D12_rare <- sum(x1_share <= 10 & x2_share <= 10)
  
  pos_r <- (x1_share <= 10 & x2_share <= 10)
  pos1_r <- (x1_share == 1 & x2_share <= 10)
  pos2_r <- (x2_share == 1 & x1_share <= 10)
  
  tmp <- sum(x1_share[pos_r] * x2_share[pos_r])
  if (tmp == 0)  # correct when number of Xi * Yi = 0 equal to 0 by Y.H. Lee
    tmp <- 1
  
  C12_rare <- 1 - (sum(x2_share[pos1_r]) + sum(x1_share[pos2_r]) - f11) / tmp
  if (C12_rare == 0 || C12_rare > 1)  # Correct when C12 = 0 or C12 > 1 !!! by c++
    C12_rare <- 1
  est <- D12 - D12_rare + D12_rare / C12_rare
  return(est)
}
