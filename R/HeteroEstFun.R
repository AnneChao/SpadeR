HeteroEstFun <-
function(x1, x2) {
  n1 <- sum(x1)
  n2 <- sum(x2)
  D1 <- sum(x1 > 0)
  D2 <- sum(x2 > 0)
  D12 <- sum(x1 > 0 & x2 > 0)
  x1_share <- x1[which(x1 > 0 & x2 > 0)]
  x2_share <- x2[which(x1 > 0 & x2 > 0)]
  f11 <- sum(x1_share == 1 & x2_share == 1)
  f1.plus <- sum(x1_share == 1 & x2_share <= 10)
  fplus.1 <- sum(x2_share == 1 & x1_share <= 10)
  f2.plus <- sum(x1_share == 2 & x2_share <= 10)
  fplus.2 <- sum(x2_share == 2 & x1_share <= 10)
  D12_rare <- sum(x1_share <= 10 & x2_share <= 10)
  
  pos <- (x1 > 0 & x2 > 0) & (x1 > 10 | x2 > 10)
  n1_rare <- n1 - sum(x1[pos])
  n2_rare <- n2 - sum(x2[pos])
  
  # correct when n1_rare = 1 or n2_rare = 0 by Y.H. Lee
  if (n1_rare == 1) 
    n1_rare <- 2
  if (n2_rare == 1) 
    n2_rare <- 2
  
  pos_r <- (x1_share <= 10 & x2_share <= 10)
  pos1_r <- (x1_share == 1 & x2_share <= 10)
  pos2_r <- (x2_share == 1 & x1_share <= 10)
  
  tmp <- sum(x1_share[pos_r] * x2_share[pos_r])
  if (tmp == 0)  # correct when number of Xi * Yi = 10 equal to 0 by Y.H. Lee
    tmp <- 1
  
  C12_rare <- 1 - (sum(x2_share[pos1_r]) + sum(x1_share[pos2_r]) - f11) / tmp
  if (C12_rare == 0 || C12_rare > 1)  # Correct when C12 = 0 or C12 > 1 !!! by c++
    C12_rare <- 1
  #   C12_rare <- round(C12_rare, 4)
  
  T10 <- sum(x1_share[x1_share <= 10 & x2_share <= 10])
  T01 <- sum(x2_share[x1_share <= 10 & x2_share <= 10])
  T11 <- tmp
  T21 <- sum(x1_share[pos_r] * (x1_share - 1)[pos_r] * x2_share[pos_r])
  T12 <- sum(x1_share[pos_r] * (x2_share - 1)[pos_r] * x2_share[pos_r])
  
  T22 <- sum(x1_share[pos_r] * x2_share[pos_r] * 
               (x1_share - 1)[pos_r] * (x2_share - 1)[pos_r])
  
  if (T10 == 0)  # correct when equal to 0 by Y.H. Lee
    T10 <- 1
  if (T11 == 0)
    T11 <-1
  if (T01 == 0)
    T01 <- 1
  
  S12_0 <- D12_rare / C12_rare
  CCV_1 <- S12_0 * n1_rare * T21 / (n1_rare - 1) / T10 / T11 - 1
  CCV_2 <- S12_0 * n2_rare * T12 / (n2_rare - 1) / T01 / T11 - 1
  CCV_12 <- n1_rare * n2_rare * S12_0^2 * T22 / 
    ((n1_rare - 1) * (n2_rare - 1) * T10 * T01 * T11) - 
    S12_0 * T11 / T10 / T01 - CCV_1 - CCV_2
  
  tmp1 <- D12 - D12_rare + D12_rare / C12_rare 
  tmp2 <- (f1.plus * CCV_1 + fplus.1 * CCV_2 + f11 * CCV_12) / C12_rare
  est <- tmp1 + tmp2
  return(est)
}
