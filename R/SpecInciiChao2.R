SpecInciiChao2 <- function(data, k, conf){
  data <- as.numeric(data)
  z <- -qnorm((1 - conf)/2)
  t <- data[1]
  dat <- data[-1]
  x <- dat[which(dat != 0)]
  Q <- function(i, data){length(data[which(data == i)])}
  
  basicInci <- function(data, k){
    data <- as.numeric(data)
    t <- data[1]
    dat <- data[-1]
    x <- dat[which(dat != 0)]
    Q <- function(i, data){length(data[which(data == i)])}
    
    D <- length(x)
    D_infreq <- length(x[which(x <= k)])
    
    if (Q(1, x) > 0 & Q(2, x) > 0){
      A <- 2*Q(2, x)/((t-1)*Q(1, x) + 2*Q(2, x))
    } else if (Q(1, x) > 0 & Q(2, x) == 0){
      A <- 2/((t-1)*(Q(1, x) - 1) + 2)
    } else {
      A <- 1
    }
    C_infreq <- 1 - Q(1, x)/sum(x[which(x <= k)])*(1-A)
    
    j <- c(1:k)
    b1 <- sum(sapply(j, function(j)j*(j-1)*Q(j, x)))
    b2 <- sum(sapply(j, function(j)j*Q(j, x)))
    gamma_infreq_square <- max(D_infreq/C_infreq*t/(t - 1)*b1/b2/(b2) - 1, 0)
    CV_infreq <- sqrt(gamma_infreq_square)
    D_freq <- length(x[which(x > k)])
    
    BASIC.DATA <- matrix(paste(c("D", "t", "k", "D_infreq", "C_infreq", "CV_infreq", "D_freq"),
                               c(D,t,k,D_infreq,C_infreq,CV_infreq,D_freq),
                               sep = "="), ncol=1)
    colnames(BASIC.DATA)=c("Value")
    rownames(BASIC.DATA)=c("Number of observed species","Number of sample/quadrats","Cut-off point",
                           "Number of observed species for infrequent species","Estimated sample coverage for infrequent species",
                           "Estimated CV for infrequent species",
                           "Number of observed species for frequent species")
    return(list(BASIC.DATA, D, t, D_infreq, C_infreq, CV_infreq, D_freq))
  }
  D <- basicInci(data, k)[[2]]
  D_infreq <- basicInci(data, k)[[4]]
  C_infreq <- basicInci(data, k)[[5]]
  CV_infreq <- basicInci(data, k)[[6]]
  D_freq <- basicInci(data, k)[[7]]
  
  q1 <- Q(1, x); q2 <- Q(2, x); q3 <- Q(3, x); q4 <- Q(4, x)
  if (q1 > 0 & q2 != 0){
    s_Chao2 <- D + (t - 1)/t*q1^2/(2*q2)
    var_Chao2 <- q2*((t - 1)/t*(q1/q2)^2/2 + ((t - 1)/t)^2*(q1/q2)^3 + ((t - 1)/t)^2*(q1/q2)^4/4)
  } else if (q1 > 1 & q2 == 0){
    s_Chao2 <- D + (t - 1)/t*q1*(q1 - 1)/(2*(q2 + 1))
    var_Chao2=(t-1)/t*q1*(q1 - 1)/2 + ((t - 1)/t)^2*q1*(2*q1-1)^2/4-((t-1)/t)^2*q1^4/4/s_Chao2
  } else {
    s_Chao2 <- D
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)Q(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*Q(i, x))))^2/t
    var_Chao2 <- var_obs
  }
  
  if (q4 != 0){
    s_iChao2 <- s_Chao2 + (t - 3)/t*q3/4/q4*max(q1 - (t - 3)/(t - 1)*q2*q3/2/q4, 0)
  } else {
    s_iChao2 <- s_Chao2 + (t - 3)/t*q3/4/(q4 + 1)*max(q1 - (t - 3)/(t - 1)*q2*q3/2, 0)  
  }
  
  diff <- function(q, x){ # fq
    q1 <- Q(1, x); q2 <- Q(2, x); q3 <- Q(3, x); q4 <- Q(4, x)
    if (q1 > 0 & q2 != 0){
      if (q == 1){
        d <- (t - 1)/t*q1/q2 - (t - 3)/t*q3/4/q4
      } else if (q == 2){
        d <- (t - 1)/t*q1^2/2/q2^2 - (t - 3)^2/t/(t - 1)*q3^2/8/q4^2
      } else if (q == 3){
        d <- (t - 3)/t*q1/4/q4
      } else {
        d <- -(t - 3)/t*q1*q3/4/q4^2 + (t - 3)^2/t/(t - 1)*q2*q3^2/4/q4^3
      }
    } else if (q1 > 1 & q2 == 0){
      if (q == 1){
        d <- (t - 1)/t*(2*q1 - 1)/2/(q2 + 1) + (t - 3)/t*q3/4/q4
      } else if (q == 2){
        d <- -(t - 1)/t*q1*(q1 - 1)/2/(q2 + 1)^2
      } else if (q == 3){
        d <- (t - 3)/t*q1/4/q4
      } else {
        d <- -(t - 3)/t*q1*q3/4/q4^2
      }    
    } else {
      d=0
    }
    return(d)
  }
  COV.q <- function(i,j){
    if (i == j){
      cov.q <- Q(i, x)*(1 - Q(i, x)/s_iChao2)
    } else {
      cov.q <- -Q(i, x)*Q(j, x)/s_iChao2
    }     
    return(cov.q)
  }
  
  ind <- 1:4
  i <- rep(sort(unique(ind)),each = length(unique(ind)))
  j <- rep(sort(unique(ind)),length(unique(ind)))       # all combination
  
  #    if (q1 - q2*q3/2/q4 > 0 & q3 != 0){
  if (q1 - (t - 3)/(t - 1)*q2*q3/2/q4 > 0 | 
        q1 - (t - 3)/(t - 1)*q2*q3/2 > 0){
    var_iChao2 <- sum(mapply(function(i, j)diff(i, x)*diff(j, x)*COV.q(i, j), i, j))
  } else {
    var_iChao2 <- var_Chao2
  }
  
  if (var_iChao2 > 0){
    var_iChao2 <- var_iChao2
  } else {
    var_iChao2 <- NA
  }
  
  m <- round(s_iChao2 - D, 5)
  if (is.nan(m) == F){
    if (m != 0){
      C <- exp(z*sqrt(log(1 + var_iChao2/(s_iChao2 - D)^2)))
      CI_iChao2 <- c(D + (s_iChao2 - D)/C, D + (s_iChao2 - D)*C)
    } else {
      i <- c(1:max(x))
      i <- i[unique(x)]
      var_obs <- sum(sapply(i, function(i)Q(i, x)*(exp(-i) - exp(-2*i)))) - 
        (sum(sapply(i, function(i)i*exp(-i)*Q(i, x))))^2/t
      var_iChao2 <- var_obs
      P <- sum(sapply(i, function(i)Q(i, x)*exp(-i)/D))
      CI_iChao2 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
    }
  }else{
    CI_iChao2 <- c(NaN, NaN)
  }
  
  table <- matrix(c(s_iChao2, sqrt(var_iChao2), CI_iChao2), ncol = 4)
  colnames(table) <- c("Estimate", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
  rownames(table) <- "improved Chao2 (Chao, 1987)"
  return(table)
  
}
