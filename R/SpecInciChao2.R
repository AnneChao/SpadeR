SpecInciChao2 <-
function(data, k, conf){
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
  
  if (Q(1, x)>0 & Q(2, x) > 0){
    S_Chao2 <- D + (t - 1)/t*Q(1, x)^2/(2*Q(2, x))
    var_Chao2 <- Q(2, x)*((t - 1)/t*(Q(1, x)/Q(2, x))^2/2 + ((t - 1)/t)^2*(Q(1, x)/Q(2, x))^3 + ((t - 1)/t)^2*(Q(1, x)/Q(2, x))^4/4)
    
    tt <- S_Chao2 - D
    K <- exp(z*sqrt(log(1 + var_Chao2/tt^2)))
    CI_Chao2 <- c(D + tt/K, D + tt*K)
  } else if (Q(1, x)>1 & Q(2, x) == 0){
    S_Chao2 <- D+(t-1)/t*Q(1,x)*(Q(1,x)-1)/(2*(Q(2,x)+1))
    var_Chao2=(t-1)/t*Q(1,x)*(Q(1,x)-1)/2+((t-1)/t)^2*Q(1,x)*(2*Q(1,x)-1)^2/4-((t-1)/t)^2*Q(1,x)^4/4/S_Chao2
    
    tt=S_Chao2-D
    K=exp(z*sqrt(log(1+var_Chao2/tt^2)))
    CI_Chao2=c(D+tt/K,D+tt*K)
  } else {
    S_Chao2 <- D
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)Q(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*Q(i, x))))^2/t
    var_Chao2 <- var_obs
    P <- sum(sapply(i, function(i)Q(i, x)*exp(-i)/D))
    CI_Chao2<- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  table <- matrix(c(S_Chao2, sqrt(var_Chao2), CI_Chao2), ncol = 4)
  colnames(table) <- c("Estimate", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
  rownames(table) <- "Chao2 (Chao, 1987)"
  return(table)
  
}
