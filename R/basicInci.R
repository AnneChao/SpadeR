basicInci <-
function(data, k){
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
  gamma_infreq_square <- max(D_infreq/C_infreq*t/(t - 1)*b1/b2/(b2 -1) - 1, 0)
  CV_infreq <- sqrt(gamma_infreq_square)
  D_freq <- length(x[which(x > k)])

 
#   BASIC.DATA <- matrix(paste(c("D", "t", "k", "D_infreq", "C_infreq", "CV_infreq", "D_freq"),
#                              round(c(D,t,k,D_infreq,C_infreq,CV_infreq,D_freq), 3),
#                              sep = "="), ncol = 1)
  BASIC.DATA <- matrix(round(c(D,t,k,D_infreq,C_infreq,CV_infreq,D_freq), 3), ncol = 1)
  nickname <- c("D", "t", "k", "D_infreq", "C_infreq", "CV_infreq", "D_freq")
  BASIC.DATA <- cbind(nickname, BASIC.DATA)
  colnames(BASIC.DATA)=c("Variable", "Value")
  rownames(BASIC.DATA)=c("Number of observed species","Number of sample/quadrats","Cut-off point",
                         "Number of observed species for infrequent species","Estimated sample coverage for infrequent species",
                         "Estimated CV for infrequent species",
                         "Number of observed species for frequent species")
  BASIC.DATA <- data.frame(BASIC.DATA)
  return(list(BASIC.DATA, D, t, D_infreq, C_infreq, CV_infreq, D_freq))
}
