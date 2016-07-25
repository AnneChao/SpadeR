Chat.Sam <- function(x, t)
{
  nT <- x[1]
  y <- x[-1]
  y <- y[y>0]
  U <- sum(y)
  Q1 <- sum(y == 1)
  Q2 <- sum(y == 2)
  Q0.hat <- ifelse(Q2 == 0, (nT - 1) / nT * Q1 * (Q1 - 1) / 2, (nT - 1) / nT * Q1 ^ 2/ 2 / Q2)  #estimation of unseen species via Chao2
  A <- ifelse(Q1>0, nT*Q0.hat/(nT*Q0.hat+Q1), 1)
  Sub <- function(t){
    if(t < nT) {
      yy <- y[(nT-y)>=t]
      out <- 1 - sum(yy / U * exp(lgamma(nT-yy+1)-lgamma(nT-yy-t+1)-lgamma(nT)+lgamma(nT-t)))
    }
    #if(t < nT) out <- 1 - sum(y / U * exp(lchoose(nT - y, t) - lchoose(nT - 1, t)))
    if(t == nT) out <- 1 - Q1 / U * A
    if(t > nT) out <- 1 - Q1 / U * A^(t - nT + 1)
    out
  }
  sapply(t, Sub)
}

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

  U<-sum(x)
  C<-Chat.Sam(data, t)
  CV_squre<-max( D/C*t/(t-1)*sum(x*(x-1))/U^2-1, 0)
  CV<-CV_squre^0.5
  U_infreq<-sum(x[x<=k])
  gamma_infreq_square_1 <- max(gamma_infreq_square*(1 + Q(1, x)/C_infreq*t/(t - 1)*b1/b2/(b2 - 1)), 0)
  CV1_infreq <- sqrt(gamma_infreq_square_1)
#   BASIC.DATA <- matrix(paste(c("D", "t", "k", "D_infreq", "C_infreq", "CV_infreq", "D_freq"),
#                              round(c(D,t,k,D_infreq,C_infreq,CV_infreq,D_freq), 3),
#                              sep = "="), ncol = 1)
  BASIC.DATA <- matrix(round(c(D,t,U,C,CV,k,U_infreq,D_infreq,C_infreq,CV_infreq,CV1_infreq,D_freq), 3), ncol = 1)
  nickname <- c("D", "T", "U", "C", "CV", "k", "U_infreq", "D_infreq", "C_infreq", "CV_infreq","CV1_infreq", "D_freq")
  BASIC.DATA <- cbind(nickname, BASIC.DATA)
  colnames(BASIC.DATA)=c("Variable", "Value")
  rownames(BASIC.DATA)=c("    Number of observed species","    Number of sampling units", "    Total number of incidences",
                         "    Coverage estimate for entire dataset", "    CV for entire dataset", "    Cut-off point","    Total number of incidences in infrequent group",
                         "    Number of observed species for infrequent group","    Estimated sample coverage for infrequent group",
                         "    Estimated CV for infrequent group in ICE",
                         "    Estimated CV1 for infrequent group in ICE-1",
                         "    Number of observed species for frequent group")
  BASIC.DATA <- data.frame(BASIC.DATA)
  return(list(BASIC.DATA, D, t, D_infreq, C_infreq, CV_infreq,CV1_infreq, D_freq))
}
