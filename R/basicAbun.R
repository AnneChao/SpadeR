Chat.Ind <- function(x, m)
{
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if(f1>0 & f2>0)
  {
    a=(n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
  }
  if(f1>1 & f2==0)
  {
    a=(n-1)*(f1-1) / ( (n-1)*(f1-1) + 2 )
  }
  if(f1==1 & f2==0) {a=0}
  if(f1==0) {a=0}

  Sub <- function(m){
    if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m == n) out <- 1-f1/n*a
    if(m > n) out <- 1-f1/n*a^(m-n+1)
    out
  }
  sapply(m, Sub)
}
basicAbun <- function(data, k){
  f <- function(i, data){length(data[which(data == i)])}
  data <- as.numeric(data)

  x <- data[which(data != 0)]
  n <- sum(x)
  D <- length(x)
  n_rare <- sum(x[which(x <= k)])
  D_rare <- length(x[which(x <= k)])
  if (n_rare != 0){
    C_rare <- 1 - f(1, x)/n_rare
  } else {
    C_rare = 1
  }
  n_abun <- n - n_rare
  D_abun <- length(x[which(x > k)])

  j <- c(1:k)
  a1 <- sum(sapply(j, function(j)j*(j - 1)*f(j, x)))
  a2 <- sum(sapply(j, function(j)j*f(j, x)))
  if (C_rare != 0){
    gamma_rare_hat_square <- max(D_rare/C_rare*a1/a2/(a2 - 1) - 1, 0)
    gamma_rare_1_square <- max(gamma_rare_hat_square*(1 + (1 - C_rare)/C_rare*a1/(a2 - 1)), 0)
  }else{
    gamma_rare_hat_square <- 0
    gamma_rare_1_square <- 0
  }
  CV_rare <- sqrt(gamma_rare_hat_square)
  CV1_rare <- sqrt(gamma_rare_1_square)

  C<-Chat.Ind(x,n)
  CV<-CV.Ind(x)

#   BASIC.DATA <- matrix(paste(c("n", "D", "k", "n_rare", "D_rare", "C_rare", "CV_rare", "CV1_rare", "n_abun", "D_abun"),
#                              round(c(n, D, k, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun),3),
#                              sep = "="), ncol = 1)
  BASIC.DATA <- matrix(round(c(n, D, C, CV, k, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun),3),
                       ncol = 1)
  #
  nickname <- matrix(c("n", "D", "C", "CV", "k", "n_rare", "D_rare", "C_rare", "CV_rare", "CV1_rare", "n_abun", "D_abun"),
                     ncol = 1)
  BASIC.DATA <- cbind(nickname, BASIC.DATA)
  colnames(BASIC.DATA) <- c("Variable", "Value")
  rownames(BASIC.DATA) <- c("    Sample size", "    Number of observed species",
                            "    Coverage estimate for entire dataset","    CV for entire dataset","    Cut-off point",
                            "    Number of observed individuals for rare group", "    Number of observed species for rare group",
                            "    Estimate of the sample coverage for rare group",
                            "    Estimate of CV for rare group in ACE", "    Estimate of CV1 for rare group in ACE-1",
                            "    Number of observed individuals for abundant group", "    Number of observed species for abundant group")
#   rownames(BASIC.DATA) <- c("(Number of observed individuals)                       n        =",
#                             "(Number of observed species)                           D        =",
#                             "(Cut-off point)                                        k        =",
#                             "(Number of observed individuals for rare species)      n_rare   =",
#                             "(Number of observed species for rare species)          D_rare   =",
#                             "(Estimation of the sample coverage for rare species)   C_rare   =",
#                             "(Estimation of CV for rare species in ACE)             CV_rare  =",
#                             "(Estimation of CV1 for rare species in ACE-1)          CV1_rare =",
#                             "(Number of observed individuals for abundant species)  n_abun   =",
#                             "(Number of observed species for abundant species)      D_abun   =")
  BASIC.DATA <- data.frame(BASIC.DATA)
  return(list(BASIC.DATA, n, D, C, CV, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun))
}
