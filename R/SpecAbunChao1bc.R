SpecAbunChao1bc <-
function(data, k, conf){
  data <- as.numeric(data)
  f <- function(i, data){length(data[which(data == i)])}
  basic <- function(data, k){
    
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
    
    BASIC.DATA <- matrix(paste(c("n", "D", "k", "n_rare", "D_rare", "C_rare", "CV_rare", "CV1_rare", "n_abun", "D_abun"),
                               round(c(n, D, k, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun), 1),
                               sep = "="), ncol = 1)
    colnames(BASIC.DATA) <- c("Value")
    rownames(BASIC.DATA) <- c("Number of observed individuals", "Number of observed species","Cut-off point",
                              "Number of observed in dividuals for rare species", "Number of observed species for rare species",
                              "Estimation of the sample converage for rare species",
                              "Estimation of CV for rare species in ACE", "Estimation of CV1 for rare species in ACE-1",
                              "Number of observed species for abundant species", "Number of observed species for abundant species")
    return(list(BASIC.DATA, n, D, n_rare, D_rare, C_rare, CV_rare, CV1_rare, n_abun, D_abun))
  }
  
  z <- -qnorm((1 - conf)/2)
 
  n <- basic(data, k)[[2]]
  D <- basic(data, k)[[3]]
  n_rare <- basic(data, k)[[4]]
  D_rare <- basic(data, k)[[5]]
  C_rare <- basic(data, k)[[6]] 
  CV_rare <- basic(data, k)[[7]]
  CV1_rare <- basic(data, k)[[8]]
  n_abun <- basic(data, k)[[9]]
  D_abun <- basic(data, k)[[10]]
  x <- data[which(data != 0)]
  #############################
  S_Chao1_bc <- D + (n - 1)/n*f(1, x)*(f(1, x) - 1)/(2*(f(2, x) + 1))
  
  if (f(2, x) > 0){
    var_Chao1_bc <- (n - 1)/n*f(1, x)*(f(1, x) - 1)/2/(f(2, x) + 1) + 
      ((n - 1)/n)^2*f(1, x)*(2*f(1, x) - 1)^2/4/(f(2, x) + 1)^2 + ((n - 1)/n)^2*f(1, x)^2*f(2, x)*(f(1, x) - 1)^2/4/(f(2, x) + 1)^4
  }else{
    var_Chao1_bc <- (n - 1)/n*f(1, x)*(f(1, x) - 1)/2 + 
      ((n - 1)/n)^2*f(1, x)*(2*f(1, x) - 1)^2/4 - ((n - 1)/n)^2*f(1, x)^4/4/S_Chao1_bc
  }
  t <- round(S_Chao1_bc - D, 5)
  if (t != 0){
    K <- exp(z*sqrt(log(1 + var_Chao1_bc/t^2)))
    CI_Chao1_bc <- c(D + t/K, D + t*K)
  } else {
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
    var_Chao1_bc <- var_obs
    P <- sum(sapply(i, function(i)f(i, x)*exp(-i)/D))
    CI_Chao1_bc <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  table <- matrix(c(S_Chao1_bc, sqrt(var_Chao1_bc), CI_Chao1_bc), ncol = 4)
  colnames(table) <- c("Estimate", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
  rownames(table) <- "Chao1-bc"
  return(table)
}
