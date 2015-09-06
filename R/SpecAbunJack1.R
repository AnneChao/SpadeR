SpecAbunJack1 <-
function(data, k, conf){
  data <- as.numeric(data)
  f <- function(i, data){length(data[which(data == i)])}
  basicAbun <- function(data, k){
    
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
  
  n <- basicAbun(data, k)[[2]]
  D <- basicAbun(data, k)[[3]]
  n_rare <- basicAbun(data, k)[[4]]
  D_rare <- basicAbun(data, k)[[5]]
  C_rare <- basicAbun(data, k)[[6]] 
  CV_rare <- basicAbun(data, k)[[7]]
  CV1_rare <- basicAbun(data, k)[[8]]
  n_abun <- basicAbun(data, k)[[9]]
  D_abun <- basicAbun(data, k)[[10]]
  x <- data[which(data != 0)]
  #############################
  S_1st_JK <- function(x, k){
    S_1st_jk <- D + (n - 1)/n*f(1, x)
  }
  s_1st_jk <- S_1st_JK(x, k)
  #### differential ####
  diff <- function(q){
    if ( q == 1){
      d <- 1 + (n - 1)/n
    } else {
      d <- 1
    }
    return(d)
  }
  
  COV.f <- function(i,j){
    if (i == j){
      cov.f <- f(i, x)*(1 - f(i, x)/s_1st_jk)
    } else {
      cov.f <- -f(i, x)*f(j, x)/s_1st_jk
    }     
    return(cov.f)
  }
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_1st <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.f(i, j), i, j))
  if (var_1st > 0){
    var_1st <- var_1st
  } else {
    var_1st <- NA
    cat("Warning: In this case, it can't estimate the variance of 1st-order-jackknife estimation", "\n\n")
  }
  ######################
  t <- round(s_1st_jk - D, 5)
  if (t != 0){
    C <- exp(z*sqrt(log(1 + var_1st/(s_1st_jk - D)^2)))
    CI_1st_jk <- c(D + (s_1st_jk - D)/C, D + (s_1st_jk - D)*C)
  } else {
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
    var_1st <- var_obs      
    P <- sum(sapply(i, function(i)f(i, x)*exp(-i)/D))
    CI_1st_jk <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  table <- matrix(c(s_1st_jk, sqrt(var_1st), CI_1st_jk), ncol = 4)
  colnames(table) <- c("Estimator", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
  rownames(table) <- "1st order jackknife"
  return(table)
}
