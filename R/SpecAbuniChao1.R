SpecAbuniChao1 <-
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
    f1 <- f(1, x); f2 <- f(2, x); f3 <- f(3, x); f4 <- f(4, x)
    if (f1 > 0 & f2 > 0){
      s_Chao1 <- D + (n - 1)/n*f1^2/(2*f2)
      var_Chao1 <- f2*((n - 1)/n*(f1/f2)^2/2 + 
                              ((n - 1)/n)^2*(f1/f2)^3 + ((n - 1 )/n)^2*(f1/f2)^4/4)
    } else if (f1 > 1 & f2 == 0){
      s_Chao1 <- D + (n - 1)/n*f1*(f1 - 1)/(2*(f2 + 1))
      var_Chao1 <- (n - 1)/n*f1*(f1 - 1)/2 + 
        ((n - 1)/n)^2*f1*(2*f1 - 1)^2/4 - ((n - 1)/n)^2*f1^4/4/s_Chao1
    } else {
      s_Chao1 <- D
      i <- c(1:max(x))
      i <- i[unique(x)]
      var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
        (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
      var_Chao1 <- var_obs
    }
    
    if (f4 != 0){
      s_iChao1 <- s_Chao1 + f3/4/f4*max(f1 - f2*f3/2/f4, 0)
    } else {
      s_iChao1 <- s_Chao1 + f3/4/(f4 + 1)*max(f1 - f2*f3/2/(f4 + 1), 0)  
    }
    
    diff <- function(q, x){ # fq
      f1 <- f(1, x); f2 <- f(2, x); f3 <- f(3, x); f4 <- f(4, x) 
      if (f1 > 0 & f2 != 0){
        if (q == 1){
          d <- (n - 1)/n*f1/f2 - f3/4/f4
        } else if (q == 2){
          d <- (n - 1)/n*f1^2/2/f2^2 - f3^2/8/f4^2
        } else if (q == 3){
          d <- f1/4/f4
        } else {
          d <- -f1*f3/4/f4^2 + f2*f3^2/4/f4^3
        }
      } else if (f1 > 1 & f2 == 0){
        if (q == 1){
          d <- (n - 1)/n*(2*f1 - 1)/2/(f2 + 1) + f3/4/f4
        } else if (q == 2){
          d <- -(n - 1)/n*f1*(f1 - 1)/2/(f2 + 1)^2
        } else if (q == 3){
          d <- f1/4/f4
        } else {
          d <- -f1*f3/4/f4^2
        }    
      } else {
        d=0
      }
      return(d)
    }
    COV.f <- function(i,j){
      if (i == j){
        cov.f <- f(i, x)*(1 - f(i, x)/s_iChao1)
      } else {
        cov.f <- -f(i, x)*f(j, x)/s_iChao1
      }     
      return(cov.f)
    }
    
    xx <- 1:4
    i <- rep(sort(unique(xx)),each = length(unique(xx)))
    j <- rep(sort(unique(xx)),length(unique(xx)))       # all combination
    
    if (f1 - f2*f3/2/f4 > 0 & f3 != 0){
      var_iChao1 <- sum(mapply(function(i, j)diff(i, x)*diff(j, x)*COV.f(i, j), i, j))
    } else {
      var_iChao1 <- var_Chao1
    }
    
    if (var_iChao1 > 0){
      var_iChao1 <- var_iChao1
    } else {
      var_iChao1 <- NA
    }
    
    t <- round(s_iChao1 - D, 5)
    if (is.nan(t) == F){
      if (t != 0){
        C <- exp(z*sqrt(log(1 + var_iChao1/(s_iChao1 - D)^2)))
        CI_iChao1 <- c(D + (s_iChao1 - D)/C, D + (s_iChao1 - D)*C)
      } else {
        i <- c(1:max(x))
        i <- i[unique(x)]
        var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
          (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
        var_iChao1 <- var_obs
        P <- sum(sapply(i, function(i)f(i, x)*exp(-i)/D))
        CI_iChao1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
      }
    }else{
      CI_iChao1 <- c(NaN, NaN)
    }
    
    table <- matrix(c(s_iChao1, sqrt(var_iChao1), CI_iChao1), ncol = 4)
    colnames(table) <- c("Estimate", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
    rownames(table) <- "iChao1 (Chiu et al. 2014)"
    return(table)
  }
