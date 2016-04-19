Jaccard_Sorensen_Abundance_equ=function(datatype = c("abundance", "incidence"),X1,X2,boot)
{
     if(datatype == "incidence")
     {
        shared.species.hat=PanEstFun.Sam(X1,X2)
        alpha.species.hat=SpecInciChao2(X1, k=10, conf=0.95)[1,1]+SpecInciChao2(X2, k=10, conf=0.95)[1,1]
        Esti.Jaccard=shared.species.hat/(alpha.species.hat-shared.species.hat)
        Esti.Sorensen=2*shared.species.hat/alpha.species.hat
        w=X1[1];z=X2[1]
        X1=X1[-1];X2=X2[-1]
     }
     n1=sum(X1);n2=sum(X2)
     if(datatype == "abundance"){w=n1;z=n2}
     I=which(X1>0 & X2>0)

     MLE.Jaccard=length(I)/sum(X1+X2>0)
     MLE.Sorensen=2*length(I)/(sum(X1>0)+sum(X2>0))
     #############################################
     if(datatype == "abundance")
     {
        shared.species.hat=PanEstFun(X1,X2)
        alpha.species.hat=SpecAbunChao1(X1, k=10, conf=0.95)[1,1]+SpecAbunChao1(X2, k=10, conf=0.95)[1,1]
        Esti.Jaccard=shared.species.hat/(alpha.species.hat-shared.species.hat)
        Esti.Sorensen=2*shared.species.hat/alpha.species.hat
     }
     #############################################
     MLE.Lennon=length(I)/(min(sum(X1>0),sum(X2>0)))
     MLE.Bray_Curtis=sum(sapply(I,function(I) 2*min(X1[I],X2[I])))/(n1+n2) 
     Morisita_Horn=2*sum(X1[I]/n1*X2[I]/n2)/sum((X1/n1)^2+(X2/n2)^2)
     Morisita_Original=2*sum(X1[I]/n1*X2[I]/n2)/sum( X1*(X1-1)/n1/(n1-1) +X2*(X2-1)/n2/(n2-1))  
 
     U_tilde=sum(X1[I])/n1;V_tilde=sum(X2[I])/n2
     fplus1=sum(X1>0 & X2==1);fplus2=sum(X1>0 & X2==2)
     fplus2=ifelse(fplus2>0,fplus2,1)
     f1plus=sum(X1==1 & X2>0);f2plus=sum(X1==2 & X2>0)
     f2plus=ifelse(f2plus>0,f2plus,1)
     U_hat=U_tilde+(z-1)/z*fplus1/2/fplus2*sum(X1[I][X2[I]==1])/n1
     U_hat=ifelse(U_hat<=1,U_hat,1)
     V_hat=V_tilde+(w-1)/w*f1plus/2/f2plus*sum(X2[I][X1[I]==1])/n2
     V_hat=ifelse(V_hat<=1,V_hat,1)
     JAu=U_tilde*V_tilde/(U_tilde+V_tilde-U_tilde*V_tilde)
     JAa=U_hat*V_hat/(U_hat+V_hat-U_hat*V_hat)
     SAu=2*U_tilde*V_tilde/(U_tilde+V_tilde)
     SAa=2*U_hat*V_hat/(U_hat+V_hat)
    
     
     p1=X1/w;p2=X2/z
     boot.Jaccard=rep(0,boot)
     boot.Esti.Jaccard=rep(0,boot)
     boot.Sorensen=rep(0,boot)
     boot.Esti.Sorensen=rep(0,boot)
     boot.Lennon=rep(0,boot)
     boot.Bray_Curtis=rep(0,boot)
     boot.Morisita_Horn=rep(0,boot)
     boot.Morisita_Original=rep(0,boot)
 
     
     boot.U_hat=rep(0,boot)
     boot.V_hat=rep(0,boot)
       boot.JAu=rep(0,boot)
       boot.JAa=rep(0,boot)
       boot.SAu=rep(0,boot)
       boot.SAa=rep(0,boot)
     for(h in 1:boot)
     {
        if(datatype == "abundance")
        {
           boot.X1=rmultinom(1,w,p1)
           boot.X2=rmultinom(1,z,p2)
           boot.n1=sum( boot.X1);boot.n2=sum(boot.X2)
           boot.shared.species.hat=PanEstFun(boot.X1,boot.X2)
           boot.alpha.species.hat=SpecAbunChao1(boot.X1, k=10, conf=0.95)[1,1]+SpecAbunChao1(boot.X2, k=10, conf=0.95)[1,1]
           boot.Esti.Jaccard[h]=boot.shared.species.hat/(boot.alpha.species.hat-boot.shared.species.hat)
           boot.Esti.Sorensen[h]=2*boot.shared.species.hat/boot.alpha.species.hat
        }
        if(datatype == "incidence")
        {
           boot.X1=sapply(1:length(p1), function(i)  rbinom(1, w, p1[i]))
           boot.X2=sapply(1:length(p2), function(i)  rbinom(1, z, p2[i]))
           boot.shared.species.hat=PanEstFun.Sam(c(w,boot.X1), c(z,boot.X2) )
           boot.alpha.species.hat=SpecInciChao2(c(w,boot.X1), k=10, conf=0.95)[1,1]+SpecInciChao2(c(z,boot.X2), k=10, conf=0.95)[1,1]
           boot.Esti.Jaccard[h]=boot.shared.species.hat/(boot.alpha.species.hat-boot.shared.species.hat)
           boot.Esti.Sorensen[h]=2*boot.shared.species.hat/boot.alpha.species.hat
           n1=sum(boot.X1)
           n2=sum(boot.X2)
        }
        I=which( boot.X1>0 &  boot.X2>0)
        
        boot.Jaccard[h]= length(I)/sum(boot.X1+boot.X2>0)
        boot.Sorensen[h]=2*length(I)/(sum(boot.X1>0)+sum(boot.X2>0))
        boot.Lennon[h]=length(I)/(min(sum(boot.X1>0),sum(boot.X2>0)))
        boot.Bray_Curtis[h]=sum(sapply(I,function(I) 2*min(boot.X1[I],boot.X2[I])))/sum(boot.X1+boot.X2)
        boot.Morisita_Horn[h]=2*sum(boot.X1[I]/n1*boot.X2[I]/n2)/sum((boot.X1/n1)^2+(boot.X2/n2)^2)
        boot.Morisita_Original[h]=2*sum(boot.X1[I]/n1*boot.X2[I]/n2)/sum( boot.X1*(boot.X1-1)/n1/(n1-1) +boot.X2*(boot.X2-1)/n2/(n2-1))
  
        boot.U_tilde=sum( boot.X1[I])/n1
        boot.V_tilde=sum( boot.X2[I])/n2
        fplus1=sum(boot.X1>0 & boot.X2==1);fplus2=sum(boot.X1>0 & boot.X2==2)
        fplus2=ifelse(fplus2>0,fplus2,1)
        f1plus=sum(boot.X1==1 & boot.X2>0);f2plus=sum(boot.X1==2 & boot.X2>0)
        f2plus=ifelse(f2plus>0,f2plus,1)
        boot.U_hat[h]=boot.U_tilde+(z-1)/z*fplus1/2/fplus2*sum(boot.X1[I][boot.X2[I]==1])/n1
        boot.U_hat[h]=ifelse(boot.U_hat[h]<=1,boot.U_hat[h],1)
        boot.V_hat[h]=boot.V_tilde+(w-1)/w*f1plus/2/f2plus*sum(boot.X2[I][boot.X1[I]==1])/n2
        boot.V_hat[h]=ifelse(boot.V_hat[h]<=1,boot.V_hat[h],1)
        boot.JAu[h]=boot.U_tilde*boot.V_tilde/(boot.U_tilde+boot.V_tilde-boot.U_tilde*boot.V_tilde)
        boot.JAa[h]=boot.U_hat[h]*boot.V_hat[h]/(boot.U_hat[h]+boot.V_hat[h]-boot.U_hat[h]*boot.V_hat[h])
        boot.SAu[h]=2*boot.U_tilde*boot.V_tilde/(boot.U_tilde+boot.V_tilde)
        boot.SAa[h]=2*boot.U_hat[h]*boot.V_hat[h]/(boot.U_hat[h]+boot.V_hat[h])
     }
     a=matrix(0,12,6)
     a[1,]=c(min(MLE.Jaccard,1),sd(boot.Jaccard),rep(0,4))
     a[2,]=c(min(Esti.Jaccard,1),sd(boot.Esti.Jaccard),rep(0,4))
     a[3,]=c(min(MLE.Sorensen,1),sd(boot.Sorensen),rep(0,4))
     a[4,]=c(min(Esti.Sorensen,1),sd(boot.Esti.Sorensen),rep(0,4))
     a[5,]=c(MLE.Lennon,sd(boot.Lennon),rep(0,4))
     a[6,]=c(min(MLE.Bray_Curtis,1),sd(boot.Bray_Curtis),rep(0,4))
     a[7,]=c(min(Morisita_Horn,1),sd(boot.Morisita_Horn),rep(0,4))
     a[8,]=c(min(Morisita_Original,1),sd(boot.Morisita_Original),rep(0,4))
     a[9,]=c(min(JAu,1),sd(boot.JAu),U_tilde,V_tilde,rep(0,2))
     a[10,]=c(min(JAa,1),sd(boot.JAa),U_hat,sd(boot.U_hat),V_hat,sd(boot.V_hat))
     a[11,]=c(min(SAu,1),sd(boot.SAu),U_tilde,V_tilde,rep(0,2))
     a[12,]=c(min(SAa,1),sd(boot.SAa),U_hat,sd(boot.U_hat),V_hat,sd(boot.V_hat))
     round(a,4)
}

Two_horn_MLE_equ=function(X1,X2)
{
  
   horn_MLE_equ=function(X1,X2)
   {
       n1=sum(X1)
       n2=sum(X2)
       w1=n1/(n1+n2);w2=1-w1
       pool.X= X1+X2
       pool.n= n1+n2
       pool.phat=pool.X/pool.n;pool.phat=pool.phat[pool.phat>0]
       p1hat=X1/n1;p1hat=p1hat[p1hat>0]
       p2hat=X2/n2;p2hat=p2hat[p2hat>0]
       Hr=-sum(pool.phat*log(pool.phat))
       Ha=-w1*sum(p1hat*log(p1hat))-w2*sum(p2hat*log(p2hat))
       horn=(Hr-Ha)/(-w1*log(w1)-w2*log(w2))
       horn
   }
   n1=sum(X1)
   n2=sum(X2)
   horn=horn_MLE_equ(X1,X2)
   boot.horn=rep(0,50)
   boot.p1=X1/n1
   boot.p2=X2/n2
   for(h in 1:50)
   {
     boot.X1=rmultinom(1,n1,boot.p1)
     boot.X2=rmultinom(1,n2,boot.p2)
     boot.horn[h]=horn_MLE_equ(boot.X1,boot.X2)  
   }
   out=c(min(horn,1),sd(boot.horn));out=round(out,4)
   return(out)
}



Chao1_equ=function(x,conf=0.95)
{
  z <--qnorm((1 - conf)/2)
  x=x[x>0]
  D=sum(x>0)
  f1=sum(x==1)
  f2=sum(x==2)
  n=sum(x)
  if (f1 > 0 & f2 > 0)
  {
    S_Chao1 <- D + (n - 1)/n*f1^2/(2*f2)
    var_Chao1 <- f2*((n - 1)/n*(f1/f2)^2/2 + 
                            ((n - 1)/n)^2*(f1/f2)^3 + ((n - 1 )/n)^2*(f1/f2)^4/4)
    
    t <- S_Chao1 - D
    K <- exp(z*sqrt(log(1 + var_Chao1/t^2)))
    CI_Chao1 <- c(D + t/K, D + t*K)
  } 
  else if (f1 > 1 & f2 == 0)
  {
    S_Chao1 <- D + (n - 1)/n*f1*(f1 - 1)/(2*(f2 + 1))
    var_Chao1 <- (n - 1)/n*f1*(f1 - 1)/2 + 
      ((n - 1)/n)^2*f1*(2*f1 - 1)^2/4 - ((n - 1)/n)^2*f1^4/4/S_Chao1
    
    t <- S_Chao1 - D
    K <- exp(z*sqrt(log(1 + var_Chao1/t^2)))
    CI_Chao1 <- c(D + t/K, D + t*K)
  } 
  else 
  {
    S_Chao1 <- D
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i) sum(x==i)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*sum(x==i))))^2/n
    var_Chao1 <- var_obs
    P <- sum(sapply(i, function(i) sum(x==i)*exp(-i)/D))
    CI_Chao1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  return( c( round(c(S_Chao1,var_Chao1^0.5,CI_Chao1[1],CI_Chao1[2]),1),conf)    )
}
Chao1_bc_equ=function(x,conf=0.95)
{
  z <- -qnorm((1 - conf)/2)
  x=x[x>0]
  D=sum(x>0)
  f1=sum(x==1)
  f2=sum(x==2)
  n=sum(x)
  
  S_Chao1_bc <- D + (n - 1)/n*f1*(f1 - 1)/(2*(f2 + 1))
  var_Chao1_bc <- (n - 1)/n*f1*(f1 - 1)/2/(f2 + 1) + 
    ((n - 1)/n)^2*f1*(2*f1 - 1)^2/4/(f2 + 1)^2 + ((n - 1)/n)^2*f1^2*f2*(f1 - 1)^2/4/(f2 + 1)^4
  
  t <- round(S_Chao1_bc - D, 5)
  if (t != 0)
  {
    K <- exp(z*sqrt(log(1 + var_Chao1_bc/t^2)))
    CI_Chao1_bc <- c(D + t/K, D + t*K)
  } 
  if(t == 0)
  {
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)sum(x==i)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*sum(x==i))))^2/n
    var_Chao1_bc <- var_obs
    P <- sum(sapply(i, function(i)sum(x==i)*exp(-i)/D))
    CI_Chao1_bc <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  round(c(S_Chao1_bc,var_Chao1_bc^0.5,CI_Chao1_bc[1],CI_Chao1_bc[2]) ,1)  
}
SpecAbunAce_equ<- function(data, k=10, conf=0.95)
{
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
  S_ACE <- function(x, k){
    j <- c(1:k)
    a1 <- sum(sapply(j, function(j)j*(j - 1)*f(j, x)))
    a2 <- sum(sapply(j, function(j)j*f(j, x)))
    if (C_rare != 0){
      gamma_rare_hat_square <- max(D_rare/C_rare*a1/a2/(a2 - 1) - 1, 0)
    }else{
      gamma_rare_hat_square <- 0
    }
    S_ace <- D_abun + D_rare/C_rare + f(1, x)/C_rare*gamma_rare_hat_square
    return(list(S_ace, gamma_rare_hat_square))
  }
  s_ace <- S_ACE(x, k)[[1]]
  gamma_rare_hat_square <- S_ACE(x, k)[[2]]
  #### differential ####
  u <- c(1:k)    
  diff <- function(q){
    if (gamma_rare_hat_square != 0){
      si <- sum(sapply(u, function(u)u*(u - 1)*f(u, x)))
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(D_rare*si + f(1, x)*si) - 
             f(1, x)*D_rare*si*(-2*(1 - f(1, x)/n_rare)*(n_rare - f(1, x))/n_rare^2*n_rare*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*(2*n_rare - 1))
          )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2 - #g2
          (1 - f(1, x)/n_rare + f(1, x)*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g3
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)*(si + D_rare*q*(q - 1)) - 
             f(1, x)*D_rare*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                                  (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
          )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2 + #g2
          (q*(f(1, x))^2/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g3
      }
      return(d)
    } else {
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1 
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1
      }
      return(d)  
    }
  }
  COV.f <- function(i,j){
    if (i == j){
      cov.f <- f(i, x)*(1 - f(i, x)/s_ace)
    } else {
      cov.f <- -f(i, x)*f(j, x)/s_ace
    }     
    return(cov.f)
  }
  
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ace <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.f(i, j), i, j))
  if (var_ace > 0){
    var_ace <- var_ace
  } else {
    var_ace <- NA
  }
  ######################
  t <- round(s_ace - D, 5)
  if (is.nan(t) == F){
    if (t != 0){
      C <- exp(z*sqrt(log(1 + var_ace/(s_ace - D)^2)))
      CI_ACE <- c(D + (s_ace - D)/C, D + (s_ace - D)*C)
    } else {
      i <- c(1:max(x))
      i <- i[unique(x)]
      var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
        (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
      var_ace <- var_obs
      P <- sum(sapply(i, function(i)f(i, x)*exp(-i)/D))
      CI_ACE <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
    }
  }else{
    CI_ACE <- c(NaN, NaN)
  }
  
  table <- matrix(c(s_ace, sqrt(var_ace), CI_ACE), ncol = 4)
  #colnames(table) <- c("Estimate", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
  #rownames(table) <- "ACE (Chao & Lee, 1992)" 
  return(round(table,1))
}
SpecAbunAce1_equ<- function(data ,k=10, conf=0.95)
{
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
  S_ACE1 <- function(x, k){
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
    s_ace1 <- D_abun + D_rare/C_rare + f(1, x)/C_rare*gamma_rare_1_square
    
    return(list(s_ace1, gamma_rare_1_square))
  }
  s_ace1 <- S_ACE1(x, k)[[1]]
  gamma_rare_1_square <- S_ACE1(x, k)[[2]]
  #### differential ####
  u <- c(1:k)    
  diff <- function(q){
    if (gamma_rare_1_square != 0){
      u <- c(1:k)
      si <- sum(sapply(u, function(u)u*(u-1)*f(u, x)))
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(D_rare*si + f(1, x)*si) - 
             f(1, x)*D_rare*si*(-2*(1 - f(1, x)/n_rare)*(n_rare - f(1, x))/n_rare^2*n_rare*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*(2*n_rare - 1))
          )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2 - #g2
          (1 - f(1, x)/n_rare + f(1, x)*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g3
          ((1 - f(1, x)/n_rare)^3*(n_rare*(n_rare - 1))^2*(2*f(1, x)*D_rare*si^2 + f(1, x)^2*si^2) - #g4
             f(1, x)^2*D_rare*si^2*(3*(1 - f(1, x)/n_rare)^2*(f(1, x) - n_rare)/(n_rare)^2*(n_rare*(n_rare - 1))^2 + 
                                      (1 - f(1, x)/n_rare)^3*2*n_rare*(n_rare - 1)^2 + (1 - f(1, x)/n_rare)^3*n_rare^2*2*(n_rare - 1)) 
          )/(1 - f(1, x)/n_rare)^6/n_rare^4/(n_rare - 1)^4 - 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*(2*f(1, x)*si) - #g5
             f(1, x)^2*si*(2*(1 - f(1, x)/n_rare)*(f(1, x) - n_rare)/n_rare^2*n_rare*(n_rare - 1) + 
                             (1 - f(1, x)/n_rare)^2*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare) 
          )/(1 - f(1, x)/n_rare)^4/n_rare^2/(n_rare - 1)^2
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g1
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)*(si + D_rare*q*(q - 1)) - 
             f(1, x)*D_rare*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                                  (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
          )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2 + #g2
          (q*(f(1, x))^2/n_rare^2)/(1 - f(1, x)/n_rare)^2 + #g3
          ((1 - f(1, x)/n_rare)^3*n_rare^2*(n_rare - 1)^2*f(1, x)^2*(si^2 + 2*D_rare*si*q*(q - 1)) - #g4
             f(1, x)^2*D_rare*si^2*(3*(1 - f(1, x)/n_rare)^2*(f(1, x)*q/n_rare^2)*(n_rare*(n_rare - 1))^2 + 
                                      2*(1 - f(1, x)/n_rare)^3*n_rare*q*(n_rare - 1)^2 + 2*(1 - f(1, x)/n_rare)^3*n_rare^2*(n_rare - 1)*q)   
          )/(1 - f(1, x)/n_rare)^6/(n_rare)^4/(n_rare - 1)^4 - 
          ((1 - f(1, x)/n_rare)^2*n_rare*(n_rare - 1)*f(1, x)^2*q*(q - 1) - #g5
             f(1, x)^2*si*(2*(1 - f(1, x)/n_rare)*f(1, x)*q/n_rare^2*n_rare*(n_rare - 1) + 
                             (1 - f(1, x)/n_rare)^2*q*(n_rare - 1) + (1 - f(1, x)/n_rare)^2*n_rare*q)
          )/(1 - f(1, x)/n_rare)^4/(n_rare)^2/(n_rare - 1)^2
      }
      return(d)
    } else {
      u <- c(1:k)
      si <- sum(sapply(u, function(u)u*(u-1)*f(u, x)))
      if ( q == 1){
        d <- (1 - f(1, x)/n_rare + D_rare*(n_rare - f(1, x))/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1 
      } else if(q > k){
        d <- 1
      } else {
        d <- (1 - f(1, x)/n_rare - D_rare*q*f(1, x)/n_rare^2)/(1 - f(1, x)/n_rare)^2 #g1
      }
      return(d)
    }
  }
  
  COV.f <- function(i,j){
    if (i == j){
      cov.f <- f(i, x)*(1 - f(i, x)/s_ace1)
    } else {
      cov.f <- -f(i, x)*f(j, x)/s_ace1
    }     
    return(cov.f)
  }
  
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ace1 <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.f(i, j), i, j))
  if (var_ace1 > 0){
    var_ace1 <- var_ace1
  } else {
    var_ace1 <- NA
  }
  ######################
  t <- round(s_ace1 - D, 5)
  if (is.nan(t) == F){
    if (t != 0){
      C <- exp(z*sqrt(log(1 + var_ace1/(s_ace1 - D)^2)))
      CI_ACE1 <- c(D + (s_ace1 - D)/C, D + (s_ace1 - D)*C)
    } else {
      i <- c(1:max(x))
      i <- i[unique(x)]
      var_obs <- sum(sapply(i, function(i)f(i, x)*(exp(-i) - exp(-2*i)))) - 
        (sum(sapply(i, function(i)i*exp(-i)*f(i, x))))^2/n
      var_ace1 <- var_obs
      P <- sum(sapply(i, function(i)f(i, x)*exp(-i)/D))
      CI_ACE1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
    }
  }else{
    CI_ACE1 <- c(NaN, NaN)
  }
  
  table <- matrix(c(s_ace1, sqrt(var_ace1), CI_ACE1), ncol = 4)
  #colnames(table) <- c("Estimate", "Est_s.e.", paste(conf*100,"% Lower Bound"), paste(conf*100,"% Upper Bound"))
  #rownames(table) <- "ACE-1 (Chao & Lee, 1992)"
  return(round(table,1))
}
SpecInciChao2 <-function(data, k=10, conf=0.95)
{
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
  return(table)
}
SpecInciChao2bc <-function(data, k=10, conf=0.95)
{
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
  
  
  S_Chao2_bc <- D + (t - 1)/t*Q(1, x)*(Q(1, x) - 1)/(2*(Q(2, x) + 1))
  var_Chao2_bc <- (t - 1)/t*Q(1, x)*(Q(1, x) - 1)/2/(Q(2, x) + 1) + ((t - 1)/t)^2*Q(1, x)*(2*Q(1, x) - 1)^2/4/(Q(2, x) + 1)^2 + ((t - 1)/t)^2*Q(1, x)^2*Q(2, x)*(Q(1, x) - 1)^2/4/(Q(2, x) + 1)^4
  
  tt <- S_Chao2_bc - D
  if (tt != 0){
    K <- exp(z*sqrt(log(1 + var_Chao2_bc/tt^2)))
    CI_Chao2_bc <- c(D + tt/K, D + tt*K)
  } else {
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)Q(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*Q(i, x))))^2/t
    var_Chao2_bc <- var_obs
    P <- sum(sapply(i, function(i)Q(i, x)*exp(-i)/D))
    CI_Chao2_bc <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  table <- matrix(c(S_Chao2_bc, sqrt(var_Chao2_bc), CI_Chao2_bc), ncol = 4)
  return(table)
}
SpecInciModelh <-function(data, k=10, conf=0.95)
{  
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
  
  S_ICE <- function(x, k){
    j <- c(1:k)
    a1 <- sum(sapply(j, function(j)j*(j - 1)*Q(j, x)))
    a2 <- sum(sapply(j, function(j)j*Q(j, x)))
    gamma_infreq_square <- max(D_infreq/C_infreq*t/(t - 1)*a1/a2/(a2 - 1) - 1,0)          
    s_ice <- D_freq + D_infreq/C_infreq + Q(1, x)/C_infreq*gamma_infreq_square
    CV_infreq_h <- sqrt(gamma_infreq_square)
    return(c(s_ice, CV_infreq_h))
  }
  s_ice <- S_ICE(x, k)[1]
  CV_infreq_h <-  S_ICE(x, k)[2]
  #### differential ####
  u <- c(1:k)    
  diff <- function(q){
    if (CV_infreq_h != 0){
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Q(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x))*2*Q(1, x)*(t - 1) - 
                           (t - 1)*Q(1, x)^2*((t - 1)*(Q(1, x) + n_infreq) + 2*Q(2, x)))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*(D_infreq*si + Q(1, x)*si) - #g3
                       Q(1, x)*D_infreq*si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*(n_infreq - 1) + C_infreq^2*n_infreq)             
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          (C_infreq - Q(1, x)*dc_infreq)/C_infreq^2 #g4
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*(2*(t - 1)*Q(1, x) + 2*(n_infreq + 2*Q(2, x))))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Q(1, x)*(si + 2*D_infreq) - Q(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*2*(n_infreq - 1) + C_infreq^2*n_infreq*2)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Q(1, x)*dc_infreq)/C_infreq^2 #g4
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*((t - 1)*Q(1, x)*q + 2*Q(2, x)*q))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Q(1, x)*(si + q*(q - 1)*D_infreq) - Q(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*q*(n_infreq - 1) + C_infreq^2*n_infreq*q)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Q(1, x)*dc_infreq)/C_infreq^2 #g4
      }
      return(d)
    }else{
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Q(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x))*2*Q(1, x)*(t - 1) - 
                           (t - 1)*Q(1, x)^2*((t - 1)*(Q(1, x) + n_infreq) + 2*Q(2, x)))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*(2*(t - 1)*Q(1, x) + 2*(n_infreq + 2*Q(2, x))))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*((t - 1)*Q(1, x)*q + 2*Q(2, x)*q))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }
      return(d)
    }
  }  
  
  COV.q <- function(i,j){
    if (i == j){
      cov.q <- Q(i, x)*(1 - Q(i, x)/s_ice)
    } else {
      cov.q <- -Q(i, x)*Q(j, x)/s_ice
    }     
    return(cov.q)
  }
  
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ice <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.q(i, j), i, j))
  if (var_ice > 0){
    var_ice <- var_ice
  } else {
    var_ice <- NA
    cat("Warning: In this case, it can't estimate the variance of Model(h) estimation", "\n\n")
  }
  ######################
  if (round(s_ice - D, 5) != 0){
    C <- exp(z*sqrt(log(1 + var_ice/(s_ice - D)^2)))
    CI_Model_h <- c(D + (s_ice - D)/C, D + (s_ice - D)*C)
  }else{
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)Q(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*Q(i, x))))^2/t
    var_ice <- var_obs
    P <- sum(sapply(i, function(i)Q(i, x)*exp(-i)/D))
    CI_Model_h <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  table <- matrix(c(s_ice, sqrt(var_ice), CI_Model_h), ncol = 4)
  return(table)
}
SpecInciModelh1 <-function(data, k=10, conf=0.95)
{
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
  
  S_Model_H1 <- function(x, k){
    j <- c(1:k)
    a1 <- sum(sapply(j, function(j)j*(j - 1)*Q(j, x)))
    a2 <- sum(sapply(j, function(j)j*Q(j, x)))
    gamma_infreq_square <- max(D_infreq/C_infreq*t/(t - 1)*a1/a2/(a2 - 1) - 1,0)      
    gamma_infreq_square_1 <- max(gamma_infreq_square*(1 + Q(1, x)/C_infreq*t/(t - 1)*a1/a2/(a2 - 1)), 0)
    s_Model_h1 <- D_freq + D_infreq/C_infreq + Q(1, x)/C_infreq*gamma_infreq_square_1
    CV_infreq_h1 <- sqrt(gamma_infreq_square_1)
    return(c(s_Model_h1, CV_infreq_h1))
  }
  s_Model_h1 <- S_Model_H1(x, k)[1]
  CV_infreq_h1 <- S_Model_H1(x, k)[2]
  #### differential ####
  u <- c(1:k)    
  diff <- function(q){
    if (CV_infreq_h1 != 0){
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Q(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x))*2*Q(1, x)*(t - 1) - 
                           (t - 1)*Q(1, x)^2*((t - 1)*(Q(1, x) + n_infreq) + 2*Q(2, x)))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*(D_infreq*si + Q(1, x)*si) - #g3
                       Q(1, x)*D_infreq*si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*(n_infreq - 1) + C_infreq^2*n_infreq)             
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          (C_infreq - Q(1, x)*dc_infreq)/C_infreq^2 + #g4
          (t/(t - 1))^2*(C_infreq^3*n_infreq^2*(n_infreq - 1)^2*(2*Q(1, x)*D_infreq*si^2 + Q(1, x)^2*si^2) - #g5
                           Q(1, x)^2*D_infreq*si^2*(3*C_infreq^2*dc_infreq*n_infreq^2*(n_infreq - 1)^2 + C_infreq^3*2*n_infreq*(n_infreq - 1)^2 + C_infreq^3*n_infreq^2*2*(n_infreq - 1))
          )/C_infreq^6/n_infreq^4/(n_infreq - 1)^4 - 
          (t/(t - 1))*si*(C_infreq^2*n_infreq*(n_infreq - 1)*2*Q(1, x) - Q(1, x)^2*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*(n_infreq - 1) + C_infreq^2*n_infreq) #g6
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*(2*(t - 1)*Q(1, x) + 2*(n_infreq + 2*Q(2, x))))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Q(1, x)*(si + 2*D_infreq) - Q(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*2*(n_infreq - 1) + C_infreq^2*n_infreq*2)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Q(1, x)*dc_infreq)/C_infreq^2 + #g4
          (t/(t - 1))^2*Q(1, x)^2*(C_infreq^3*n_infreq^2*(n_infreq - 1)^2*(si^2 + D_infreq*2*si*2) - #g5
                                     D_infreq*si^2*(3*C_infreq^2*dc_infreq*n_infreq^2*(n_infreq - 1)^2 + C_infreq^3*2*n_infreq*2*(n_infreq - 1)^2 + C_infreq^3*n_infreq^2*2*(n_infreq - 1)*2)
          )/C_infreq^6/n_infreq^4/(n_infreq - 1)^4 - 
          t/(t - 1)*Q(1, x)^2*(C_infreq^2*n_infreq*(n_infreq - 1)*2 - si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*2*(n_infreq - 1) + C_infreq^2*2*n_infreq)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*((t - 1)*Q(1, x)*q + 2*Q(2, x)*q))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 + #g2
          t/(t - 1)*(C_infreq^2*n_infreq*(n_infreq - 1)*Q(1, x)*(si + q*(q - 1)*D_infreq) - Q(1, x)*D_infreq*si*( #g3
            2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*q*(n_infreq - 1) + C_infreq^2*n_infreq*q)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2 - 
          ( - Q(1, x)*dc_infreq)/C_infreq^2 + #g4
          (t/(t - 1))^2*Q(1, x)^2*(C_infreq^3*n_infreq^2*(n_infreq - 1)^2*(si^2 + D_infreq*2*si*q*(q - 1)) - #g5
                                     D_infreq*si^2*(3*C_infreq^2*dc_infreq*n_infreq^2*(n_infreq - 1)^2 + C_infreq^3*2*n_infreq*q*(n_infreq - 1)^2 + C_infreq^3*n_infreq^2*2*(n_infreq - 1)*q)
          )/C_infreq^6/n_infreq^4/(n_infreq - 1)^4 -
          t/(t - 1)*Q(1, x)^2*(C_infreq^2*n_infreq*(n_infreq - 1)*q*(q - 1) - #g6
                                 si*(2*C_infreq*dc_infreq*n_infreq*(n_infreq - 1) + C_infreq^2*q*(n_infreq - 1) + C_infreq^2*n_infreq*q)
          )/C_infreq^4/n_infreq^2/(n_infreq - 1)^2                                 
      }
      return(d)
    }else{
      n_infreq <- sum(x[which(x <= k)])
      si <- sum(sapply(u, function(u)u*(u-1)*Q(u, x)))
      if ( q == 1){
        dc_infreq <-  - (n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x))*2*Q(1, x)*(t - 1) - 
                           (t - 1)*Q(1, x)^2*((t - 1)*(Q(1, x) + n_infreq) + 2*Q(2, x)))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      } else if (q == 2){
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*(2*(t - 1)*Q(1, x) + 2*(n_infreq + 2*Q(2, x))))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }else if(q > k){
        d <- 1
      } else {
        dc_infreq <-  - ( - (t - 1)*Q(1, x)^2*((t - 1)*Q(1, x)*q + 2*Q(2, x)*q))/(n_infreq*((t - 1)*Q(1, x) + 2*Q(2, x)))^2
        d <- (C_infreq - D_infreq*dc_infreq)/C_infreq^2 #g2
      }
      return(d)
    }
  }  
  
  COV.q <- function(i,j){
    if (i == j){
      cov.q <- Q(i, x)*(1 - Q(i, x)/s_Model_h1)
    } else {
      cov.q <- -Q(i, x)*Q(j, x)/s_Model_h1
    }     
    return(cov.q)
  }
  
  i <- rep(sort(unique(x)),each = length(unique(x)))
  j <- rep(sort(unique(x)),length(unique(x)))       # all combination
  
  var_ice1 <- sum(mapply(function(i, j)diff(i)*diff(j)*COV.q(i, j), i, j))
  if (var_ice1 > 0){
    var_ice1 <- var_ice1
  } else {
    var_ice1 <- NA
    cat("Warning: In this case, it can't estimate the variance of Model(h)-1 estimation", "\n\n")
  }
  ######################   
  if (round(s_Model_h1 - D, 5) != 0){
    C <- exp(z*sqrt(log(1 + var_ice1/(s_Model_h1 - D)^2)))
    CI_Model_h1 <- c(D + (s_Model_h1 - D)/C, D + (s_Model_h1 - D)*C)
  } else {
    i <- c(1:max(x))
    i <- i[unique(x)]
    var_obs <- sum(sapply(i, function(i)Q(i, x)*(exp(-i) - exp(-2*i)))) - 
      (sum(sapply(i, function(i)i*exp(-i)*Q(i, x))))^2/t
    var_ice1 <- var_obs
    P <- sum(sapply(i, function(i)Q(i, x)*exp(-i)/D))
    CI_Model_h1 <- c(max(D, D/(1 - P) - z*sqrt(var_obs)/(1 - P)), D/(1 - P) + z*sqrt(var_obs)/(1 - P))  
  }
  table <- matrix(c(s_Model_h1, sqrt(var_ice1), CI_Model_h1), ncol = 4)
  return(table)
}

correct_obspi<- function(X)
{
  Sobs <- sum(X > 0)   
  n <- sum(X)		  	
  f1 <- sum(X == 1) 	
  f2 <- sum(X == 2)
  if(f1>0 & f2>0)
  {
    a=(n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n
  }
  if(f1>1 & f2==0)
  {
    a=(n-1)*(f1-1) / ( (n-1)*(f1-1) + 2 )*f1/n
  } 
  if(f1==1 & f2==0) {a=0}
  if(f1==0 ) {a=0} 	
  b <- sum(X / n * (1 - X / n) ^ n)
  w <- a / b  			
  Prob.hat <- X / n * (1 - w * (1 - X / n) ^ n)	
  Prob.hat
}

p1bar_equ=function(X)
{
  n=sum(X)
  f1=sum(X==1)
  f2=sum(X==2)
  if(f1>0 & f2>0)
  {
    a=2*f2/( (n-1)*f1+2*f2  )
  }
  if(f1>1 & f2==0)
  {
    a=2/( (n-1)*(f1-1)+2      )
  }
  if(f1==1 &  f2==0){a=0}
  if(f1==0){a=0}
  return(a)
}

KH_Bray_curtis_equ=function(X1,X2,w1)
{
  w2=1-w1
  p1bar_1=p1bar_equ(X1)
  p1bar_2=p1bar_equ(X2)
  I=which(X1*X2>0)
  Y1=X1[I];Y2=X2[I]
  n1=sum(X1);n2=sum(X2)
  f.1=sum(X1>0 & X2==1)
  f.2=sum(X1>0 & X2==2);f.2=ifelse(f.2>0,f.2,1)
  f1.=sum(X1==1 & X2>0)
  f2.=sum(X1==2 & X2>0);f2.=ifelse(f2.>0,f2.,1)
  
  temp=0
  for(i in 1:length(I) )
  {
    a1=w1^2*ifelse(Y1[i]>1,Y1[i]*(Y1[i]-1)/n1/(n1-1),p1bar_1^2)-
      2*w1*w2*ifelse(Y1[i]>1,Y1[i]/n1,p1bar_1)*ifelse(Y2[i]>1,Y2[i]/n2,p1bar_2)+
      w2^2*ifelse(Y2[i]>1,Y2[i]*(Y2[i]-1)/n2/(n2-1),p1bar_2^2)
    a1=max(0,a1)
    if(a1==0)
    {
      a1=(abs(w1*ifelse(Y1[i]>1,Y1[i]/n1,p1bar_1)-w2*ifelse(Y2[i]>1,Y2[i]/n2,p1bar_2)))^2
    }
    temp=temp+a1^0.5
  }
  I=which(X1>0 & X2==1)
  Y1=X1[I];Y2=X2[I]
  if(length(I)>0)
  {
    for(i in 1:length(I) )
    {
      a1=w1^2*ifelse(Y1[i]>1,Y1[i]*(Y1[i]-1)/n1/(n1-1),p1bar_1^2)-
        2*w1*w2*ifelse(Y1[i]>1,Y1[i]/n1,p1bar_1)*p1bar_2+
        w2^2*p1bar_2^2
      a1=max(0,a1)
      if(a1==0)
      {
        a1=(abs(w1*ifelse(Y1[i]>1,Y1[i]/n1,p1bar_1)-w2*ifelse(Y2[i]>1,Y2[i]/n2,p1bar_2)))^2
        
      }
      temp=temp+a1^0.5*f.1/2/f.2 #(1-p1bar_2)/(n2*p1bar_2)
    }
  }
  I=which(X1==1 & X2>0)
  Y1=X1[I];Y2=X2[I]
  if(length(I)>0)
  {
    for(i in 1:length(I) )
    {
      a1=w2^2*ifelse(Y2[i]>1,Y2[i]*(Y2[i]-1)/n2/(n2-1),p1bar_2^2)-
        2*w1*w2*ifelse(Y2[i]>1,Y2[i]/n2,p1bar_2)*p1bar_1+
        w1^2*p1bar_1^2
      a1=max(0,a1)
      if(a1==0)
      {
        a1=(abs(w1*ifelse(Y1[i]>1,Y1[i]/n1,p1bar_1)-w2*ifelse(Y2[i]>1,Y2[i]/n2,p1bar_2)))^2
        
      }
      temp=temp+a1^0.5*f1./2/f2.#(1-p1bar_1)/(n1*p1bar_1)
    }
  }
  f11=sum(X1==1 & X2==1)
  sumtemp=f11*abs(w1*p1bar_1-w2*p1bar_2)*(1-p1bar_1)/(n1*p1bar_1)*(1-p1bar_2)/(n2*p1bar_2)
  if(p1bar_1==0 | p1bar_2==0)
  {
    sumtemp=0
  }
  temp=temp+sumtemp
  
  ##########################
  da1=X1  
  da2=X2
  I=which(da1*da2>0); sda1=da1[I];sda2=da2[I]
  
  
  U1=sum(sda1)/n1; V1=sum(sda2)/n2;
  ff1=sum(sda2==1);#ff1=ifelse(ff1==0, 1,ff1);
  ff2=sum(sda2==2);ff2=ifelse(ff2==0, 1,ff2);
  
  f1=sum(sda1==1);#f1=ifelse(f1==0, 1,f1);
  f2=sum(sda1==2);f2=ifelse(f2==0, 1,f2);  
  
  U2=(ff1/(2*ff2))*sum(sda1[sda2==1])/n1;uC1=f1/sum(sda1);
  V2=(f1/(2*f2))*sum(sda2[sda1==1])/n2;uC2=ff1/sum(sda2);
  U=U2+U1;U=min(U,1)
  V=V1+V2;V=min(V,1)
  ##########################
  mle=sum(abs(w1*X1/n1-w2*X2/n2))
  out=w1*(1-U)+w2*(1-V)+temp
  if(out<0 | out>1)
  {
    out=mle
  } 
  return(out)
}

MLE_Braycurtis_equ=function(X1,X2,w1)
{
   w2=1-w1
   n1=sum(X1)
   n2=sum(X2)
   mle=1-sum(abs(w1*X1/n1-w2*X2/n2))
   p1hat=X1/n1
   p2hat=X2/n2
   boot.BC=rep(0,50)
   for(h in 1:50)
   {
     boot.X1=rmultinom(1,n1,p1hat)
     boot.X2=rmultinom(1,n2,p2hat)
     boot.BC[h]=sum(abs(w1*boot.X1/n1-w2*boot.X2/n2))
   }
   out=c(min(mle,1),sd(boot.BC));out=round(out,4)
   return(out)
}

KH_Braycurtis_equ=function(X1,X2,w1)
{
   BC=1-KH_Bray_curtis_equ(X1,X2,w1)
   n1=sum(X1);p1hat=X1/n1
   n2=sum(X2);p2hat=X2/n2
   boot.BC=rep(0,50)
   for(h in 1:50)
   {
       boot.X1=rmultinom(1,n1,p1hat)
       boot.X2=rmultinom(1,n2,p2hat)
       boot.BC[h]=KH_Bray_curtis_equ(boot.X1,boot.X2,w1)  
   }
   out=c(min(BC,1),sd(boot.BC));out=round(out,4)
   return(out)
}

entropy_MEE_equ=function(X)
{
  x=X
  x=x[x>0]
  n=sum(x)
  UE <- sum(x/n*(digamma(n)-digamma(x)))
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if(f1>0)
  {
     A <-1-ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
     B=sum(x==1)/n*(1-A)^(-n+1)*(-log(A)-sum(sapply(1:(n-1),function(k){1/k*(1-A)^k})))
  }
  if(f1==0){B=0}
  if(f1==1 & f2==0){B=0}
  UE+B
}
Two_com_correct_obspi=function(X1,X2)
{
   n1=sum(X1)
   n2=sum(X2)
   f11=sum(X1==1)
   f12=sum(X1==2)
   f21=sum(X2==1)
   f22=sum(X2==2)
   C1=1-f11/n1*(n1 - 1) * f11 / ((n1 - 1) * f11 + 2 * f12)
   C2=1-f21/n2*(n2 - 1) * f21 / ((n2 - 1) * f21 + 2 * f22)
   
   PP1=correct_obspi(X1)
   PP2=correct_obspi(X2)
   D12=which(X1>0 & X2>0)
  
   f0hat_1=ceiling( ifelse(f12 == 0,  f11 * (f11 - 1) / 2,  f11 ^ 2/ 2 / f12)   )   
   f0hat_2=ceiling( ifelse(f22 == 0,  f21 * (f21 - 1) / 2,  f21 ^ 2/ 2 / f22)   )
   #-----------------------------------------------------------------------------
 
   r1=which(X1>0 & X2==0)
   f.1=length(which(X1>0 & X2==1))
   f.2=length(which(X1>0 & X2==2))
   f.0=ceiling(ifelse(f.2>0,f.1^2/2/f.2,f.1*(f.1-1)/2))
   #------------------------------------------------------------------------------
   r2=which(X1==0 & X2>0)
   f1.=length(which(X1==1 & X2>0))
   f2.=length(which(X1==2 & X2>0))
   f0.=ceiling(ifelse(f2.>0,f1.^2/2/f2.,f1.*(f1.-1)/2))
   #------------------------------------------------------------------------------
   t11=length(which(X1==1 & X2==1))
   t22=length(which(X1==2 & X2==2))
   f00hat=ceiling( ifelse(t22 == 0,  t11 * (t11 - 1) / 4,  t11 ^ 2/ 4 / t22)   )
   #------------------------------------------------------------------------------
   temp1=max(length(r1),f.0)-length(r1)+f0.+f00hat
   temp2=max(length(r2),f0.)-length(r2)+f.0+f00hat
   p0hat_1=(1-C1)/max(f0hat_1,temp1)
   p0hat_2=(1-C2)/max(f0hat_2,temp2)
   #------------------------------------------------------------------------------
   P1=PP1[D12]
   P2=PP2[D12]
   if(length(r1)> f.0)
   {
      P1=c(P1,PP1[r1])
      Y=c(rep(p0hat_2, f.0), rep(0,length(r1)-f.0))
      P2=c(P2,sample(Y,length(Y)) )
   }
   if(length(r1)< f.0)
   {
      P1=c(P1,PP1[r1],rep( p0hat_1,f.0-length(r1)))
      P2=c(P2,rep(p0hat_2, f.0) )
   }
   #----------------------------------------------------------------------------   
   if(length(r2)> f0.)
   {
      Y=c(rep(p0hat_1,f0.),rep(0,length(r2)- f0.))
      P1=c(P1,sample(Y,length(Y)))
      P2=c(P2,PP2[r2] )
   }
   if(length(r2)< f0.)
   {
      P1=c(P1,rep(p0hat_1,f0.))
      P2=c(P2,PP2[r2],rep( p0hat_2,f0.-length(r2)) )
   }
   P1=c(P1,rep( p0hat_1,f00hat))
   P2=c(P2,rep( p0hat_2,f00hat))
   P1=c(P1, rep(p0hat_1,max( f0hat_1-temp1,0)) , rep(    0  ,max( f0hat_2-temp2,0))         )
   P2=c(P2, rep(     0 ,max( f0hat_1-temp1,0)) , rep(p0hat_2,max( f0hat_2-temp2,0))         ) 
   #------------------------------------------------------------------------------------
   a=cbind(P1,P2)
   return(a)
}


  
C1n_equ=function(method=c("relative","absolute"),X,boot=200)
{
  X=as.matrix(X)
  n=colSums(X)
  min_n=min(n)
  temp=sum(ifelse(n==min_n,0,1))
  no.community=length(X[1,])
  if(method=="relative")
  {
    if(temp>0)
    {   
       D1_alpha=sum( sapply(1:no.community,function(k) entropy_MEE_equ(X[,k])) )/no.community
       Horn=rep(0,200)
       rarefy.X=matrix(0,length(X[,1]),no.community)
       for(h in 1:200)
       {
           for(j in 1:no.community)
           {
               if(n[j]==min_n){rarefy.X[,j]=X[,j]}
               if(n[j]>min_n)
               {
                  Y=X[,j]
                  total=n[j]
                  y=0
                  for(i in 1:min_n)
                  {
                      P=Y/total
                      z=rmultinom(1,1,P)
                      Y=Y-z
                      total=total-1
                      y=y+z
                  }
                  rarefy.X[,j]=y
               } 
           }
           Horn[h]=( entropy_MEE_equ(rowSums(rarefy.X))-D1_alpha)/log(no.community)
       }
       C1n=1-mean(Horn)
       C1n_se=sd(Horn)
       a=c( C1n, C1n_se, max(0,C1n-1.96*C1n_se), min(1,C1n+1.96*C1n_se))
       return(a)
     }
     if(temp==0)
     {
        
        D1_alpha=sum( sapply(1:no.community,function(k) entropy_MEE_equ(X[,k])) )/(no.community)
        D1_gamma=entropy_MEE_equ(rowSums(X))
        Horn=(D1_gamma-D1_alpha)/log(no.community)
        if(no.community==2)
        {
           p_hat=Two_com_correct_obspi(X[,1],X[,2])
           boot.Horn=rep(0,boot)
           for(h in 1:boot)
           {
               boot.X=cbind(rmultinom(1,n[1],p_hat[,1]),rmultinom(1,n[2],p_hat[,2]) )
               boot.D1_alpha=sum( sapply(1:2,function(k) entropy_MEE_equ(boot.X[,k])) )/no.community
               boot.Horn[h]=( entropy_MEE_equ(rowSums(boot.X))-boot.D1_alpha)/log(no.community)
           }
           C1n=1-Horn
           C1n_se=sd(boot.Horn)
           a=c( C1n, C1n_se, max(0,C1n-1.96*C1n_se), min(1,C1n+1.96*C1n_se))
           return(a)
        }
        if(no.community>2)
        {
           boot.Horn=rep(0,boot)
           for(h in 1:boot)
           {
               boot.X=sapply(1:no.community,function(k)  rmultinom(1,n[k],X[,k]/n[k])   )
               boot.D1_alpha=sum( sapply(1:no.community,function(k) entropy_MEE_equ(boot.X[,k])) )/no.community
               boot.Horn[h]=( entropy_MEE_equ(rowSums(boot.X))-boot.D1_alpha)/log(no.community)
           }
           C1n=1-Horn
           C1n_se=sd(boot.Horn)
           a=c( C1n, C1n_se, max(0,C1n-1.96*C1n_se), min(1,C1n+1.96*C1n_se))
           return(a)
        }
      }
   }
   if(method=="absolute")
   {
      pool.size=sum(n)
      w=n/pool.size
      D1_alpha=sum( sapply(1:no.community,function(k) w[k]*entropy_MEE_equ(X[,k])) )
      D1_gamma=entropy_MEE_equ(rowSums(X))
      Horn=(D1_gamma-D1_alpha)/sum( -w*log(w))
      if(no.community==2)
      {
         p_hat=Two_com_correct_obspi(X[,1],X[,2])
         boot.Horn=rep(0,boot)
         for(h in 1:boot)
         {
             boot.X=cbind(rmultinom(1,n[1],p_hat[,1]),rmultinom(1,n[2],p_hat[,2]) )
             boot.D1_alpha=sum( sapply(1:2,function(k) w[k]*entropy_MEE_equ(boot.X[,k])) )
             boot.Horn[h]=( entropy_MEE_equ(rowSums(boot.X))-boot.D1_alpha)/sum( -w*log(w))
         }
         C1n=1-Horn
         C1n_se=sd(boot.Horn)
         a=c( C1n, C1n_se, max(0,C1n-1.96*C1n_se), min(1,C1n+1.96*C1n_se))
         return(a)
      }
      if(no.community>2)
      {
         boot.Horn=rep(0,boot)
         for(h in 1:boot)
         {
             boot.X=sapply(1:no.community,function(k)  rmultinom(1,n[k],X[,k]/n[k])   )
             boot.D1_alpha=sum( sapply(1:no.community,function(k) w[k]*entropy_MEE_equ(boot.X[,k])) )
             boot.Horn[h]=( entropy_MEE_equ(rowSums(boot.X))-boot.D1_alpha)/sum( -w*log(w))
         }
         C1n=1-Horn
         C1n_se=sd(boot.Horn)
         a=c( C1n, C1n_se, max(0,C1n-1.96*C1n_se), min(1,C1n+1.96*C1n_se))
         return(a)
      }
   }
}




print.spadeTwo <- function(x, ...){
	# abundance
	if(x$datatype=="abundance"){
    
	cat('(1) The loaded set includes abundance (or frequency) data from 2 assemblages\n')
    #cat('    (indexed by 1, 2 according to data entry order) and a total of',x$info1[1],'distinct species.\n\n')
    cat('(2) The two chosen assemblages for comparison: 1 vs. 2\n\n')
    cat('(3) Basic Data Information:\n\n')
    cat('         ( Number of observed individuals in Assemblage 1 )  n1=',x$info1[2],'\n')
    cat('         ( Number of observed individuals in Assemblage 2 )  n2=',x$info1[3],'\n')
    cat('             ( Number of observed species in Assemblage 1 )  D1=',x$info1[4],'\n')
    cat('             ( Number of observed species in Assemblage 2 )  D2=',x$info1[5],'\n')
    cat('   ( Number of observed shared species in two assemblages ) D12=',x$info1[6],'\n')
    cat('                      (Bootstrap replications for s.e. estimate)',x$info1[7],'\n\n')
    cat('    Some Statistics:\n')
    cat('          f[11]=',x$info2[1],'; f[1+]=',x$info2[2],
         '; f[+1]=',x$info2[3],'; f[2+]=',x$info2[4],
         '; f[+2]=',x$info2[5],'\n\n')
	
    cat('(4) Estimation Results of Some Similarity Indices: \n\n')
    cat('                                      Estimate       s.e.        U_hat* (s.e.)     V_hat** (s.e.)\n')
    cat('    Incidence-based:\n')
    cat('    ===============\n')
	temp <- apply(as.matrix(x$similarity), 2, as.numeric)
    cat('      Jaccard incidence  (observed)   ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'\n')
	  cat('      Jaccard incidence  (estimated)  ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'\n')
    cat('      Sorensen incidence (observed)   ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'\n')
	  cat('      Sorensen incidence (estimated)  ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'\n\n')
    #cat('      Lennon et al (2001)             ',sprintf("%.4f",temp[5,1]),'     ',sprintf("%.4f",temp[5,2]),'\n\n')
    cat('    Abundance-based:\n')
    cat('    ===============\n')
    #cat('      Bray-Curtis (observed)          ',sprintf("%.4f",temp[6,1]),'     ',sprintf("%.4f",temp[6,2]),'\n')
    cat('      Morisita-Horn                   ',sprintf("%.4f",temp[7,1]),'     ',sprintf("%.4f",temp[7,2]),'\n')
    cat('      Morisita Original               ',sprintf("%.4f",temp[8,1]),'     ',sprintf("%.4f",temp[8,2]),'\n') 
   # cat('      Horn (relative)                 ',sprintf("%.4f",1-temp[6,1]),'     ',sprintf("%.4f",temp[6,2]),'\n')
    #cat('      Horn (estimated)                 ',sprintf("%.4f",temp[9,1]),'     ',sprintf("%.4f",temp[9,2]),'\n')
    cat('      Jaccard Abundance  (unadjusted) ',sprintf("%.4f",temp[10,1]),'     ',sprintf("%.4f",temp[10,2]),'     ',sprintf("%.4f",temp[10,3]),'          ',sprintf("%.4f",temp[10,4]),'\n')
    cat('      Jaccard Abundance  (  adjusted) ',sprintf("%.4f",temp[11,1]),'     ',sprintf("%.4f",temp[11,2]),'     ',sprintf("%.4f",temp[11,3]),
               '(',sprintf("%.4f",temp[11,4]),')', sprintf("%.4f",temp[11,5]),'(',sprintf("%.4f",temp[11,6]),')','\n')
    cat('      Sorensen Abundance (unadjusted) ',sprintf("%.4f",temp[12,1]),'     ',sprintf("%.4f",temp[12,2]),'     ',sprintf("%.4f",temp[12,3]),
               '          ',sprintf("%.4f",temp[12,4]),'\n')
    cat('      Sorensen Abundance (  adjusted) ',sprintf("%.4f",temp[13,1]),'     ',sprintf("%.4f",temp[13,2]),
        '     ',sprintf("%.4f",temp[13,3]),'(',sprintf("%.4f",temp[13,4]),')', sprintf("%.4f",temp[13,5]),'(',sprintf("%.4f",temp[13,6]),')','\n')
    cat('     
      *  U denotes the total relative abundances of the shared species in the first assemblage;
           U_hat is an estimate of U.
      ** V denotes the total relative abundances of the shared species in the second assemblage;
           V_hat is an estimate of V.\n\n')
	  cat('    Equal-effort:\n')
   	cat('    ===============\n')
	  cat('      Bray-Curtis (observed)          ',sprintf("%.4f",temp[6,1]),'     ',sprintf("%.4f",temp[6,2]),'\n')
  	cat('      Bray-Curtis (estimated)         ',sprintf("%.4f",temp[15,1]),'     ',sprintf("%.4f",temp[15,2]),'\n')
	  cat('      Horn (observed)                 ',sprintf("%.4f",1-temp[14,1]),'     ',sprintf("%.4f",temp[14,2]),'\n')
	  cat('      Horn (estimated)                ',sprintf("%.4f",1-temp[9,1]),'     ',sprintf("%.4f",temp[9,2]),'\n')
	  cat('\n')
    cat('      In this part, we assume that the sample size ratio is equal to the population size ratio.')
	  cat('\n\n')
	  cat('    Equal-weight:\n')
	  cat('    ===============\n')
	  cat('      Bray-Curtis (observed)          ',sprintf("%.4f",temp[16,1]),'     ',sprintf("%.4f",temp[16,2]),'\n')
	  cat('      Bray-Curtis (estimated)         ',sprintf("%.4f",temp[17,1]),'     ',sprintf("%.4f",temp[17,2]),'\n')
  	cat('\n\n')
    cat('(5) References:\n')
    cat('
    Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2005). A new statistical approach for assessing
    similarity of species composition with incidence and abundance data. Ecology Letters, 8, 148-159.

    Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2006). Abundance-based similarity indices and 
    their estimation when there are unseen species in samples. Biometrics, 62, 361-371.\n\n')
    #cat('(6) ESTIMATION RESULTS OF THE NUMBER OF SPECIES FOR EACH ASSEMBLAGE:\n\n')
    #cat('     Model           Estimate      s.e.            95% CI\n\n')
    #cat('     Assemblage 1:\n')

    #temp=x$assemblage1
    #cat('       Chao1          ',sprintf("%.1f",temp[1,1]),'      ',sprintf("%.1f",temp[1,2]),'        (',sprintf("%.1f",temp[1,3]),',',sprintf("%.1f",temp[1,4]),')\n')
    #cat('       Chao1-bc       ',sprintf("%.1f",temp[2,1]),'      ',sprintf("%.1f",temp[2,2]),'        (',sprintf("%.1f",temp[2,3]),',',sprintf("%.1f",temp[2,4]),')\n')
    #cat('       ACE            ',sprintf("%.1f",temp[3,1]),'      ',sprintf("%.1f",temp[3,2]),'        (',sprintf("%.1f",temp[3,3]),',',sprintf("%.1f",temp[3,4]),')\n')
    #cat('       ACE-1          ',sprintf("%.1f",temp[4,1]),'      ',sprintf("%.1f",temp[4,2]),'        (',sprintf("%.1f",temp[4,3]),',',sprintf("%.1f",temp[4,4]),')\n')
    
	  #cat('     Assemblage 2:\n')
    #temp=x$assemblage2
    #cat('       Chao1          ',sprintf("%.1f",temp[1,1]),'      ',sprintf("%.1f",temp[1,2]),'        (',sprintf("%.1f",temp[1,3]),',',sprintf("%.1f",temp[1,4]),')\n')
    #cat('       Chao1-bc       ',sprintf("%.1f",temp[2,1]),'      ',sprintf("%.1f",temp[2,2]),'        (',sprintf("%.1f",temp[2,3]),',',sprintf("%.1f",temp[2,4]),')\n')
    #cat('       ACE            ',sprintf("%.1f",temp[3,1]),'      ',sprintf("%.1f",temp[3,2]),'        (',sprintf("%.1f",temp[3,3]),',',sprintf("%.1f",temp[3,4]),')\n')
    #cat('       ACE-1          ',sprintf("%.1f",temp[4,1]),'      ',sprintf("%.1f",temp[4,2]),'        (',sprintf("%.1f",temp[4,3]),',',sprintf("%.1f",temp[4,4]),')\n\n\n')
   # cat('(7) DESCRIPTIONS OF MODELS FOR SINGLE ASSEMBLAGE ANALYSIS:\n')
    #cat('
    # Chao1 (Chao, 1984): This approach uses the numbers of singletons and doubletons to
    # estimate the number of missing species because missing species information is
    # mostly concentrated on those low frequency counts; see Chao (1984), Shen, Chao and Lin (2003)
    # and Chao, Shen and Hwang (2006).

    # Chao1-bc: a bias-corrected form for the Chao1; see Chao (2005).

    # ACE (Abundance-based Coverage Estimator): A non-parametric estimator proposed by Chao and Lee (1992)
    # and Chao, Ma and Yang (1993).  The observed species are separated as rare and abundant groups;
    # only the rare group is used to estimate the number of missing species.
    # The estimated CV is used to characterize the degree of heterogeneity among species
    # discovery probabilities.  See Eq.(2.14) in Chao and Lee (1992) or Eq.(2.2) of Chao et al. (2000).

   #  ACE-1: A modified ACE for highly heterogeneous communities; See Eq.(2.15) of Chao and Lee (1992).
    # \n')
	}else if(x$datatype=="incidence"){
	cat('(1) The loaded set includes multiple-sample incidence data from 2 assemblages\n')
    #cat('    (indexed by 1, 2 according to data entry order) and a total of',x$info1[1],'distinct species.\n\n')
    cat('(2) The two chosen assemblages for comparison: 1 vs. 2\n\n')
    cat('(3) Basic Data Information:\n\n')
    cat('                    ( Number of samples from assemblage 1)   w=',x$info1[2],'\n')
    cat('                    ( Number of samples from assemblage 2)   z=',x$info1[3],'\n')
    cat('             ( Number of observed species in assemblage 1)  D1=',x$info1[4],'\n')
    cat('             ( Number of observed species in assemblage 2)  D2=',x$info1[5],'\n')
    cat('   ( Number of observed shared species in two assemblages ) D12=',x$info1[6],'\n')
    cat('                      (Bootstrap replications for s.e. estimate)',x$info1[7],'\n\n')
    cat('    Some Statistics:\n')
    cat('          Q[11]=',x$info2[1],
	            '; Q[1+]=',x$info2[2], '; Q[+1]=',x$info2[3],
				'; Q[2+]=',x$info2[4], '; Q[+2]=',x$info2[5],'\n\n')
    cat('(4) Estimation Results of Some Similarity Indices: \n\n')
    cat('                                            Estimate       s.e.        U_hat* (s.e.)     V_hat** (s.e.)\n')
    cat('    Incidence-based:\n')
    cat('    ===============\n')
    
	temp <- apply(as.matrix(x$similarity), 2, as.numeric)
    cat('      Jaccard incidence  (observed)         ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'\n')
	  cat('      Jaccard incidence  (estimated)        ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'\n')
    cat('      Sorensen incidence (observed)         ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'\n')
	  cat('      Sorensen incidence (estimated)        ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'\n')
    cat('      Lennon et al (2001)                   ',sprintf("%.4f",temp[5,1]),'     ',sprintf("%.4f",temp[5,2]),'\n\n')
    cat('    Multiple incidence-based:\n')
    cat('    ===============\n')
    cat('      Bray-Curtis                           ',sprintf("%.4f",temp[6,1]),'     ',sprintf("%.4f",temp[6,2]),'\n')
    cat('      Morisita-Horn                         ',sprintf("%.4f",temp[7,1]),'     ',sprintf("%.4f",temp[7,2]),'\n')
    cat('      Morisita Original                     ',sprintf("%.4f",temp[8,1]),'     ',sprintf("%.4f",temp[8,2]),'\n')
    cat('      Incidence-based Jaccard  (unadjusted) ',sprintf("%.4f",temp[9,1]),'     ',sprintf("%.4f",temp[9,2]),'     ',sprintf("%.4f",temp[9,3]),'          ',sprintf("%.4f",temp[9,4]),'\n')
    cat('      Incidence-based Jaccard  (  adjusted) ',sprintf("%.4f",temp[10,1]),'     ',sprintf("%.4f",temp[10,2]),'     ',sprintf("%.4f",temp[10,3]),
               '(',sprintf("%.4f",temp[10,4]),')', sprintf("%.4f",temp[10,5]),'(',sprintf("%.4f",temp[10,6]),')','\n')
    cat('      Incidence-based Sorensen (unadjusted) ',sprintf("%.4f",temp[11,1]),'     ',sprintf("%.4f",temp[11,2]),'     ',sprintf("%.4f",temp[11,3]),
               '          ',sprintf("%.4f",temp[11,4]),'\n')
    cat('      Incidence-based Sorensen (  adjusted) ',sprintf("%.4f",temp[12,1]),'     ',sprintf("%.4f",temp[12,2]),
        '     ',sprintf("%.4f",temp[12,3]),'(',sprintf("%.4f",temp[12,4]),')', sprintf("%.4f",temp[12,5]),'(',sprintf("%.4f",temp[12,6]),')','\n')
    cat('     
      *  U denotes the total relative incidences of the shared species in the first assemblage;
           U_hat is an estimate of U.
      ** V denotes the total relative incidences of the shared species in the second assemblage;
           V_hat is an estimate of V.\n\n')
    cat('(5) References:\n')
    cat('
    Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2005). A new statistical approach for assessing
    similarity of species composition with incidence and abundance data. Ecology Letters, 8, 148-159.

    Chao, A., Chazdon, R. L., Colwell, R. K. and Shen, T.-J. (2006). Abundance-based similarity indices and 
    their estimation when there are unseen species in samples. Biometrics, 62, 361-371.\n\n')
    cat('(6) ESTIMATION RESULTS OF THE NUMBER OF SPECIES FOR EACH ASSEMBLAGE:\n\n')
    cat('     Model           Estimate       s.e.            95% CI\n\n')
    cat('     Assemblage 1:\n')
    temp=x$assemblage1
    cat('       Chao2          ',sprintf("%.1f",temp[1,1]),'      ',sprintf("%.1f",temp[1,2]),'        (',sprintf("%.1f",temp[1,3]),',',sprintf("%.1f",temp[1,4]),')\n')
    cat('       Chao2-bc       ',sprintf("%.1f",temp[2,1]),'      ',sprintf("%.1f",temp[2,2]),'        (',sprintf("%.1f",temp[2,3]),',',sprintf("%.1f",temp[2,4]),')\n')
    cat('       ICE            ',sprintf("%.1f",temp[3,1]),'      ',sprintf("%.1f",temp[3,2]),'        (',sprintf("%.1f",temp[3,3]),',',sprintf("%.1f",temp[3,4]),')\n')
    cat('       ICE-1          ',sprintf("%.1f",temp[4,1]),'      ',sprintf("%.1f",temp[4,2]),'        (',sprintf("%.1f",temp[4,3]),',',sprintf("%.1f",temp[4,4]),')\n')
    cat('     Assemblage 2:\n')
    temp=x$assemblage2
    cat('       Chao2          ',sprintf("%.1f",temp[1,1]),'      ',sprintf("%.1f",temp[1,2]),'        (',sprintf("%.1f",temp[1,3]),',',sprintf("%.1f",temp[1,4]),')\n')
    cat('       Chao2-bc       ',sprintf("%.1f",temp[2,1]),'      ',sprintf("%.1f",temp[2,2]),'        (',sprintf("%.1f",temp[2,3]),',',sprintf("%.1f",temp[2,4]),')\n')
    cat('       ICE            ',sprintf("%.1f",temp[3,1]),'      ',sprintf("%.1f",temp[3,2]),'        (',sprintf("%.1f",temp[3,3]),',',sprintf("%.1f",temp[3,4]),')\n')
    cat('       ICE-1          ',sprintf("%.1f",temp[4,1]),'      ',sprintf("%.1f",temp[4,2]),'        (',sprintf("%.1f",temp[4,3]),',',sprintf("%.1f",temp[4,4]),')\n')
    cat('(7) DESCRIPTIONS OF MODELS FOR SINGLE ASSEMBLAGE ANALYSIS:\n')
    cat('
    Chao2 (Chao, 1987): This approach uses the frequencies of uniques and duplicates to
    estimate the number of missing species; see Chao (1987).

    Chao2-bc: a bias-corrected form for the Chao2; see Chao (2005).

    ICE (Incidence-based Coverage Estimator): Model(h) assumes that the detection
    probabilities are heterogeneous among species.  The estimator given here is an improved
    version of Eq.(3.18) in Lee and Chao (1994) by using an improved estimated sample coverage
    given in Shen (2003) and the SPADE User Guide; see Eq.(3.23) of Lee and Chao (1994)
    for the estimated squared CV.

    ICE-1:  A modified ICE for highly-heterogeneous cases.
     \n')
	}
}
