CV.Ind=function(x)
{
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  C.hat=1-f1/n
  Sobs=sum(x>0)
  S0=Sobs/C.hat
  r.square=max(S0*sum(x*(x-1))/n/(n-1)-1,0  )
  r.square^0.5
}
Chao1=function(x,conf=0.95)
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

Chao1_bc=function(x,conf=0.95)
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




EstiBootComm.Ind <- function(Spec)
{
	Sobs <- sum(Spec > 0) 	#observed species
	n <- sum(Spec)		  	#sample size
	f1 <- sum(Spec == 1) 	#singleton
	f2 <- sum(Spec == 2) 	#doubleton
	a <- ifelse(f1 == 0, 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n)
	b <- sum(Spec / n * (1 - Spec / n) ^ n)
	w <- a / b  			#adjusted factor for rare species in the sample
	f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))	#estimation of unseen species via Chao1
	Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
	Prob.hat.Unse <- rep(2 * f2/((n - 1) * f1 + 2 * f2), f0.hat)		#estimation of relative abundance of unseen species in the sample
	return(c(Prob.hat, Prob.hat.Unse))				#Output: a vector of estimated relative abundance
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
entropy_HT_equ<-function(X)
{
    x=X
    x=x[x>0]
    n=sum(x)
    f1=sum(x==1)
    C_head=1-f1/n
    a=-sum(C_head*(x/n)*log(C_head*(x/n))/(1-(1-C_head*(x/n))^n))
    a
}
entropy_J1_equ=function(X)
{
    X=X[X>0]
    Y=X[X>1]
    n=sum(X)
    -n*sum(X/n*log(X/n))-(n-1)/n*sum( (n-X)*(-X/(n-1)*log(X/(n-1))) )-(n-1)/n*sum(-Y*(Y-1)/(n-1)*log((Y-1)/(n-1)))
}
entropy_MLE_equ=function(X)
{
    X=X[X>0]
    n=sum(X)
    -sum(X/n*log(X/n))
}
entropy_MLE_bc_equ=function(X)
{
    entropy_MLE_equ(X)+(SpecAbunChao1(X,k=10,conf=0.95)[1]-1)/2/sum(X)
}
Shannon_index=function(x,boot=50)
{
  x=x[x>0]
  n=sum(x)
  MLE=entropy_MLE_equ(x)
  MLE_bc=entropy_MLE_bc_equ(x)
  J1=entropy_J1_equ(x)
  HT=entropy_HT_equ(x)
  MEE=entropy_MEE_equ(x)
  p_hat=EstiBootComm.Ind(x)
  Boot.X=rmultinom(boot,n,p_hat)
  temp1=apply(Boot.X,2,entropy_MLE_equ)
  temp2=apply(Boot.X,2,entropy_MLE_bc_equ)
  temp3=apply(Boot.X,2,entropy_J1_equ)
  temp4=apply(Boot.X,2,entropy_HT_equ)
  temp5=apply(Boot.X,2,entropy_MEE_equ)
     MLE_sd=sd(temp1)
  MLE_bc_sd=sd(temp2)
      J1_sd=sd(temp3)
      HT_sd=sd(temp4)
     MEE_sd=sd(temp5)

     MLE_exp_sd=sd(exp(temp1))
  MLE_bc_exp_sd=sd(exp(temp2))
      J1_exp_sd=sd(exp(temp3))
      HT_exp_sd=sd(exp(temp4))
     MEE_exp_sd=sd(exp(temp5))

  a=matrix(0,10,4)
  a[1,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
  a[2,]=c(MLE_bc,MLE_bc_sd,MLE_bc-1.96*MLE_bc_sd,MLE_bc+1.96*MLE_bc_sd)
  a[3,]=c(J1,J1_sd,J1-1.96*J1_sd,J1+1.96*J1_sd)
  a[4,]=c(HT,HT_sd,HT-1.96*HT_sd,HT+1.96*HT_sd)
  a[5,]=c(MEE,MEE_sd,MEE-1.96*MEE_sd,MEE+1.96*MEE_sd)
  a[6,]=c(exp(MLE),MLE_exp_sd,exp(MLE)-1.96*MLE_exp_sd,exp(MLE)+1.96*MLE_exp_sd)
  a[7,]=c(exp(MLE_bc),MLE_bc_exp_sd,exp(MLE_bc)-1.96*MLE_bc_exp_sd,exp(MLE_bc)+1.96*MLE_bc_exp_sd)
  a[8,]=c(exp(J1),J1_exp_sd,exp(J1)-1.96*J1_exp_sd,exp(J1)+1.96*J1_exp_sd)
  a[9,]=c(exp(HT),HT_exp_sd,exp(HT)-1.96*HT_exp_sd,exp(HT)+1.96*HT_exp_sd)
 a[10,]=c(exp(MEE),MEE_exp_sd,exp(MEE)-1.96*MEE_exp_sd,exp(MEE)+1.96*MEE_exp_sd)
  return(a)
}

simpson_MLE_equ=function(X)
{
   X=X[X>0]
   n=sum(X)
   a=sum((X/n)^2)
   a
}
simpson_MVUE_equ=function(X)
{
   X=X[X>0]
   n=sum(X)
   a=sum(X*(X-1))/n/(n-1)
   a
}

Simpson_index=function(x,boot=50)
{
   x=x[x>0]
   n=sum(x)
   MVUE=simpson_MVUE_equ(x)
   MLE=simpson_MLE_equ(x)

   #ACE=SpecAbunAce(x)[1]
   #AA=sum(  ( x*(x-1)/n/(n-1)-x*(2*n-1)/n/(n-1)*MVUE  )^2  )
   #BB=sum( x*(x-1)/n/(n-1)-x*(2*n-1)/n/(n-1)*MVUE   )
   #MVUE_sd=(AA-BB^2/ACE)^0.5

   #AA=sum(  ( (x/n)^2-2*x/n*MLE  )^2  )
   #BB=sum( (x/n)^2-2*x/n*MLE   )
   #MLE_sd=(AA-BB^2/ACE)^0.5
   #MVUE_recip_sd=MVUE_sd/MVUE
    #MLE_recip_sd=MLE_sd/MLE
   p_hat=EstiBootComm.Ind(x)
   Boot.X=rmultinom(boot,n,p_hat)
   temp1=apply(Boot.X,2,simpson_MVUE_equ)
   temp2=apply(Boot.X,2,simpson_MLE_equ)
   MVUE_sd=sd(temp1)
   MVUE_recip_sd=sd(1/temp1)
   MLE_sd=sd(temp2)
   MLE_recip_sd=sd(1/temp2)

   a=matrix(0,4,4)
   a[1,]=c(MVUE,MVUE_sd,MVUE-1.96*MVUE_sd,MVUE+1.96*MVUE_sd)
   a[2,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
   a[3,]=c(1/MVUE,MVUE_recip_sd,1/MVUE-1.96*MVUE_recip_sd,1/MVUE+1.96*MVUE_recip_sd)
   a[4,]=c(1/MLE,MLE_recip_sd,1/MLE-1.96*MLE_recip_sd,1/MLE+1.96*MLE_recip_sd)
   return(a)
}
######################################################2015.09.14
SpecInci <- function(data, k=10, conf=0.95)
{
  Chao2   <- SpecInciChao2(data, k = k, conf = conf)
  Chao2bc <- SpecInciChao2bc(data, k = k, conf = conf)
  iChao2  <- SpecInciiChao2(data, k=k, conf=conf)[1,]
  Modelh  <- SpecInciModelh(data, k = k, conf = conf)[-c(5)]
  Modelh1 <- SpecInciModelh1(data, k = k, conf = conf)[-c(5)]
  table   <- rbind(Chao2, Chao2bc, iChao2, Modelh, Modelh1)
  table   <- round(table,1)
  colnames(table) <- c("Estimate", "s.e.", "95%Lower", "95%Upper")
  rownames(table) <- c("Chao2 (Chao, 1987)", "Chao2-bc","iChao2", "ICE (Lee & Chao, 1994)", "ICE-1 (Lee & Chao, 1994)")
  return(table)
}

EstiBootComm.Sam <- function(data)
{
  data = data[data>0]
  T1 <- data[1]
  X = data[-c(1)]
  U = sum(X)   #observed species
  Q1 = length(X[X==1])   #singleton
  Q2 = length(X[X==2]) 	#doubleton
  if(Q2>0)
  {
    A = 2*Q2/((T1-1)*Q1+2*Q2)
  }
  else if(Q1>1)
  {
    A=2/((T1-1)*(Q1-1)+2)
  }else
  {
    A=0
  }
  C1 = 1 - Q1/U*(1 - A)
  W = U/T1*(1 - C1)/sum(X/T1*(1-X/T1)^T1)  			#adjusted factor for rare species in the sample
  Q0 = ceiling(ifelse(Q2>0, (T1-1)/T1*Q1^2/2/Q2, (T1-1)/T1*Q1*(Q1-1)/2))	#estimation of unseen species via Chao2
  Prob.hat = X/T1*(1-W*(1-X/T1)^T1)					#estimation of detection probability of observed species in the sample
  Prob.hat.Unse <- rep(U/T1*(1-C1)/Q0, Q0)		#estimation of detection probability of unseen species in the sample
  return(c(Prob.hat,  Prob.hat.Unse))									#Output: a vector of estimated detection probability
}
entropy_MLE_Inci_equ <- function(X)
{
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  H_MLE <- -sum(X/U*log(X/U))
  return(H_MLE)
}
entropy_MLE_bc_Inci_equ <- function(X)
{
  t <- X[1]
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  X_freq <- X[X > 10]
  X_infreq <- X[X <= 10]
  D_freq <- length(X_freq)
  D_infreq <- length(X_infreq)
  Q1 <- sum(X == 1)
  Q2 <- sum(X == 2)
  if(Q1 > 0 & Q2 > 0)
  {
    A <- 2*Q2/((t-1)*Q1 + 2*Q2)
  }
  else if (Q1 > 0 & Q2 == 0)
  {
    A <- 2/((t-1)*(Q1 - 1) + 2)
  }
  else
  {
    A <- 1
  }
  C_infreq <- 1 - Q1/sum(X_infreq)*(1-A)

  j <- c(1:10)
  b1 <- sum(sapply(j, function(j){j*(j-1)*sum(X == j)}))
  b2 <- sum(sapply(j, function(j){j*sum(X == j)}))
  gamma_infreq_square <- max(D_infreq/C_infreq*t/(t-1)*b1/b2/(b2-1) - 1, 0)

  ICE <- D_freq + D_infreq/C_infreq + Q1/C_infreq*gamma_infreq_square

  H_MLE <- -sum(X/U*log(X/U))
  H_MLE_bc <- H_MLE + (ICE/U + 1/t)/1

  return(H_MLE_bc)
}
entropy_HT_Inci_equ <- function(X)
{
  t <- X[1]
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  Q1 <- sum(X == 1)
  Q2 <- sum(X == 2)
  if(Q1 > 0 & Q2 > 0){
    A <- 2*Q2/((t-1)*Q1 + 2*Q2)
  } else if (Q1 > 0 & Q2 == 0){
    A <- 2/((t-1)*(Q1 - 1) + 2)
  } else {
    A <- 1
  }
  C <- 1 - Q1/U*(1-A)
  H_HT <- t/U*(-sum(C*X/t*log(C*X/t)/(1-(1-C*X/t)^t))) + log(U/t)
  return(H_HT)
}
entropy_MEE_Inci_equ <- function(X)
{
  t <- X[1]
  X <- X[-1]
  X <- X[X > 0]
  U <- sum(X)
  Q1 <- sum(X == 1)
  Q2 <- sum(X == 2)
  if(Q1 > 0 & Q2 > 0){
    A <- 2*Q2/((t-1)*Q1 + 2*Q2)
  } else if (Q1 > 0 & Q2 == 0){
    A <- 2/((t-1)*(Q1 - 1) + 2)
  } else {
    A <- 1
  }

  UE <- sum(X/t*(digamma(t)-digamma(X)))
  if(Q1 > 0 & A!=1){
    B <- Q1/t*(1-A)^(-t+1)*(-log(A)-sum(sapply(1:(t-1), function(k){1/k*(1-A)^k})))
    H_MEE <- t/U*(UE + B) + log(U/t)
  }else{
    H_MEE <- t/U*UE + log(U/t)
  }
  return(H_MEE)
}
Shannon_Inci_index=function(x,boot=50)
{
  x = unlist(x)
  t = x[1]
  MLE=entropy_MLE_Inci_equ(x)
  #MLE_bc=entropy_MLE_bc_Inci_equ(x)
  #HT=entropy_HT_Inci_equ(x)
  MEE=entropy_MEE_Inci_equ(x)
  p_hat=EstiBootComm.Sam(x)
  Boot.X = sapply(1:length(p_hat), function(i){
    rbinom(boot,t,p_hat[i])})
  Boot.X = cbind(rep(t,boot), Boot.X)
  temp1=apply(Boot.X,1,entropy_MLE_Inci_equ)
  #temp2=apply(Boot.X,1,entropy_MLE_bc_Inci_equ)
  #temp4=apply(Boot.X,1,entropy_HT_Inci_equ)
  temp5=apply(Boot.X,1,entropy_MEE_Inci_equ)
  MLE_sd=sd(temp1)
  #MLE_bc_sd=sd(temp2)
  #HT_sd=sd(temp4)
  MEE_sd=sd(temp5)

  MLE_exp_sd=sd(exp(temp1))
  #MLE_bc_exp_sd=sd(exp(temp2))
  #HT_exp_sd=sd(exp(temp4))
  MEE_exp_sd=sd(exp(temp5))

  a=matrix(0,8,4)
  a[1,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
  #a[2,]=c(MLE_bc,MLE_bc_sd,MLE_bc-1.96*MLE_bc_sd,MLE_bc+1.96*MLE_bc_sd)
  #a[3,]=c(HT,HT_sd,HT-1.96*HT_sd,HT+1.96*HT_sd)
  a[4,]=c(MEE,MEE_sd,MEE-1.96*MEE_sd,MEE+1.96*MEE_sd)
  a[5,]=c(exp(MLE),MLE_exp_sd,exp(MLE)-1.96*MLE_exp_sd,exp(MLE)+1.96*MLE_exp_sd)
  #a[6,]=c(exp(MLE_bc),MLE_bc_exp_sd,exp(MLE_bc)-1.96*MLE_bc_exp_sd,exp(MLE_bc)+1.96*MLE_bc_exp_sd)
  #a[7,]=c(exp(HT),HT_exp_sd,exp(HT)-1.96*HT_exp_sd,exp(HT)+1.96*HT_exp_sd)
  a[8,]=c(exp(MEE),MEE_exp_sd,exp(MEE)-1.96*MEE_exp_sd,exp(MEE)+1.96*MEE_exp_sd)
  return(a)
}
simpson_Inci_MVUE_equ=function(Y)
{
  t=Y[1]
  Y=Y[-1]
  Y=Y[Y>0]
  U=sum(Y)
  a=(sum(Y*(Y-1))/U^2/(1-1/t))
}
simpson_Inci_MLE_equ=function(Y)
{
  t=Y[1]
  Y=Y[-1]
  Y=Y[Y>0]
  a=(sum(Y^2)/sum(Y)^2)
}
Simpson_Inci_index=function(x,boot=200)
{
  x=x[x>0]
  t = x[1]
  MVUE=simpson_Inci_MVUE_equ(x)
  MLE=simpson_Inci_MLE_equ(x)

  p_hat=EstiBootComm.Sam(x)
  #set.seed(1)
  Boot.X = sapply(1:length(p_hat), function(i){
    rbinom(boot,t,p_hat[i])})
  Boot.X = cbind(rep(t,boot), Boot.X)
  temp1=apply(Boot.X,1,simpson_Inci_MVUE_equ)
  temp2=apply(Boot.X,1,simpson_Inci_MLE_equ)

  MVUE_sd=sd(temp1)
  MLE_sd=sd(temp2)

  #MVUE_recip_sd=MVUE_sd/MVUE
  #MLE_recip_sd=MLE_sd/MLE
  MVUE_recip_sd=sd(1/temp1)
  MLE_recip_sd=sd(1/temp2)

  a=matrix(0,4,4)
  a[1,]=c(MVUE,MVUE_sd,MVUE-1.96*MVUE_sd,MVUE+1.96*MVUE_sd)
  a[2,]=c(MLE,MLE_sd,MLE-1.96*MLE_sd,MLE+1.96*MLE_sd)
  a[3,]=c(1/MVUE,MVUE_recip_sd,1/MVUE-1.96*MVUE_recip_sd,1/MVUE+1.96*MVUE_recip_sd)
  a[4,]=c(1/MLE,MLE_recip_sd,1/MLE-1.96*MLE_recip_sd,1/MLE+1.96*MLE_recip_sd)
  return(a)
}
######################################################2015.09.14
conf.reg=function(x,LCL,UCL,...) polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)


#X=read.table("Data4a.txt")
#Y=read.table("Data4b1_t.txt")
#Diversity(datatype="Abundance",X)
#Diversity(datatype="Frequencies_of_Frequencies",Y)

print.spadeDiv <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

  if(x$datatype=="abundance"){

    cat("\n(1) BASIC DATA INFORMATION:\n")
    print(x$Basic_data)
    cat("\n(2) ESTIMATION OF SPECIES RICHNESS (DIVERSITY OF ORDER 0):\n\n")
    print(x$Species_richness)
    cat("
        Descriptions of richness estimators (See Species Part)
       ")
    cat("\n(3a) SHANNON ENTROPY:\n\n")
    print(x$Shannon_index)
    #cat("\n")
    #cat(" For a review of the four estimators, see Chao and Shen (2003).\n")
    #MLE_bc: bias-corrected empirical estimator.
    cat("
        MLE: empirical or observed entropy.
        Jackknife: see Zahl (1977).
        Chao & Shen: based on the Horvitz-Thompson estimator and sample coverage method; see Chao and Shen (2003).
        see Chao and Shen (2003).
  	    Chao et al. (2013): A nearly optimal estimator of Shannon entropy; see Chao et al. (2013).
  	    Estimated standard error is computed based on a bootstrap method.
        \n")

    cat("(3b) SHANNON DIVERSITY (EXPONENTIAL OF SHANNON ENTROPY):\n\n")
    print(x$Shannon_diversity)

    cat("\n(4a) SIMPSON CONCENTRATION INDEX:\n\n")
    print(x$Simpson_index)

    cat("
        MVUE: minimum variance unbiased estimator; see Eq. (2.27) of Magurran (1988).
        MLE: maximum likelihood estimator or empirical index; see Eq. (2.26) of Magurran (1988).
       ")

    cat("\n(4b) SIMPSON DIVERSITY (INVERSE OF SIMPSON CONCENTRATION):\n\n")
    print(x$Simpson_diversity)

    cat("\n(5) CHAO AND JOST (2015) ESTIMATES OF HILL NUMBERS \n\n")
    print(x$Hill_numbers)

    cat("
        ChaoJost: diversity profile estimator derived by Chao and Jost (2015).
  	    Empirical: maximum likelihood estimator (observed index).
       ")
  }else{
    cat("\n(1) BASIC DATA INFORMATION:\n")
    print(x$Basic_data)
    cat("\n(2) ESTIMATION OF SPECIES RICHNESS (DIVERSITY OF ORDER 0):\n\n")
    print(x$Species_richness)
    cat("
         Descriptions of richness estimators (See Species Part)
        ")
    cat("\n(3a) SHANNON INDEX:\n\n")
    print(x$Shannon_index)

    cat("\n(3b) EXPONENTIAL OF SHANNON INDEX (DIVERSITY OF ORDER 1):\n\n")
    print(x$Shannon_diversity)
    cat("\n(4a) SIMPSON INDEX:\n\n")
    print(x$Simpson_index)
    cat("\n(4b) INVERSE OF SIMPSON INDEX (DIVERSITY OF ORDER 2):\n\n")
    print(x$Simpson_diversity)
    cat("\n(5) Chao and Jost (2015) estimates of Hill numbers of order q from 0 to 3\n\n")
    print(x$Hill_numbers)

    cat("
        ChaoJost: diversity profile estimator derived by Chao and Jost (2015).
        Empirical: maximum likelihood estimator (observed index).
      ")

  }
  Lower=min(x$Hill_numbers[,3],x$Hill_numbers[,6])
  Upper=max(x$Hill_numbers[,4],x$Hill_numbers[,7])
  plot(0,type="n",xlim=c(min(x$Hill_numbers[,1]),max(x$Hill_numbers[,1])),ylim=c(Lower,Upper),xlab="Order  q",ylab="Hill  numbers")
  conf.reg(x$Hill_numbers[,1],x$Hill_numbers[,3],x$Hill_numbers[,4], col=adjustcolor(2, 0.2), border=NA)
  conf.reg(x$Hill_numbers[,1],x$Hill_numbers[,6],x$Hill_numbers[,7], col=adjustcolor(4, 0.2), border=NA)
  lines(x$Hill_numbers[,1],x$Hill_numbers[,2],col=2,lwd=3)
  lines(x$Hill_numbers[,1],x$Hill_numbers[,5],col=4,lty=3,lwd=3)
  legend("topright", c("ChaoJost","Empirical"),col=c(2,4),lwd=c(3,3),lty=c(1,3),bty="n",cex=0.8)
}



Chao_Hill_abu = function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))

  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      r <- 1:(n-1)
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum((1-p1)^r/r)))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,A^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      r <- 0:(n-1)
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}



Chao_Hill_inc = function(x,q){
  n = x[1]
  x = x[-1];x = x[x>0]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/U*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/U*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)*U/n
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,((n/U)^q*A)^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      ((n/U)^q*(A+B))^(1/(1-q))
    }
  }
  sapply(q, Sub)
}



Chao_Hill = function(x,q,datatype = c("abundance","incidence_freq")){
  datatype = match.arg(datatype,c("abundance","incidence_freq"))
  if(datatype == "abundance"){
    est = Chao_Hill_abu(x,q)
  }else{
    est = Chao_Hill_inc(x,q)
  }
  return(est)
}



Hill <- function(x,q,datatype = c("abundance","incidence_freq")){
  if(datatype=="incidence_freq"){x = x[-1]}
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}



Bt_prob_abu = function(x){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
  W = (1-C)/sum(x/n*(1-x/n)^n)

  p.new = x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0 = (1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}



Bt_prob_inc = function(x){
  n = x[1]
  x = x[-1]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  C=1-f1/U*(1-A)
  W=U/n*(1-C)/sum(x/n*(1-x/n)^n)

  p.new=x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0=U/n*(1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}



Bt_prob = function(x,datatype = c("abundance","incidence_freq")){
  datatype = match.arg(datatype,c("abundance","incidence_freq"))
  if(datatype == "abundance"){
    prob = Bt_prob_abu(x)
  }else{
    prob = Bt_prob_inc(x)
  }
  return(prob)
}


Bootstrap.CI = function(x,q,B = 200,datatype = c("abundance","incidence_freq"),conf = 0.95){
  datatype = match.arg(datatype,c("abundance","incidence_freq"))
  p.new = Bt_prob(x,datatype)
  n = ifelse(datatype=="abundance",sum(x),x[1])
  # set.seed(456)
  if(datatype=="abundance"){
    data.bt = rmultinom(B,n,p.new)
  }else{
    data.bt = rbinom(length(p.new)*B,n,p.new)
    data.bt = matrix(data.bt,ncol=B)
    #data.bt = rbind(rep(n,B),data.bt)
  }

  mle = apply(data.bt,2,function(x)Hill(x,q,datatype))

  ################ Abundance ###############
  if(datatype == "abundance"){
    d1 = sort(unique(data.bt[data.bt>0]))
    M = max(d1)
    Sub = function(q){
      d2 = sapply(1:length(d1),function(i){
        k=0:(n-d1[i])
        sum(choose(k-q,k)*exp(lchoose(n-k-1,d1[i]-1)-lchoose(n,d1[i])))
      })
      d3 = rep(0,M)
      d3[d1] = d2
      bt.pro = sapply(1:B,function(b){
        f1=sum(data.bt[,b]==1)
        f2=sum(data.bt[,b]==2)
        if(f2>0){
          A=2*f2/((n-1)*f1+2*f2)
        }else if(f1!=0){
          A=2/((n-1)*(f1-1)+2)
        }else{
          A=1
        }
        if(q!=1){
          t1=table(data.bt[,b][data.bt[,b]>0])
          t2=as.numeric(names(t1))
          aa=d3[t2]
          e1=sum(t1*aa)

          if(A==1){e2=0}else{
            r=0:(n-1)
            e2=f1/n*(1-A)^(-n+1)*(A^(q-1)-sum(choose(q-1,r)*(A-1)^r))
          }

          if(e1+e2!=0){
            e=(e1+e2)^(1/(1-q))
          }else{e=NA}

        }else{
          y2=data.bt[,b][which(data.bt[,b]>0 & data.bt[,b]<=(n-1))]
          e1=sum(y2/n*(digamma(n)-digamma(y2)))

          if(A==1){e2=0}else{
            r=1:(n-1)
            e2=f1/n*(1-A)^(-n+1)*(-log(A)-sum((1-A)^r/r))
          }
          e=exp(e1+e2)
        }
        e
      })
      bt.pro
    }
    pro = t(sapply(q,Sub))
  }else{
    d1 = sort(unique(data.bt[data.bt>0]))
    M = max(d1)
    Sub = function(q){
      d2 = sapply(1:length(d1),function(i){
        k = 0:(n-d1[i])
        sum(choose(k-q,k)*exp(lchoose(n-k-1,d1[i]-1)-lchoose(n,d1[i])))
      })
      d3 = rep(0,M)
      d3[d1] = d2

      bt.pro = sapply(1:B,function(b){
        y2=data.bt[,b];y2 = y2[y2>0]
        U = sum(y2)
        Q1=sum(y2==1)
        Q2=sum(y2==2)
        if(Q2>0){
          A=2*Q2/((n-1)*Q1+2*Q2)
        }else if(Q1!=0){
          A=2/((n-1)*(Q1-1)+2)
        }else{
          A=1
        }
        if(q!=1){
          t1=table(data.bt[,b][data.bt[,b]>0])
          t2=as.numeric(names(t1))
          aa=d3[t2]
          e1=sum(t1*aa)

          if(A==1){e2=0}else{
            r=0:(n-1)
            e2=Q1/n*(1-A)^(-n+1)*(A^(q-1)-sum(choose(q-1,r)*(A-1)^r))
          }

          if(e1+e2!=0){
            e=((n/U)^q*(e1+e2))^(1/(1-q))
          }else{e=NA}

        }else{

          e1=sum(y2/U*(digamma(n)-digamma(y2)))
          r =1:(n-1)
          e2 = ifelse(Q1==0|A==1,0,Q1/U*(1-A)^(1-n)*(-log(A)-sum((1-A)^r/r)))
          e = exp(e1+e2)*U/n
                  }
        e
      })
      bt.pro
    }


    pro = t(sapply(q,Sub))

  }


  #pro = apply(data.bt,2,function(x)Chao_Hill(x,q,datatype))

  mle.mean = rowMeans(mle)
  pro.mean = rowMeans(pro)

  LCI.mle =  -apply(mle,1,function(x)quantile(x,probs = (1-conf)/2)) + mle.mean
  UCI.mle = apply(mle,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - mle.mean

  LCI.pro =  -apply(pro,1,function(x)quantile(x,probs = (1-conf)/2)) + pro.mean
  UCI.pro = apply(pro,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - pro.mean

  LCI = rbind(LCI.mle,LCI.pro)
  UCI = rbind(UCI.mle,UCI.pro)

  sd.mle = apply(mle,1,sd)
  sd.pro = apply(pro,1,function(x)sd(x,na.rm = T))
  se = rbind(sd.mle,sd.pro)

  return(list(LCI=LCI,UCI=UCI,se=se))

}


ChaoHill <- function(dat, datatype=c("abundance", "incidence_freq"),q=NULL, from=0, to=3, interval=0.1, B=1000, conf=0.95){
  datatype = match.arg(datatype,c("abundance","incidence_freq"))
  # for real data estimation

  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)
  if(is.null(q)){q <- seq(from, to, by=interval)}
  if(!is.null(q)){q <- q}
  #-------------
  #Estimation
  #-------------
  MLE=Hill(dat,q,datatype)

  qD_pro=Chao_Hill(dat,q,datatype)

  CI_bound = Bootstrap.CI(dat,q,B,datatype,conf)
  se = CI_bound$se
  #-------------------
  #Confidence interval
  #-------------------
  tab.est=data.frame(rbind(MLE,qD_pro))

  LCI <- tab.est - CI_bound$LCI
  UCI <- tab.est + CI_bound$UCI

  colnames(tab.est) <- colnames(se) <- colnames(LCI) <- colnames(UCI) <- paste("q = ", q, sep="")
  rownames(tab.est) <- rownames(se) <- rownames(LCI) <- rownames(UCI) <- c("Observed", "Chao_2015")
  return(list(EST = tab.est,
              SD = se,
              LCI = LCI,
              UCI = UCI))

}




conf.reg=function(x_axis,LCL,UCL,...) {
  x.sort <- order(x_axis)
  x <- x_axis[x.sort]
  LCL <- LCL[x.sort]
  UCL <- UCL[x.sort]
  polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
}


reshapeChaoHill <- function(out){

  tab <- data.frame(q=as.numeric(substring(colnames(out$EST), 5)),
                    method=rep(rownames(out$EST), each=ncol(out$EST)),
                    est=c(t(out$EST)[,1],t(out$EST)[,2]),
                    se=c(t(out$SD)[,1],t(out$SD)[,2]),
                    qD.95.LCL=c(t(out$LCI)[,1],t(out$LCI)[,2]),
                    qD.95.UCL=c(t(out$UCI)[,1],t(out$UCI)[,2]))
  tab$est <- round(tab$est,3)
  tab$se <- round(tab$se,3)
  tab$qD.95.LCL <- round(tab$qD.95.LCL,3)
  tab$qD.95.UCL <- round(tab$qD.95.UCL,3)

  tab
}


Diversity_Inc=function(X)
{
  X=X[,1]

  Hill <- reshapeChaoHill(ChaoHill(X, datatype = "incidence_freq", from=0, to=3, interval=0.25, B=50, conf=0.95))
  #df$method <- factor(df$method, c("Observed", "Chao_2013"), c("Empirical", "Estimation"))
  Hill<-cbind(Hill[1:13,1],Hill[14:26,3],Hill[1:13,3],Hill[14:26,4],Hill[1:13,4])
  Hill<-round(Hill,3)
  Hill <- data.frame(Hill)
  colnames(Hill)<-c("q","Chao","Empirical","Chao(s.e.)","Empirical(s.e.)")


  z <- list("HILL.NUMBERS"= Hill)
  class(z) <- c("spadeDiv_Inc")
  return(z)

  #cat("\n")
  #cat("(5)  FISHER ALPHA INDEX:\n\n")
  #table_alpha=round(alpha(X),3)
  #colnames(table_alpha)<-c("Estimator", "Est_s.e.", paste("95% Lower Bound"), paste("95% Upper Bound"))
  #rownames(table_alpha)<-c(" alpha")
  #print( table_alpha)
  #cat("\n")
  #cat(" See Eq. (2.9) of Magurran (1988) for a definition of Fisher's alpha index.\n")
}
#X=read.table("Data4a.txt")
#Y=read.table("Data4b1_t.txt")
#Diversity(datatype="Abundance",X)
#Diversity(datatype="Frequencies_of_Frequencies",Y)

print.spadeDiv_Inc <- function(x, digits = max(3L, getOption("digits") - 3L), ...){


  cat("\n(5)  The estimates of Hill's number at order q from 0 to 3\n\n")
  print(x$HILL.NUMBERS)

  cat("
      Chao: see Chao and Jost (2015).
      Empirical: maximum likelihood estimator.
      ")

}
