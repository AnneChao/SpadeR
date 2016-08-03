###########################################2016.07.24-(P.L.Lin)
Horn_MLE_Multi = function(X, method=c("equal effort", "equal weight"))
{
  N <- ncol(X)
  n <- apply(X = X, MARGIN = 2, FUN = sum)
  p <- sapply(1:N, FUN = function(i){
    X[, i]/n[i]
  })
  if(method == "equal effort"){
    w <- n/sum(n)
  }else{w <- rep(1/N,N)}
  pbar <- w%*%t(p)
  pbar <- pbar[pbar>0]
  Hr <- sum(- pbar*log( pbar))
  Ha <- sum(sapply(1:N, FUN = function(i){
    p <- p[, i][p[, i]>0]
    w[i]*sum(-p*log(p))
  }))
  Ch <- 1-(Hr-Ha)/sum(-w*log(w))
  out=c(Ch)
  return(out)
}
Multi_Esti_GammaEntropy <- function(Mat, method="H-T")
{
  n <- apply(Mat, 2, sum)
  N <- ncol(Mat)
  pij <- sweep(Mat, 2, colSums(Mat), "/")
  pi. <- rowMeans(pij)
  
  if(method=="H-T"){
    Cj <- 1 - apply(Mat,2,function(x)sum(x==1)/sum(x))
    tmpFun <- function(Mat){
      tmp <- rowSums(Mat)==1
      if(sum(tmp)==0){
        rep(0,ncol(Mat))
      } else if(sum(tmp)==1){
        Mat[tmp,]
      } else {
        colSums(Mat[tmp,])
      }
    }
    
    Ct <- 1 - mean(tmpFun(Mat)/n)
    A <- matrix(0,ncol=ncol(Mat), nrow=nrow(Mat))
    for(j in 1:N){
      A[,j] <- (1-Cj[j]*pij[,j])^n[j]
    }
    
    -sum(Ct*pi.*log(Ct*pi.)/(1-apply(A,1,prod)), na.rm=TRUE)
  } else{
    pi. <- pi.[pi.>0]
    -sum(pi.*log(pi.))  #MLE
  }
}
Equal_weight_Horn_Esti_equ=function(X, datatype="abundance")
{
  nboot=50
  boot.Horn=rep(0,nboot)
  for(i in 1:nboot)
  {
    if(datatype=="abundance"){
      p <- Boots.pop(X)
      boot.X=sapply(1:dim(X)[2],function(k) rmultinom(1, sum(X[,k]), p[,k]))
    }else{
      p <- Boots.pop_inc(X)
      boot.X=sapply(1:dim(X)[2],function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, X[1,k], p[i,1]) ))
    }
    boot.Hr=Multi_Esti_GammaEntropy(boot.X, method="H-T")
    boot.Ha=mean(sapply(1:dim(X)[2],function(j)   entropy_MEE_equ(boot.X[,j])))
    boot.Horn[i]=1-(boot.Hr-boot.Ha)/log(dim(X)[2])
  }
  if(datatype=="incidence") X <- X[-1, ]
  Hr=Multi_Esti_GammaEntropy(X, method="H-T")
  Ha=mean(sapply(1:dim(X)[2],function(j)   entropy_MEE_equ(X[,j])))
  Horn=1-(Hr-Ha)/log(dim(X)[2])
  se_hat=sd(boot.Horn)
  out=c(Horn,se_hat,Horn-1.96*se_hat,Horn+1.96*se_hat)
  return(out)
}
Boots.pop=function(data)
{
  N=ncol(data);n=colSums(data);
  pool=rowSums(data);OBS=length(pool[pool>0]);
  data=data[pool>0,];
  obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]));
  F1=sum(pool==1);F2=sum(pool==2);
  F0=round(ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2)));
  
  f1=sapply(1:N,function(k) sum(data[,k]==1));
  f2=sapply(1:N,function(k) sum(data[,k]==2));
  C=sapply(1:N,function(k) 1-f1[k]/n[k]);
  
  f0=round(sapply(1:N,function(k) ifelse(f2[k]==0,f1[k]*(f1[k]-1)/2,f1[k]^2/(2*f2[k]))));
  r.data=sapply(1:N,function(k) data[,k]/n[k]);
  W=sapply(1:N,function(k) (1-C[k])/sum(r.data[,k]*(1-r.data[,k])^n[k]))
  
  if(F0>0){ boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=F0))
  }else{boots.pop=r.data}
  
  for(i in 1:N)
  {
    if(f0[i]>0)
    {
      f0[i]=ifelse(f0[i]+obs[i]>OBS+F0, OBS+F0-obs[i],f0[i])
      boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^n[i])   #
      I=which(boots.pop[,i]==0);II=sample(I,f0[i])
      boots.pop[II,i]=rep((1-C[i])/f0[i],f0[i])
    }
  }
  return(boots.pop)
}
Boots.pop_inc=function(data)
{ 
  data <- as.matrix(data)
  t=data[1,];Tt=sum(t);data=data[-1,];
  N=ncol(data);
  pool=rowSums(data);OBS=length(pool[pool>0]);
  data=data[pool>0,]; 
  obs=sapply(1:N,function(k) length(data[,k][data[,k]>0]));
  Q1=sum(pool==1);Q2=sum(pool==2);
  Q0=round(((Tt-1)/Tt)*ifelse(Q2==0,Q1*(Q1-1)/2,Q1^2/(2*Q2)));
  
  q1=sapply(1:N,function(k) sum(data[,k]==1));
  q2=sapply(1:N,function(k) sum(data[,k]==2));
  P1=sapply(1:N,function(k) ifelse(q1[k]+q2[k]==0,0,2*q2[k]/((t[k]-1)*q1[k]+2*q2[k])));
  
  q0=round(sapply(1:N,function(k) ((t[k]-1)/t[k])*ifelse(q2[k]==0,q1[k]*(q1[k]-1)/2,q1[k]^2/(2*q2[k]))));
  r.data=sapply(1:N,function(k) data[,k]/t[k]);
  W=sapply(1:N,function(k) (1-P1[k])*(q1[k]/t[k])/sum(r.data[,k]*(1-r.data[,k])^t[k]));
  
  if(Q0>0){ boots.pop=rbind(r.data,matrix(0,ncol=N,nrow=Q0))
  }else{boots.pop=r.data}
  
  for(i in 1:N){
    if(q0[i]>0){
      q0[i]=ifelse(q0[i]+obs[i]>OBS+Q0, OBS+Q0-obs[i],q0[i])
      boots.pop[,i][1:OBS]=boots.pop[,i][1:OBS]*(1-W[i]*(1-boots.pop[,i][1:OBS])^t[i])   #
      I=which(boots.pop[,i]==0);II=sample(I,q0[i])
      boots.pop[II,i]=rep((q1[i]/t[i])/q0[i],q0[i])
    }
  }
  return(boots.pop)
}
Horn_Multi_equ <- function(X, datatype="abundance", nboot=50,method=c("equal", "unequal"))
{
  boot=matrix(0,2,nboot)
  for(i in 1:nboot)
  {
    if(datatype=="abundance"){
      p <- Boots.pop(X)
      boot.X=sapply(1:dim(X)[2],function(k) rmultinom(1, sum(X[,k]), p[,k]))
    }else{
      p <- Boots.pop_inc(X)
      boot.X=sapply(1:dim(X)[2],function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, X[1,k], p[i,1]) ))
    }
    boot[,i]=Horn.Est(boot.X, method)
  }
  se <- apply(boot, MARGIN = 1, FUN = sd)
  if(datatype=="incidence") X <- X[-1, ]
  value <- Horn.Est(as.matrix(X), method)
  out <- c(value[1], se[1], max(0,value[1]-1.96*se[1]), min(1,value[1]+1.96*se[1]))
  out2 <- c(value[2], se[2],max(0,value[2]-1.96*se[2]), min(1,value[2]+1.96*se[2]))
  return(list(est=out,mle=out2))
  return(out)
}
SimilarityMul=function(X ,q, nboot=50, datatype="abundance", method=c("equal weight","unequal weight"))
{ 
  if(datatype=="incidence"){
    Y <- X ; X <- X[-1, ]
  }
  N=ncol(X);ni=colSums(X);n=sum(X);
  pool=rowSums(X);
  bX=apply(X,2,function(x) x/sum(x));pool.x=rowSums(bX)/N;
  
  if(q==0){
    f1=apply(X,2,function(x) sum(x==1));
    f2=apply(X,2,function(x) sum(x==2));
    Sobs=apply(X,2,function(x) sum(x>0)); 
    Si=Sobs+sapply(1:N, function(k) ifelse(f2[k]==0, f1[k]*(f1[k]-1)/2,f1[k]^2/(2*f2[k]))); 
    Sa=mean(Si);
    UqN.mle=(1/N-mean(Sobs)/sum(pool>0))/(1/N-1);
    CqN.mle=(N-sum(pool>0)/mean(Sobs))/(N-1);
    
    F1=sum(pool==1);F2=sum(pool==2);
    Sg=sum(pool>0)+ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2));
    UqN=min(1,(1/N-Sa/Sg)/(1/N-1));UqN=max(0,UqN);
    CqN=min(1,(N-Sg/Sa)/(N-1));CqN=max(0,CqN);
    
    b.UqN=numeric(nboot); b.UqN.mle=numeric(nboot);
    b.CqN=numeric(nboot); b.CqN.mle=numeric(nboot);
    for(i in 1:nboot){
      if(datatype=="abundance"){
        p <- Boots.pop(X)
        XX=sapply(1:N,function(k) rmultinom(1, ni[k], p[,k]))
      }else{
        p <- Boots.pop_inc(Y)
        XX=sapply(1:N,function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, Y[1,k], p[i,1]) ))
      }
      f1=apply(XX,2,function(x) sum(x==1));
      f2=apply(XX,2,function(x) sum(x==2));
      Sobs=apply(XX,2,function(x) sum(x>0)); 
      Si=Sobs+sapply(1:N,function(k) ifelse(f2[k]==0,f1[k]*(f1[k]-1)/2,f1[k]^2/(2*f2[k])))
      Sa=mean(Si);
      pool=rowSums(XX);
      b.UqN.mle[i]=(1/N-mean(Sobs)/sum(pool>0))/(1/N-1);
      b.CqN.mle[i]=(N-sum(pool>0)/mean(Sobs))/(N-1);
      
      F1=sum(pool==1);F2=sum(pool==2);
      Sg=sum(pool>0)+ifelse(F2==0,F1*(F1-1)/2,F1^2/(2*F2));
      b.UqN[i]=min(1,(1/N-Sa/Sg)/(1/N-1)); b.UqN[i]=max(0,b.UqN[i]);
      b.CqN[i]=min(1,(N-Sg/Sa)/(N-1));b.CqN[i]=max(0,b.CqN[i]);
    }
    se.U=sd(b.UqN);se.U.mle=sd(b.UqN.mle); #standard deviations of UqN est. and mle
    se.C=sd(b.CqN);se.C.mle=sd(b.CqN.mle); #standard deviations of  CqN est. and mle
    
    out1=rbind(c(UqN.mle,se.U.mle,min(1,UqN.mle+1.96*se.U.mle),max(0,UqN.mle-1.96*se.U.mle)),
               c(UqN,se.U,min(1,UqN+1.96*se.U),max(0,UqN-1.96*se.U)));
    out2=rbind(c(CqN.mle,se.C.mle,min(1,CqN.mle+1.96*se.C.mle),max(0,CqN.mle-1.96*se.C.mle)),
               c(CqN,se.C,min(1,CqN+1.96*se.C),max(0,CqN-1.96*se.C)));
  }
  
  if(q==2){
    if(method=="equal weight"){
      a.mle=N/sum(bX^2)
      g.mle=1/sum(pool.x^2);
      b.mle=g.mle/a.mle;
      UqN.mle=(N-b.mle)/(N-1);
      CqN.mle=(1/N-1/b.mle)/(1/N-1);
      
      Ai=sapply(1:N,function(k) sum(X[,k]*(X[,k]-1)/(ni[k]*(ni[k]-1))));
      bX.1=apply(X,2,function(x) (x-1)/(sum(x)-1));
      temp=sapply(1:nrow(X),function(j) (sum(bX[j,]%*%t(bX[j,]))-sum(bX[j,]^2))+sum(bX[j,]*bX.1[j,]));
      G=1/(sum(temp)/N^2);
      
      B=G/(1/mean(Ai));
      UqN=min(1,(N-B)/(N-1));UqN=max(0,UqN);
      CqN=min(1,(1/N-1/B)/(1/N-1));CqN=max(0,CqN);
      
      b.UqN=numeric(nboot);b.UqN.mle=numeric(nboot); b.CqN=numeric(nboot);b.CqN.mle=numeric(nboot);
      for(i in 1:nboot){
        if(datatype=="abundance"){
          p <- Boots.pop(X)
          XX=sapply(1:N,function(k) rmultinom(1, ni[k], p[,k]))
        }else{
          p <- Boots.pop_inc(Y)
          XX=sapply(1:N,function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, Y[1,k], p[i,1]) ))
        }
        bXX=apply(XX,2,function(x) x/sum(x));  
        pool.x=rowSums(bXX)/N;
        a.mle=N/sum(bXX^2);g.mle=1/sum(pool.x^2);b.mle=g.mle/a.mle;
        b.UqN.mle[i]=(N-b.mle)/(N-1);
        b.CqN.mle[i]=(1/N-1/b.mle)/(1/N-1);
        
        Ai=sapply(1:N,function(k) sum(XX[,k]*(XX[,k]-1)/(ni[k]*(ni[k]-1))));
        bXX.1=apply(XX,2,function(x) (x-1)/(sum(x)-1));
        temp=sapply(1:nrow(XX),function(j) (sum(bXX[j,]%*%t(bXX[j,]))-sum(bXX[j,]^2))+sum(bXX[j,]*bXX.1[j,]));
        G=1/(sum(temp)/N^2);
        
        B=G/(1/mean(Ai));
        b.UqN[i]=(N-B)/(N-1);
        b.CqN[i]=(1/N-1/B)/(1/N-1);
      }
      se.U=sd(b.UqN);se.U.mle=sd(b.UqN.mle);se.C=sd(b.CqN);se.C.mle=sd(b.CqN.mle);
      out1=rbind(c(UqN.mle,se.U.mle,min(1,UqN.mle+1.96*se.U.mle),max(0,UqN.mle-1.96*se.U.mle)),
                 c(UqN,se.U,min(1,UqN+1.96*se.U),max(0,UqN-1.96*se.U)));
      out2=rbind(c(CqN.mle,se.C.mle,min(1,CqN.mle+1.96*se.C.mle),max(0,CqN.mle-1.96*se.C.mle)),
                 c(CqN,se.C,min(1,CqN+1.96*se.C),max(0,CqN-1.96*se.C)));
    }
    if(method=="unequal weight"){
      a.mle=1/(N*sum((X/n)^2));g.mle=1/sum((pool/n)^2);b.mle=g.mle/a.mle;
      UqN.mle=(N-b.mle)/(N-1);
      CqN.mle=(1/N-1/b.mle)/(1/N-1);
      
      A=(1/N)*(1/sum(X*(X-1)/(n*(n-1))));
      G=1/sum(pool*(pool-1)/(n*(n-1)));
      B=G/A;
      UqN=min(1,(N-B)/(N-1));UqN=max(0,UqN);    
      CqN=min(1,(1/N-1/B)/(1/N-1));CqN=max(0,CqN);
      
      b.UqN=numeric(nboot);b.UqN.mle=numeric(nboot);b.CqN=numeric(nboot);b.CqN.mle=numeric(nboot);
      for(i in 1:nboot){
        if(datatype=="abundance"){
          p <- Boots.pop(X)
          XX=sapply(1:N,function(k) rmultinom(1, ni[k], p[,k]))
        }else{
          p <- Boots.pop_inc(Y)
          XX=sapply(1:N,function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, Y[1,k], p[i,1]) ))
        }
        pool=rowSums(XX);
        a.mle=1/(N*sum((XX/n)^2));g.mle=1/sum((pool/n)^2);b.mle=g.mle/a.mle;
        b.UqN.mle[i]=(N-b.mle)/(N-1);
        b.CqN.mle[i]=(1/N-1/b.mle)/(1/N-1);
        
        A=(1/N)*(1/sum(XX*(XX-1)/(n*(n-1))));
        G=1/sum(pool*(pool-1)/(n*(n-1)));
        B=G/A;
        b.UqN[i]=(N-B)/(N-1);
        b.CqN[i]=(1/N-1/B)/(1/N-1);
      }
      se.U=sd(b.UqN);se.U.mle=sd(b.UqN.mle);se.C=sd(b.CqN);se.C.mle=sd(b.CqN.mle);
      out1=rbind(c(UqN.mle,se.U.mle,min(1,UqN.mle+1.96*se.U.mle),max(0,UqN.mle-1.96*se.U.mle)),
                 c(UqN,se.U,min(1,UqN+1.96*se.U),max(0,UqN-1.96*se.U)));
      out2=rbind(c(CqN.mle,se.C.mle,min(1,CqN.mle+1.96*se.C.mle),max(0,CqN.mle-1.96*se.C.mle)),
                 c(CqN,se.C,min(1,CqN+1.96*se.C),max(0,CqN-1.96*se.C)));
    }
  }
  out1 <- cbind(out1[,c(1, 2)], out1[, 4], out1[, 3])
  colnames(out1)=c("UqN","se","95%.Lower","95%.Upper")
  rownames(out1)=c("Emperical","Estimate")
  out2 <- cbind(out2[,c(1, 2)], out2[, 4], out2[, 3])
  colnames(out2)=c("CqN","se","95%.Lower","95%.Upper")
  rownames(out2)=c("Emperical","Estimate")
  return(list(UqN=out1,CqN=out2));
}
Cq2_est_equ <- function(X, q, boot, datatype="abundance" ,method=c("equal effort", "equal weight"))
{ 
  N <- ncol(X)
  if(datatype=="abundance"){
    n <- apply(X = X, MARGIN = 2, FUN = sum)
  }else{
    n <- apply(X = X[-1,], MARGIN = 2, FUN = sum)
  }
  weight <- n/sum(n)
  weight <- - sum(weight*log(weight)) / log(N)
  plus_CI <-function(x){
    if(x[1] >= 1) x[1] <- 1
    if(x[1] <= 0) x[1] <- 0
    c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
  }
  if(q == 0){
    mat <- Jaccard_Sorensen_Abundance_equ(datatype, X[, 1], X[, 2], boot)[, c(1, 2)]
    out1 <- plus_CI(c(mat[4,1],mat[4,2]))
    out2 <- plus_CI(c(mat[2,1],mat[2,2]))
    out=rbind(out1, out2) 
  }
  if(method == "equal effort"){
    if(q == 1){
      out2 = Two_Horn_equ(X[, 1], X[, 2], weight = "unequal", datatype, method = "est", boot)
      out1 = plus_CI(c(weight*out2[1],out2[2]))
      out=rbind(out1, out2) 
    }
    if(q == 2){
      if(datatype=="incidence"){
        out=C2N_ee_se_inc(X, boot)
        out<-rbind(plus_CI(out[3,]),plus_CI(out[4,]))
      }else{
        out=SimilarityTwo(X, q, boot, datatype, method="unequal weight")
        out<-rbind(out$CqN[2,], out$UqN[2,])
      }
    }    
  }
  if(method == "equal weight"){
    if(q == 1){
      out = Two_Horn_equ(X[, 1], X[, 2], weight = "equal", datatype, method = "est", boot) ; out=rbind(out, out)
    }
    if(q == 2){
      out=SimilarityTwo(X, q, boot, datatype, method="equal weight") 
      out<-rbind(out$CqN[2,], out$UqN[2,])
    }
  }
  colnames(out) <- c("Est","se","95%.Lower","95%.Upper")
  return(out)
}
BC_equ <- function(X, datatype="abundance", nboot)
{
  boot=matrix(0,2,nboot)
  for(i in 1:nboot)
  {
    if(datatype=="abundance"){
      p <- Boots.pop(X)
      boot.X=sapply(1:dim(X)[2],function(k) rmultinom(1, sum(X[,k]), p[,k]))
    }else{
      p <- Boots.pop_inc(X)
      boot.X=sapply(1:dim(X)[2],function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, X[1,k], p[i,1]) ))
    }
    boot[,i]=BC.Est(boot.X)
  }
  se <- apply(boot, MARGIN = 1, FUN = sd)
  if(datatype=="incidence") X <- X[-1, ]
  value <- BC.Est(as.matrix(X))
  out <- c(value[1], se[1], max(0,value[1]-1.96*se[1]), min(1,value[1]+1.96*se[1]))
  out2 <- c(value[2], se[2],max(0,value[2]-1.96*se[2]), min(1,value[2]+1.96*se[2]))
  return(list(est=out,mle=out2))
}
###########################################
C0n_equ=function(X)
{
  X=as.matrix(X)
  I=which(X>0)
  X[I]=rep(1,length(I))
  no.community=length(X[1,])
  C0n_num=sum(rowSums(X)>0)
  C0n_dem=sum(X)
  C0n=no.community/(1-no.community)*( C0n_num-C0n_dem)/C0n_dem
  return( C0n)
}
correct_obspi<- function(X)
{
  Sobs <- sum(X > 0) 	
  n <- sum(X)		  	
  f1 <- sum(X == 1) 	
  f2 <- sum(X == 2)
  if(f1>0 & f2>0){a= (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n}
  if(f1>0 & f2==0){a= f1 / n*(n-1)*(f1-1)/((n-1)*(f1-1)+2)}
  if(f1==1 & f2==0){a=0}
  if(f1==0){a=0}
  b <- sum(X / n * (1 - X / n) ^ n)
  w <- a / b  			
  Prob.hat <- X / n * (1 - w * (1 - X / n) ^ n)	
  Prob.hat
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
Horn_MLE_Multi_equ <- function(X, datatype="abundance", nboot=50,method=c("equal effort", "equal weight"))
{
  boot.Horn=rep(0,nboot)
  for(i in 1:nboot)
  {
    if(datatype=="abundance"){
      p <- Boots.pop(X)
      boot.X=sapply(1:dim(X)[2],function(k) rmultinom(1, sum(X[,k]), p[,k]))
    }else{
      p <- Boots.pop_inc(X)
      boot.X=sapply(1:dim(X)[2],function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, X[1,k], p[i,1]) ))
    }
    boot.Horn[i]=Horn_MLE_Multi(boot.X, method)
  }
  if(datatype=="incidence") X <- X[-1, ]
  Horn=Horn_MLE_Multi(X, method)
  se_hat=sd(boot.Horn)
  out=c(Horn,se_hat,Horn-1.96*se_hat,Horn+1.96*se_hat)
  return(out)
}
C1n_equ=function(method=c("absolute","relative"), X, datatype="abundance" , boot=200)
{
  X=as.matrix(X)
  if(datatype=="incidence"){
    Y <- X ; X <- X[-1, ]
  }
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
        if(datatype=="abundance"){
          p <- Boots.pop(X)
          boot.X=sapply(1:dim(X)[2],function(k) rmultinom(1, sum(X[,k]), p[,k]))
        }else{
          p <- Boots.pop_inc(Y)
          boot.X=sapply(1:dim(Y)[2],function(k) sapply(1:nrow(p),FUN = function(i) rbinom(1, Y[1,k], p[i,1]) ))
        }
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
C2n_equ=function(method=c("MVUE","MLE"),X)
{ 
  X=as.matrix(X)
  no.community=length(X[1,])
  sobs=sum(rowSums(X)>0)
  n=colSums(X)
  temp1=X%*%matrix(c(1/n),no.community,1)
  temp2=sum(sapply(1:no.community,function(k) sum((X[,k]/n[k])^2) ))
  if(method=="MVUE")
  {
    c2n_dem=sum(sapply(1:no.community,function(k)  sum(X[,k]*(X[,k]-1)/n[k]/(n[k]-1))))
  }
  if(method=="MLE")
  {
    c2n_dem=sum(sapply(1:no.community,function(k)  sum(X[,k]*(X[,k])/n[k]/(n[k]))))
  }
  c2n_num=(sum(temp1^2)-temp2)/(no.community-1)
  c2n=c2n_num/c2n_dem
  return(c2n)
}
Cqn_se_equ=function(X,q=2,boot=200,method=c("relative","absolute"))
{
  X=as.matrix(X)
  no.community=length(X[1,])
  n=colSums(X)
  p_hat=sapply(1:no.community,function(j) X[,j]/n[j])
  
  if(q==0)
  {
    C0n=C0n_equ(X)
    boot.C0n=rep(0,boot)
    for(h in 1:boot)
    {
      boot.X=sapply(1:no.community,function(j) rmultinom(1,n[j],p_hat[,j]) )
      boot.C0n[h]=C0n_equ(boot.X)
    }
    C0n_se=sd(boot.C0n)
    a=c(C0n,C0n_se,max(0,C0n-1.96*C0n_se),min(1,C0n+1.96*C0n_se) )
    return(a)
  }
  if(q==1)
  {
    a=C1n_equ(method,X,boot)
    return(a)
  }
  if(q==2)
  {
    C2n=C2n_equ(method="MVUE",X)
    boot.C2n=rep(0,boot)
    if(C2n>1)
    {
      C2n=C2n_equ(method="MLE",X)       
      for(h in 1:boot)
      {
        boot.X=sapply(1:no.community,function(j) rmultinom(1,n[j],p_hat[,j]) )
        boot.C2n[h]=C2n_equ(method="MLE",boot.X)
      }
    }
    if(C2n<=1)
    {
      for(h in 1:boot)
      {
        boot.X=sapply(1:no.community,function(j) rmultinom(1,n[j],p_hat[,j]) )
        boot.C2n[h]=C2n_equ(method="MVUE",boot.X)
      }
    }
    C2n_se=sd(boot.C2n)
    Bootmean=mean(boot.C2n)
    a=c(C2n,C2n_se,max(0,C2n-Bootmean+quantile(boot.C2n, probs = 0.025)),min(1,C2n-Bootmean+quantile(boot.C2n, probs = 0.975)),
        max(0,1-C2n-(1-Bootmean)+(1-quantile(boot.C2n, probs = 0.975))),min(1,1-C2n-(1-Bootmean)+(1-quantile(boot.C2n, probs = 0.025)))
    )
    return(a)
  }
}
C33_equ=function(method=c("MVUE","MLE"),X)
{
  X=as.matrix(X)
  n=colSums(X)
  if(method=="MVUE")
  {
    C33_num=sum(  3*X[,1]*(X[,1]-1)/n[1]/(n[1]-1)*(X[,2]/n[2]+X[,3]/n[3])+    
                    3*X[,2]*(X[,2]-1)/n[2]/(n[2]-1)*(X[,1]/n[1]+X[,3]/n[3])+
                    3*X[,3]*(X[,3]-1)/n[3]/(n[3]-1)*(X[,1]/n[1]+X[,2]/n[2])+
                    6*X[,1]/n[1]*X[,2]/n[2]*X[,3]/n[3])/8
    C33_dem=sum(sapply(1:3,function(k)  sum( X[,k]*(X[,k]-1)*(X[,k]-2)/n[k]/(n[k]-1)/(n[k]-2))))
  }
  if(method=="MLE")
  {
    C33_num=sum(  3*X[,1]*(X[,1])/n[1]/(n[1])*(X[,2]/n[2]+X[,3]/n[3])+    
                    3*X[,2]*(X[,2])/n[2]/(n[2])*(X[,1]/n[1]+X[,3]/n[3])+
                    3*X[,3]*(X[,3])/n[3]/(n[3])*(X[,1]/n[1]+X[,2]/n[2])+
                    6*X[,1]/n[1]*X[,2]/n[2]*X[,3]/n[3])/8
    C33_dem=sum(sapply(1:3,function(k)  sum( X[,k]*(X[,k])*(X[,k])/n[k]/(n[k])/(n[k]))))
  }
  return(C33_num/C33_dem)                      
}
C33_se_equ=function(X,boot=200)
{
  X=as.matrix(X)
  no.community=length(X[1,])
  C33=C33_equ(method="MVUE",X)
  n=colSums(X)
  p_hat=sapply(1:no.community,function(j) X[,j]/n[j])
  boot.C33=rep(0,boot)
  if(C33>1)
  {
    C33=C33_equ(method="MLE",X)       
    for(h in 1:boot)
    {
      boot.X=sapply(1:no.community,function(j) rmultinom(1,n[j],p_hat[,j]) )
      boot.C33[h]=C33_equ(method="MLE",boot.X)
    }
  }
  if(C33<=1)
  {
    for(h in 1:boot)
    {
      boot.X=sapply(1:no.community,function(j) rmultinom(1,n[j],p_hat[,j]) )
      boot.C33[h]=C33_equ(method="MVUE",boot.X)
    }
  }
  C33_se=sd(boot.C33)
  Bootmean=mean(boot.C33)
  a=c(C33,C33_se,max(0,C33-Bootmean+quantile(boot.C33, probs = 0.025)),min(1,C33-Bootmean+quantile(boot.C33, probs = 0.975)),
      max(0,1-C33-(1-Bootmean)+(1-quantile(boot.C33, probs = 0.975))),min(1,1-C33-(1-Bootmean)+(1-quantile(boot.C33, probs = 0.025)))
  )
  return(a) 
}

#Multiple_Community_Measure=function(X,q=2,boot=200,method=c("relative","absolute"))

print.spadeMult <- function(x, ...){
  if(x$datatype=="abundance"){
    cat('\n(1) BASIC DATA INFORMATION:\n\n')
    cat('    The loaded set includes abundance/incidence data from',x$info[1],'communities\n')
    cat('    and a total of',x$info[2],'species.\n\n')
    cat('    Sample size in each community                            n1   =', x$info[3],'\n')
    N <- x$info[1]
    q <- x$q
    method <- x$goal
    for(j in 2:N){
      cat('                                                            ','n')
      cat(j,'  =',x$info[2+j],'\n')   
    }
    cat('\n')
    cat('    Number of observed species in one community              D1   =', x$info[N+3],'\n')
    
    for(j in 2:N){
      cat('                                                            ','D')
      cat(j,'  =',x$info[N+2+j],'\n')
    }
    cat('\n')
    cat('    Number of observed shared species in two communities     D12  =', x$info[3+2*N], '\n')
    
    if(N>2){
      k <- 1
      for(i in 1:(N-1)){     
        for(j in (i+1):N){
          if(i==1 & j==2) next
          cat('                                                            ','D')
          cat(i,j,'  = ', x$info[3+2*N+k], '\n', sep="")
          k <- k + 1
        }
      }
    }
    cat('\n')
    if(N==3)
    {
      cat('    Number of observed shared species in three communities   D123 =',rev(x$info)[2],'\n\n') 
    }
    cat('    Number of bootstrap replications for s.e. estimate             ',rev(x$info)[1],'\n\n')
    cat('(2) EMPIRICAL SIMILARITY INDICES: \n\n')
    cat('                                       Estimate       s.e.       95%Lower     95%Upper\n')
    cat('    (a) Classic richness-based similarity\n\n')
    temp <- apply(as.matrix(x$Empirical_richness), 2, as.numeric)
    cat('        C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')
    cat('        U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')
    cat('    (b) Measures for comparing species relative abundances\n\n')
    temp <- apply(as.matrix(x$Empirical_relative), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1, Horn)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (q=2, Regional overlap)    ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('    (c) Measures for comparing size-weighted species relative abundances\n\n')
    temp <- x$Empirical_WtRelative
    cat('        Horn size-weighted (q=1)       ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('    (d) Measures for comparing species absolute abundances\n\n')
    temp <- apply(as.matrix(x$Empirical_absolute), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1)                  ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (Morisita-Horn)            ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (Regional overlap)         ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('        Bray-Curtis                    ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
    cat('(3) ESTIMATED SIMILARITY INDICES: \n\n')
    cat('                                       Estimate       s.e.       95%Lower     95%Upper\n')
    cat('    (a) Classic richness-based similarity\n\n')
    temp <- apply(as.matrix(x$estimated_richness), 2, as.numeric)
    if(temp[1,1]>1) {cat('        C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",1),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
    if(temp[1,1]<=1){cat('        C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
    if(temp[2,1]>1) {cat('        U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",1) ,'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
    if(temp[2,1]<=1){cat('        U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
    #cat('      Lennon et al (2001)             ',sprintf("%.4f",temp[5,1]),'     ',sprintf("%.4f",temp[5,2]),'\n\n')
    cat('    (b) Measures for comparing species relative abundances\n\n')
    temp <- apply(as.matrix(x$estimated_relative), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1, Horn)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (q=2, Regional overlap)    ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('    (c) Measures for comparing size-weighted species relative abundances\n\n')
    temp <- x$estimated_WtRelative
    cat('        Horn size-weighted (q=1)       ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('    (d) Measures for comparing species absolute abundances\n\n')
    temp <- apply(as.matrix(x$estimated_absolute), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1)                  ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (Morisita-Horn)            ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (Regional overlap)         ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('        Bray-Curtis                    ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
   
  }else{
    
    cat('\n(1) BASIC DATA INFORMATION:\n\n')
    cat('    The loaded set includes abundance/incidence data from',x$info[[1]],'communities\n')
    cat('    and a total of',x$info[2],'species.\n\n')
    
    cat('    Number of sample units in each community                 T1   =', x$info[3],'\n')
    N <- x$info[1]
    q <- x$q
    method <- x$goal
    for(j in 2:N){
      cat('                                                            ','T')
      cat(j,'  =',x$info[2+j],'\n')   
    }
    cat('\n')
    cat('    Number of total incidences in each community             U1   =', x$info[N+3],'\n')
    N <- x$info[1]
    q <- x$q
    method <- x$goal
    for(j in 2:N){
      cat('                                                            ','U')
      cat(j,'  =',x$info[N+2+j],'\n')   
    }
    cat('\n')
    cat('    Number of observed species in one community              D1   =', x$info[2*N+3],'\n')
    for(j in 2:N){
      cat('                                                            ','D')
      cat(j,'  =',x$info[2*N+2+j],'\n')
    }
    cat('\n')
    cat('    Number of observed shared species in two communities     D12  =', x$info[3+3*N], '\n')
    
    if(N>2){
      k <- 1
      for(i in 1:(N-1)){     
        for(j in (i+1):N){
          if(i==1 & j==2) next
          cat('                                                            ','D')
          cat(i,j,'  = ', x$info[3+3*N+k], '\n', sep="")
          k <- k + 1
        }
      }
    }
    cat('\n')
    if(N==3)
    {
      cat('    Number of observed shared species in three communities   D123 =',rev(x$info)[2],'\n\n') 
    }
    cat('    Number of bootstrap replications for s.e. estimate             ',rev(x$info)[1],'\n\n')
    cat('(2) EMPIRICAL SIMILARITY INDICES: \n\n')
    cat('                                       Estimate       s.e.       95%Lower     95%Upper\n')
    cat('    (a) Classic richness-based similarity\n\n')
    temp <- apply(as.matrix(x$Empirical_richness), 2, as.numeric)
    cat('        C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')
    cat('        U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')
    cat('    (b) Measures for comparing species relative abundances\n\n')
    temp <- apply(as.matrix(x$Empirical_relative), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1, Horn)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (q=2, Regional overlap)    ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('    (c) Measures for comparing size-weighted species relative abundances\n\n')
    temp <- x$Empirical_WtRelative
    cat('        Horn size-weighted (q=1)       ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('    (d) Measures for comparing species absolute abundances\n\n')
    temp <- apply(as.matrix(x$Empirical_absolute), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1)                  ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (Morisita-Horn)            ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (Regional overlap)         ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('        Bray-Curtis                    ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
    cat('(3) ESTIMATED SIMILARITY INDICES: \n\n')
    cat('                                       Estimate       s.e.       95%Lower     95%Upper\n')
    cat('    (a) Classic richness-based similarity\n\n')
    temp <- apply(as.matrix(x$estimated_richness), 2, as.numeric)
    if(temp[1,1]>1) {cat('        C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",1),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
    if(temp[1,1]<=1){cat('        C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
    if(temp[2,1]>1) {cat('        U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",1) ,'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
    if(temp[2,1]<=1){cat('        U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
    #cat('      Lennon et al (2001)             ',sprintf("%.4f",temp[5,1]),'     ',sprintf("%.4f",temp[5,2]),'\n\n')
    cat('    (b) Measures for comparing species relative abundances\n\n')
    temp <- apply(as.matrix(x$estimated_relative), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1, Horn)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (q=2, Regional overlap)    ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('    (c) Measures for comparing size-weighted species relative abundances\n\n')
    temp <- x$estimated_WtRelative
    cat('        Horn size-weighted (q=1)       ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('    (d) Measures for comparing species absolute abundances\n\n')
    temp <- apply(as.matrix(x$estimated_absolute), 2, as.numeric)
    cat('        C1');cat(N);cat('=U1');cat(N);cat(' (q=1)                  ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
    cat('        C2');cat(N);cat(' (Morisita-Horn)            ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
    cat('        U2');cat(N);cat(' (Regional overlap)         ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
    cat('        Bray-Curtis                    ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
  }
  cat('(4) ESTIMATED PAIRWISE SIMILARITY:\n\n')
  if(q == 0){
    cat('    -------------------------Measure C02---------------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    Cqn_PC <- x$pairwise$C02
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    C02(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','     ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise similarity=',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise similarity matrix: \n\n')
    C_SM=x$similarity.matrix$C02
    cat('    C02(i,j)  ')
    for(i in 1:N)
    {
      cat(i,"      ")
    }
    cat('\n')
    for(i in 1:N)
    {
      cat('      ',i,'     ')
      for(j in 1:N)
      {
        if(i>j){cat('        ')}
        if(i<=j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    cat('    -------------------------Measure U02---------------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    Cqn_PC <- x$pairwise$U02
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    U02(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','     ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    
    cat('\n')
    cat('    Average pairwise similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise similarity matrix: \n\n')
    C_SM=x$similarity.matrix$U02
    cat('    U02(i,j)  ')
    for(i in 1:N)
    {
      cat(i,"      ")
    }
    cat('\n')
    for(i in 1:N)
    {
      cat('      ',i,'     ')
      for(j in 1:N)
      {
        if(i>j){cat('        ')}
        if(i<=j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
  }
  if(q == 1){
    cat('    ----------------------Measure C12 (=U12)------------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    ###################################################CqN_ Equal weight
    Cqn_PC <- x$pairwise$C12
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    C12(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','     ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise similarity=',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise similarity matrix: \n\n')
    C_SM=x$similarity.matrix$C12
    cat('    C12(i,j)  ')
    for(i in 1:N)
    {
      cat(i,"      ")
    }
    cat('\n')
    for(i in 1:N)
    {
      cat('      ',i,'     ')
      for(j in 1:N)
      {
        if(i>j){cat('        ')}
        if(i<=j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    if(method=="relative"){
      cat('    -------------------Measure Horn size-weighted--------------------\n\n')
      cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
      Cqn_PC <- x$pairwise$Horn
      no.temp=1
      for(i in 1:(N-1))
      {
        for(j in (i+1):N)
        {
          temp=Cqn_PC[no.temp,]
          cat('    Horn(')
          cat(i)
          cat(',')
          cat(j)
          if(temp[1]>1)
          {cat(')','    ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
          }
          if(temp[1]<=1)
          {cat(')','    ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
          no.temp=no.temp+1
        }
      }
      
      cat('\n')
      cat('    Average pairwise similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
      cat('    Pairwise similarity matrix: \n\n')
      C_SM=x$similarity.matrix$Horn
      
      cat('    Horn(i,j) ')
      for(i in 1:N)
      {
        cat(i,"      ")
      }
      cat('\n')
      for(i in 1:N)
      {
        cat('      ',i,'     ')
        for(j in 1:N)
        {
          if(i>j){cat('        ')}
          if(i<=j) {
            if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
            else cat(sprintf("%.3f#",1),'  ')
          }
        }
        cat('\n')
      }
      cat('\n')
    }
  }
  if(q == 2){
    cat('    -------------------------Measure C22---------------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    Cqn_PC <- x$pairwise$C22
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    C22(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','     ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise similarity=',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise similarity matrix: \n\n')
    C_SM=x$similarity.matrix$C22
    cat('    C22(i,j)  ')
    for(i in 1:N)
    {
      cat(i,"      ")
    }
    cat('\n')
    for(i in 1:N)
    {
      cat('      ',i,'     ')
      for(j in 1:N)
      {
        if(i>j){cat('        ')}
        if(i<=j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    ###################################################UqN_ Equal weight
    cat('    -------------------------Measure U22---------------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    Cqn_PC <- x$pairwise$U22
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    U22(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','     ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise similarity matrix: \n\n')
    C_SM=x$similarity.matrix$U22
    cat('    U22(i,j)  ')
    for(i in 1:N)
    {
      cat(i,"      ")
    }
    cat('\n')
    for(i in 1:N)
    {
      cat('      ',i,'     ')
      for(j in 1:N)
      {
        if(i>j){cat('        ')}
        if(i<=j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
  }
  
  cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
  #cat('
  #    References:
  #    
  #    Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-
  #    stage probabilistic approach to multiple-community similarity indices. 
  #    Biometrics, 64, 1178-1186.
  #    
  #    Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular
  #    Ecology, 17, 4015-4026.
  #    ')
  cat('\n')
}




#print.spadeMult <- function(x, ...){
#  cat('\n(1) BASIC DATA INFORMATION:\n\n')
#  cat('    The loaded set includes abundance (or frequency) data from',x$info[1],'communities\n')
#  cat('    and a total of',x$info[2],'distinct species.\n\n')
#  cat('   (Number of observed individuals in each community)       n1 =', x$info[3],'\n')
#  N <- x$info[1]
#  q <- x$q
#  for(j in 2:N){
#    cat('                                                           ','n')
#    cat(j,'=',x$info[2+j],'\n')   
#  }
#  cat('\n')
# cat('   (Number of observed species in one community)            D1 =', x$info[N+3],'\n')
#  
#  for(j in 2:N){
#    cat('                                                           ','D')
#    cat(j,'=',x$info[N+2+j],'\n')
#  }
#  cat('\n')
#  cat('   (Number of observed shared species in two communities)  ','D12 =', x$info[3+2*N], '\n')
#  
#  if(N>2){
#    k <- 1
#    for(i in 1:(N-1)){     
#      for(j in (i+1):N){
#        if(i==1 & j==2) next
#        cat('                                                           ','D')
#        cat(i,j,' = ', x$info[3+2*N+k], '\n', sep="")
#        k <- k + 1
#      }
#    }
#  }
#  cat('\n')
#  if(N==3)
#  {
#    cat('   (Number of observed shared species in three communities)','D123','=',rev(x$info)[2],'\n\n') 
#  }
#  cat('   (Bootstrap replications for s.e. estimate)              ',rev(x$info)[1],'\n\n')
#  cat('(2) ESTIMATION OF OVERLAP MEASURE IN',N,'COMMUNITIES:\n\n')
#  cat('    Estimator','     Estimate','     s.e.','     95% Confidence Interval\n\n')
#  temp0n=x$overlap[1,]
#  if(temp0n[1]>1)
#  {
#    cat('    C0')
#    cat(N,'(Sorensen)',sprintf("%.3f",1)        ,'#      ',sprintf("%.3f",temp0n[2]),'        (',
#        sprintf("%.3f",temp0n[3]),',',sprintf("%.3f",temp0n[4]),')\n')
#  }
#  if(temp0n[1]<=1)
#  {
#    cat('    C0')
#    cat(N,'(Sorensen)',sprintf("%.3f",temp0n[1]),'       ',sprintf("%.3f",temp0n[2]),'        (',
#        sprintf("%.3f",temp0n[3]),',',sprintf("%.3f",temp0n[4]),')\n')
#  }
#  #temp1n=x$overlap[2,]
#  #cat('    C1')
#  #cat(N,'          ',sprintf("%.3f",temp1n[1]),'       ',sprintf("%.3f",temp1n[2]),'        (',
#  #    sprintf("%.3f",temp1n[3]),',',sprintf("%.3f",temp1n[4]),')\n')
#  temp1n=x$overlap[2,]
#  if(temp1n[1]>1)
#  {
#    cat('    C1')
#    cat(N)
#    cat('*(Horn)    ',sprintf("%.3f",1)        ,'#      ',sprintf("%.3f",temp1n[2]),'        (',
#        sprintf("%.3f",temp1n[3]),',',sprintf("%.3f",temp1n[4]),')\n')
#  }
#  if(temp1n[1]<=1)
#  {
#    cat('    C1')
#    cat(N)
#    cat('*(Horn)    ',sprintf("%.3f",temp1n[1]),'       ',sprintf("%.3f",temp1n[2]),'        (',
#        sprintf("%.3f",temp1n[3]),',',sprintf("%.3f",temp1n[4]),')\n')
    
#  }
  
#  #cat('*          ',sprintf("%.3f",1-temp1n[1]),'       ',sprintf("%.3f",temp1n[2]),'        (',
#  #     sprintf("%.3f",max(1-temp1n[1]-1.96*temp1n[2],0)),',',sprintf("%.3f",min(1-temp1n[1]+1.96*temp1n[2],1)),')\n')
#  temp2n=x$overlap[3,]
#  cat('    C2')
#  cat(N,'(Morisita)',sprintf("%.3f",temp2n[1]),'       ',sprintf("%.3f",temp2n[2]),'        (',
#      sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')  
  
#  if(N==3)
#  {
#    temp33=x$overlap[4,]
#    cat('    C33','          ',sprintf("%.3f",temp33[1]),'       ',sprintf("%.3f",temp33[2]),'        (',
#        sprintf("%.3f",temp33[3]),',',sprintf("%.3f",temp33[4]),')\n')
#  }
#  cat('\n')
#  cat('    C0')
#  cat(N,': A similarity measure of comparing',N,
#      'communities using empirical method.\n')
#  #cat('    C1')
#  #cat(N,': A similarity measure of comparing',N,
#  #   'communities based on equal sample size among all communities.\n')
#  cat('    C1')
#  cat(N)
#  cat('*')
#  cat(': A similarity measure of comparing',N,
#      'communities based on equal-effort sample size among all communities.\n')
#  cat('    C2')
#  cat(N,': A similarity measure of comparing',N,'communities based on shared information between any two communities.\n')
#  if(N==3)
#  {
#    cat('    C33',': A similarity measure of comparing 3 communities using all shared information.\n')
#  }
#  cat('    
#      Confidence Interval: Based on an improved bootstrap percentile method. (recommend for use in the case when 
#      similarity is close to 0 or 1 ) \n\n')
#  cat('    # if the estimate is greater than 1, it is replaced by 1.\n\n')
#  cat('    Pairwise Comparison:\n\n')
#  cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
#  Cqn_PC <- x$pairwise
#  no.temp=1
#  for(i in 1:(N-1))
#  {
#    for(j in (i+1):N)
#    {
#      temp=Cqn_PC[no.temp,]
#      if(q==0){cat('    C02(')}
#     if(q==1 & x$method=="relative"){cat('    C12(')}
#      if(q==1 & x$method=="absolute"){cat('    C12*(')}
#      if(q==2){cat('    C22(')}
#      cat(i)
#      cat(',')
#      cat(j)
#      if(temp[1]>1)
#      {cat(')','     ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
#           sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
#      }
#      if(temp[1]<=1)
#      {cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
#           sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
#      no.temp=no.temp+1
#    }
#  }
  
#  cat('\n')
#  cat('    Average Pairwise =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n')
#  if(q==0 || q==1)
#  {
#    cat('
#        If the lower bound is less than 0, it is replaced by 0; if the upper bound
#        is greater than 1, it is replaced by 1.\n\n')
#  }
#  if(q==2)
#  {
#    cat('
#        MLE is used for replacing nearly unbiased estimate because the estimate 
#        is greater than 1.
#        If the lower bound is less than 0, it is replaced by 0; if the upper bound
#       is greater than 1, it is replaced by 1.\n\n')
#  }
#  cat('    Similarity Matrix: \n\n')
#  C_SM=x$similarity.matrix
#  
#  if(q==0){cat('    C02(i,j)   \t')}
#  if(q==1 & x$method=="relative"){cat('    C12(i,j)   \t')}
#  if(q==1 & x$method=="absolute"){cat('    C12*(i,j)   ')}
#  if(q==2){cat('    C22(i,j)   \t')}
#  for(i in 1:N)
#  {
#    cat(i,'\t')
#  }
#  cat('\n')
#  for(i in 1:N)
#  {
#    cat('       ',i,'\t')
#    for(j in 1:N)
#    {
#      if(i>j){cat('\t')}
#      if(i<=j) {
#        if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),' \t')
#        else cat(sprintf("%.3f#",1),'\t')
#      }
#    }
#    cat('\n')
#  }
#  if(q==2)
#  {
#    cat('\n')
#    cat('(3) ESTIMATION OF MORISITA DISSIMILARITY IN',N,'COMMUNITIES\n\n')
#    cat('    Estimator','     Estimate','     s.e.','    95% Confidence Interval\n\n')
#    cat('    1 - C2')
#    cat(N,'      ',sprintf("%.3f",1-temp2n[1]),'       ',sprintf("%.3f",temp2n[2]),'        (',
#        sprintf("%.3f",max(0,1-temp2n[1]-1.96*temp2n[2])),',',sprintf("%.3f",min(1,1-temp2n[1]+1.96*temp2n[2])),')\n')
#    if(N==3)
#    {
#      cat('    1 - C33','      ',sprintf("%.3f",1-temp33[1]),'       ',sprintf("%.3f",temp33[2]),'        (',
#          sprintf("%.3f",max(0,1-temp33[1]-1.96*temp33[2])),',',sprintf("%.3f",min(1,1-temp33[1]+1.96*temp33[2])),')\n')
#    }
#    cat('\n')
#    cat('    1 - C2')
#    cat(N,': This is the genetic diversity measure D defined in Jost (2008) for
#        comparing',N,'communities.\n')
#    if(N==3)
#    {
#      cat('    1 - C33',': A genetic diversity measure for comparing 3 subpopulations based on
#          all shared information.\n')
#    }
#    cat('\n')
#    cat('    Pairwise Comparison:\n\n')
#    cat('    Estimator','          Estimate','     s.e.','       95% Confidence Interval\n\n')
#    no.temp=1
#    for(i in 1:(N-1))
#    {
#      for(j in (i+1):N)
#      {
#        temp=Cqn_PC[no.temp,]
#        cat('    1-C22(')
#        cat(i)
#        cat(',')
#        cat(j)
#        cat(')','        ',sprintf("%.3f",1-temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
#            sprintf("%.3f",temp[5]),',',sprintf("%.3f",temp[6]),')\n')
#        no.temp=no.temp+1
#      }
#    }
#    cat('\n')
#    cat('    Average Pairwise =',sprintf("%.3f",1-mean(Cqn_PC[,1])),'\n')
#    cat('
#        MLE is used for replacing nearly unbiased estimate because the estimate 
#        is greater than 1.
#        If the lower bound is less than 0, it is replaced by 0; if the upper bound
#        is greater than 1, it is replaced by 1.\n\n')
#    cat('    1-C22: This is the genetic diversity measure D defined in Jost (2008) for
#        comparing 2 subpopulations.')
#    cat('\n\n')
#    cat('    Dissimilarity Matrix: \n\n')
#    cat('    1-C22(i,j)\t')
#    for(i in 1:N)
#    {
#      cat(i,'\t')
#    }
#    cat('\n')
#    for(i in 1:N)
#    {
#      cat('       ',i,'\t')
#      for(j in 1:N)
#      {
#        if(i>j){cat('\t')}
#        if(i<=j){
#          if(C_SM[i,j]<=1) cat(sprintf("%.3f",1-C_SM[i,j]),' \t')
#          else cat(sprintf("%.3f#",1),' \t')
#        }
#      }
#      cat('\n')
#    }
#    cat('\n')
#    }
#  cat('
#      References:
      
#      Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-
#      stage probabilistic approach to multiple-community similarity indices. 
#      Biometrics, 64, 1178-1186.
#      
#      Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular
#      Ecology, 17, 4015-4026.
#      ')
#  cat('\n')
#  }

	

