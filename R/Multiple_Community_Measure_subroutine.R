C0n_equ=function(X)
{
  X=as.matrix(X)
  I=which(X>0)
  X[I]=rep(1,length(I))
  no.community=length(X[1,])
  C0n_num=sum(rowSums(X)>0)
  C0n_dem=sum(X)
  C0n=no.community/(1-no.community)*( C0n_num-C0n_dem)/C0n_dem
  C0n
}
correct_obspi<- function(X)
{
   Sobs <- sum(X > 0) 	
   n <- sum(X)		  	
   f1 <- sum(X == 1) 	
   f2 <- sum(X == 2) 	
   a <- ifelse(f1 == 0, 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n)
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
C1n_equ=function(method=c("absolute"),X,boot=200)
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
   c2n
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
   C33_num/C33_dem                      
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
	cat('\n(1) BASIC DATA INFORMATION:\n\n')
	cat('    The loaded set includes abundance (or frequency) data from',x$info[1],'communities\n')
	cat('    and a total of',x$info[2],'distinct species.\n\n')
		cat('   (Number of observed individuals in each community)       n1 =', x$info[3],'\n')
	N <- x$info[1]
  q <- x$q
	for(j in 2:N){
		cat('                                                           ','n')
		cat(j,'=',x$info[2+j],'\n')   
    }
		cat('\n')
		cat('   (Number of observed species in one community)            D1 =', x$info[N+3],'\n')
   
   for(j in 2:N){
		cat('                                                           ','D')
		cat(j,'=',x$info[N+2+j],'\n')
    }
		cat('\n')
		cat('   (Number of observed shared species in two communities)  ','D12 =', x$info[3+2*N], '\n')
   
   if(N>2){
    k <- 1
    for(i in 1:(N-1)){     
          for(j in (i+1):N){
		  if(i==1 & j==2) next
        cat('                                                           ','D')
        cat(i,j,' = ', x$info[3+2*N+k], '\n', sep="")
		k <- k + 1
          }
      }
    }
   cat('\n')
   if(N==3)
   {
      cat('   (Number of observed shared species in three communities)','D123','=',rev(x$info)[2],'\n\n') 
   }
   cat('   (Bootstrap replications for s.e. estimate)              ',rev(x$info)[1],'\n\n')
   cat('(2) ESTIMATION OF OVERLAP MEASURE IN',N,'COMMUNITIES:\n\n')
   cat('    Estimator','     Estimate','     s.e.','     95% Confidence Interval\n\n')
   temp0n=x$overlap[1,]
   cat('    C0')
   cat(N,'(Sorensen)',sprintf("%.3f",temp0n[1]),'       ',sprintf("%.3f",temp0n[2]),'        (',
       sprintf("%.3f",temp0n[3]),',',sprintf("%.3f",temp0n[4]),')\n')
   
   #temp1n=x$overlap[2,]
   #cat('    C1')
   #cat(N,'          ',sprintf("%.3f",temp1n[1]),'       ',sprintf("%.3f",temp1n[2]),'        (',
   #    sprintf("%.3f",temp1n[3]),',',sprintf("%.3f",temp1n[4]),')\n')
   temp1n=x$overlap[2,]
   cat('    C1')
   cat(N)
   cat('*(Horn)    ',sprintf("%.3f",temp1n[1]),'       ',sprintf("%.3f",temp1n[2]),'        (',
      sprintf("%.3f",temp1n[3]),',',sprintf("%.3f",temp1n[4]),')\n')
	 #cat('*          ',sprintf("%.3f",1-temp1n[1]),'       ',sprintf("%.3f",temp1n[2]),'        (',
	 #     sprintf("%.3f",max(1-temp1n[1]-1.96*temp1n[2],0)),',',sprintf("%.3f",min(1-temp1n[1]+1.96*temp1n[2],1)),')\n')
   temp2n=x$overlap[3,]
   cat('    C2')
   cat(N,'(Morisita)',sprintf("%.3f",temp2n[1]),'       ',sprintf("%.3f",temp2n[2]),'        (',
       sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')
   if(N==3)
   {
      temp33=x$overlap[4,]
      cat('    C33','          ',sprintf("%.3f",temp33[1]),'       ',sprintf("%.3f",temp33[2]),'        (',
      sprintf("%.3f",temp33[3]),',',sprintf("%.3f",temp33[4]),')\n')
   }
   cat('\n')
   cat('    C0')
   cat(N,': A similarity measure of comparing',N,
      'communities using empirical method.\n')
   #cat('    C1')
   #cat(N,': A similarity measure of comparing',N,
   #   'communities based on equal sample size among all communities.\n')
   cat('    C1')
   cat(N)
   cat('*')
   cat(': A similarity measure of comparing',N,
      'communities based on equal-effort sample size among all communities.\n')
   cat('    C2')
   cat(N,': A similarity measure of comparing',N,'communities based on shared information between any two communities.\n')
   if(N==3)
   {
      cat('    C33',': A similarity measure of comparing 3 communities using all shared information.\n')
   }
   cat('    
    Confidence Interval: Based on an improved bootstrap percentile method. (recommend for use in the case when 
                         similarity is close to 0 or 1 ) \n\n')
   cat('    Pairwise Comparison:\n\n')
   cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
   Cqn_PC <- x$pairwise
   no.temp=1
   for(i in 1:(N-1))
   {
       for(j in (i+1):N)
       {
           temp=Cqn_PC[no.temp,]
           if(q==0){cat('    C02(')}
           if(q==1 & x$method=="relative"){cat('    C12(')}
           if(q==1 & x$method=="absolute"){cat('    C12*(')}
           if(q==2){cat('    C22(')}
           cat(i)
           cat(',')
           cat(j)
           cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
           sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
           no.temp=no.temp+1
       }
   }
   
   cat('\n')
   cat('    Average Pairwise =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n')
   if(q==0 || q==1)
   {
    cat('
    If the lower bound is less than 0, it is replaced by 0; if the upper bound
       is greater than 1, it is replaced by 1.\n\n')
   }
   if(q==2)
   {
    cat('
    MLE is used for replacing nearly unbiased estimate because the estimate 
	        is greater than 1.
    If the lower bound is less than 0, it is replaced by 0; if the upper bound
       is greater than 1, it is replaced by 1.\n\n')
   }
   cat('    Similarity Matrix: \n\n')
   C_SM=x$similarity.matrix
   
   if(q==0){cat('    C02(i,j)   \t')}
   if(q==1 & x$method=="relative"){cat('    C12(i,j)   \t')}
   if(q==1 & x$method=="absolute"){cat('    C12*(i,j)   ')}
   if(q==2){cat('    C22(i,j)   \t')}
   for(i in 1:N)
   {
       cat(i,'\t')
   }
   cat('\n')
   for(i in 1:N)
   {
       cat('       ',i,'\t')
       for(j in 1:N)
       {
           if(i>j){cat('\t')}
           if(i<=j)cat(sprintf("%.3f",C_SM[i,j]),'\t')
       }
       cat('\n')
   }
   if(q==2)
   {
   cat('\n')
   cat('(3) ESTIMATION OF MORISITA DISSIMILARITY IN',N,'COMMUNITIES\n\n')
   cat('    Estimator','     Estimate','     s.e.','    95% Confidence Interval\n\n')
   cat('    1 - C2')
   cat(N,'      ',sprintf("%.3f",1-temp2n[1]),'       ',sprintf("%.3f",temp2n[2]),'        (',
       sprintf("%.3f",max(0,1-temp2n[1]-1.96*temp2n[2])),',',sprintf("%.3f",min(1,1-temp2n[1]+1.96*temp2n[2])),')\n')
   if(N==3)
   {
      cat('    1 - C33','      ',sprintf("%.3f",1-temp33[1]),'       ',sprintf("%.3f",temp33[2]),'        (',
      sprintf("%.3f",max(0,1-temp33[1]-1.96*temp33[2])),',',sprintf("%.3f",min(1,1-temp33[1]+1.96*temp33[2])),')\n')
   }
   cat('\n')
   cat('    1 - C2')
   cat(N,': This is the genetic diversity measure D defined in Jost (2008) for
              comparing',N,'communities.\n')
   if(N==3)
   {
      cat('    1 - C33',': A genetic diversity measure for comparing 3 subpopulations based on
              all shared information.\n')
   }
   cat('\n')
   cat('    Pairwise Comparison:\n\n')
   cat('    Estimator','          Estimate','     s.e.','       95% Confidence Interval\n\n')
   no.temp=1
   for(i in 1:(N-1))
   {
       for(j in (i+1):N)
       {
           temp=Cqn_PC[no.temp,]
           cat('    1-C22(')
           cat(i)
           cat(',')
           cat(j)
           cat(')','        ',sprintf("%.3f",1-temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
               sprintf("%.3f",temp[5]),',',sprintf("%.3f",temp[6]),')\n')
           no.temp=no.temp+1
       }
   }
   cat('\n')
   cat('    Average Pairwise =',sprintf("%.3f",1-mean(Cqn_PC[,1])),'\n')
   cat('
    MLE is used for replacing nearly unbiased estimate because the estimate 
	        is greater than 1.
    If the lower bound is less than 0, it is replaced by 0; if the upper bound
       is greater than 1, it is replaced by 1.\n\n')
   cat('    1-C22: This is the genetic diversity measure D defined in Jost (2008) for
           comparing 2 subpopulations.')
   cat('\n\n')
   cat('    Dissimilarity Matrix: \n\n')
   cat('    1-C22(i,j)\t')
   for(i in 1:N)
   {
       cat(i,'\t')
   }
   cat('\n')
   for(i in 1:N)
   {
       cat('       ',i,'\t')
       for(j in 1:N)
       {
           if(i>j){cat('\t')}
           if(i<=j)cat(sprintf("%.3f",1-C_SM[i,j]),'\t')
       }
       cat('\n')
   }
   cat('\n')
   }
   cat('
   References:
   
   Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-
   stage probabilistic approach to multiple-community similarity indices. 
   Biometrics, 64, 1178-1186.
   
   Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular
   Ecology, 17, 4015-4026.
   ')
   cat('\n')
}
