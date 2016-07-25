GST_equ=function(X,method=c("est","mle"))
{ 
  if(method=="est"){
    X=as.matrix(X)
    no.community=length(X[1,])
    n=colSums(X)
    temp1=X%*%matrix(c(1/n),no.community,1)
    temp2=sum(sapply(1:no.community,function(k) sum((X[,k]/n[k])^2) ))
    
    Hs_hat=1-1/no.community*sum(sapply(1:no.community,function(k)  sum(X[,k]*(X[,k]-1)/n[k]/(n[k]-1))))
    Ht_hat=1-1/no.community^2*sum(sapply(1:no.community,function(k)  sum(X[,k]*(X[,k]-1)/n[k]/(n[k]-1))))-1/no.community^2*(sum(temp1^2)-temp2) 
    GST=1-Hs_hat/Ht_hat
  }else{
    X = as.matrix(X)
    N <- length(X[1,])
    n <- colSums(X)
    ps <- sapply(1:N, FUN = function(i){
      (X[,i] /n[i] )^2
    })
    Hs <- 1 - sum(ps)/N
    Ht <- 1 - sum((apply(sapply(1:N, FUN = function(i){
      (X[,i] /n[i] )
    }),MARGIN = 1, FUN = sum)/N)^2)
    GST=1-Hs/Ht
  }
  return(GST)
} 
GST_se_equ <- function(X, nboot=50){
  nboot=50
  plus_CI <-function(x){
    if(x[1] >= 1) x[1] <- 1
    c(x, max(0,x[1]-1.96*x[2]), min(1,x[1]+1.96*x[2]))
  }
  boot.gst.mle=rep(0,nboot)
  boot.gst.est=rep(0,nboot)
  gst.mle=GST_equ(X, method="mle")
  gst.est=GST_equ(X, method="est")
  for(i in 1:nboot)
  {
    p <- Boots.pop(X)
    boot.X=sapply(1:dim(X)[2],function(j)   rmultinom(1,sum(X[,j]),p[,j] ) )
    boot.gst.mle[i]=GST_equ(boot.X, method="mle")
    boot.gst.est[i]=GST_equ(boot.X, method="est")
  }
  se_mle=sd(boot.gst.mle)
  se_est=sd(boot.gst.est)
  out1= plus_CI(c(gst.mle, se_mle))
  out2= plus_CI(c(gst.est, se_est))
  out <- rbind(out1, out2)
  return(out)
}
print.spadeGenetic <- function(x, ...)
{
  cat('\n(1) BASIC DATA INFORMATION:\n\n')
  cat('    The loaded set includes abundance (or frequency) data from',x$info[1],'subpopulations\n')
  cat('    and a total of',x$info[2],'distinct alleles are found.\n\n')
  cat('    Sample size in each subpopulation                           n1   =', x$info[3],'\n')
  N <- x$info[1]
  q <- x$q
  for(j in 2:N){
    cat('                                                               ','n')
    cat(j,'  =',x$info[2+j],'\n')   
  }
  cat('\n')
  cat('    Number of observed alleles in one subpopulation             D1   =', x$info[N+3],'\n')
  
  for(j in 2:N){
    cat('                                                               ','D')
    cat(j,'  =',x$info[N+2+j],'\n')
  }
  cat('\n')
  cat('    Number of observed shared alleles in two subpopulations     D12  =', x$info[3+2*N], '\n')
  
  if(N>2){
    k <- 1
    for(i in 1:(N-1)){     
      for(j in (i+1):N){
        if(i==1 & j==2) next
        cat('                                                               ','D')
        cat(i,j,'  = ', x$info[3+2*N+k], '\n', sep="")
        k <- k + 1
      }
    }
  }
  cat('\n')
  if(N==3)
  {
    cat('    Number of shared alleles in three subpopulations            D123 =',rev(x$info)[2],'\n\n') 
  }
  cat('    Number of bootstrap replications for s.e. estimate                ',rev(x$info)[1],'\n\n')
  cat('(2) EMPIRICAL DIS-SIMILARITY INDICES: \n\n')
  cat('                                         Estimate       s.e.       95%Lower     95%Upper\n')
  cat('    (a) Classic richness-based similarity\n\n')
  temp <- apply(as.matrix(x$Empirical_incidence), 2, as.numeric)
  cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')
  cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')
  cat('    (b) Measures for comparing species relative abundances\n\n')
  temp <- apply(as.matrix(x$Empirical_ew), 2, as.numeric)
  cat('        1-C1');cat(N);cat('=1-U1');cat(N);cat(' (q=1, Horn)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('        1-C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
  cat('        1-U2');cat(N);cat(' (q=2, Regional overlap)    ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
  cat('        Gst                              ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
  cat('    (c) Measures for comparing size-weighted species relative abundances\n\n')
  temp <- x$Empirical_ee
  cat('        Horn size-weighted(q=1)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('(3) ESTIMATED DIS-SIMILARITY INDICES: \n\n')
  cat('                                         Estimate       s.e.       95%Lower     95%Upper\n')
  cat('    (a) Classic richness-based similarity\n\n')
  temp <- apply(as.matrix(x$estimated_incidence), 2, as.numeric)
  if(temp[1,1]>1) {cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",1),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
  if(temp[1,1]<=1){cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
  if(temp[2,1]>1) {cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",1) ,'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
  if(temp[2,1]<=1){cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
  cat('    (b) Measures for comparing species relative abundances\n\n')
  temp <- apply(as.matrix(x$estimated_ew), 2, as.numeric)
  cat('        1-C1');cat(N);cat('=1-U1');cat(N);cat(' (q=1, Horn)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('        1-C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
  cat('        1-U2');cat(N);cat(' (q=2, Regional overlap)    ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
  cat('        Gst                              ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
  cat('    (c) Measures for comparing size-weighted species relative abundances\n\n')
  temp <- x$estimated_ee
  cat('        Horn size-weighted (q=1)         ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('(4) ESTIMATED PAIRWISE DIS-SIMILARITY:\n\n')
  cat('    ---------------------Measure 1-C');cat(q);cat('2----------------------\n\n')
  cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
  ###################################################CqN_ Equal weight
  Cqn_PC <- x$pairwise$CqN_ew
  no.temp=1
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      temp=Cqn_PC[no.temp,]
      if(q==0){cat('    1-C02(')}
      if(q==1){cat('    1-C12(')}
      #if(q==1 & x$method=="relative"){cat('    C12(')}
      #f(q==1 & x$method=="absolute"){cat('    C12*(')}
      if(q==2){cat('    1-C22(')}
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
  cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
  cat('    Pairwise dis-similarity matrix: \n\n')
  C_SM=x$similarity.matrix$CqN_ew
  if(q==0){cat('    1-C02(i,j)')}
  if(q==1){cat('    1-C12(i,j)')}
  if(q==2){cat('    1-C22(i,j)')}
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
  if(q==1) cat('    ---------------Measure Horn size-weighted-----------------\n\n')
  if(q!=1){cat('    ---------------------Measure 1-U');cat(q);cat('2----------------------\n\n')}
  cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
  if(q!=1) Cqn_PC <- x$pairwise$UqN_ew
  if(q==1) Cqn_PC <- x$pairwise$CqN_ee
  no.temp=1
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      temp=Cqn_PC[no.temp,]
      if(q==0){cat('    1-U02(')}
      if(q==1){cat('    Horn(')}
      #if(q==1 & x$method=="relative"){cat('    C12(')}
      #f(q==1 & x$method=="absolute"){cat('    C12*(')}
      if(q==2){cat('    1-U22(')}
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
  cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
  cat('    Pairwise dis-similarity matrix: \n\n')
  if(q!=1) C_SM <- x$similarity.matrix$UqN_ew
  if(q==1) C_SM <- x$similarity.matrix$CqN_ee
  
  if(q==0){cat('    1-U02(i,j)')}
  if(q==1){cat('    Horn(i,j) ')}
  if(q==2){cat('    1-U22(i,j)')}
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
  cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
  cat('\n')
  
  
  #cat('(2) NEARLY UNBIASED ESTIMATION OF ALLELIC DIFFERENTIATION OR MORISITA DISSIMILARITY IN ',N,'SUBPOPULATIONS:\n\n')
  #cat('    Estimator','     Estimate','      s.e.','      95% Confidence Interval\n\n')
  #if ( q==0 ) temp0n=x$overlap[1,]
  #if ( q==1 ) temp0n=x$overlap[3,]
  #if ( q==2 ) temp0n=x$overlap[4,]
  #if ( temp0n[1]<=1 ){
  #  cat('    1-C')
  #  cat(q)
  #  cat(N,sprintf("         %.3f",1-temp0n[1]),'       ',sprintf("%.3f",temp0n[2]),'        (',
  #      sprintf("%.3f",1-temp0n[4]),',',sprintf("%.3f",1-temp0n[3]),')\n')
  #}
  #if ( temp0n[1]>1 ){
  #  cat('    1-C')
  #  cat(q)
  #  cat(N,sprintf("         %.3f#",0),'      ',sprintf("%.3f",temp0n[2]),'        (',
  #      sprintf("%.3f",1-temp0n[4]),',',sprintf("%.3f",1-temp0n[3]),')\n')
  #}
  # 
  #cat('\n')
  #cat('    1-C')
  #cat(q)
  #cat(N,':This is the genetic diversity measure defined in Jost (2008) for comparing subpopulations based 
  #    on allele shared information between any two subpopulations.\n')
  #cat('    
  #    Confidence Interval: Based on an improved bootstrap percentile method. (recommend for use in the case when 
  #    similarity is close to 0 or 1 ) \n\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('    Pairwise Comparison:\n\n')
  #cat('    Estimator','     Estimate','          s.e.','      95% Confidence Interval\n\n')
  #Cqn_PC <- x$pairwise  
  #no.temp=1
  #share_index = sapply(1:choose(N,2), function(i){x$info[2+2*N+i]})
  #no.temp2=1
  #for(i in 1:(N-1))
  #{
  #  for(j in (i+1):N)
  #  {
  #    temp=Cqn_PC[no.temp,]
  #    cat('    1-C')
  #    cat(q)
  #    cat('2(')
  #    cat(i)
  #    cat(',')
  #    cat(j)
  #    if ( share_index[no.temp2]!=0 )
  #    {
  #      if ( temp[1]>1 ){
  #        cat(')','        ',sprintf("%.3f #",0),'     ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",1-temp[4]),',',sprintf("%.3f",1-temp[3]),')\n')
  #      }
  #      if ( temp[1]<=1 ){
  #        cat(')','        ',sprintf("%.3f",1-temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",1-temp[4]),',',sprintf("%.3f",1-temp[3]),')\n')
  #      }
  #    }
  #    else
  #    {
  #      cat(')','        ',sprintf("%.3f",1),'       ',sprintf("%.3f",0),'        (',
  #          sprintf("%.3f",1),',',sprintf("%.3f",1),')\n')
  #      #         cat(')','        ',sprintf("%.3f##",1),'     ',sprintf("%.3f##",0),'      (',
  #      #             sprintf("%.3f",1),',',sprintf("%.3f",1),')##\n')
  #    }
  #    no.temp=no.temp+1
  #    no.temp2=no.temp2+1
  #  }
  #}
  #cat('\n')
  #Cqn_PC <- x$pairwise
  #cat('    Average Pairwise =',sprintf("%.3f",1-mean(Cqn_PC[,1])),'\n')
  #cat('    ## There are no shared species, thus estimated similarity is zero and should be used for caution.\n\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('    1-C')
  #cat(q)
  #cat('2: This is the genetic diversity measure D defined in Jost (2008) for comparing 2 subpopulations.')
  #cat('\n\n')
  #cat('    Dissimilarity Matrix: \n\n')
  #C_SM=x$similarity.matrix
  #cat('    1-C')
  #cat(q)
  #cat('2(i,j)\t')
  #for(i in 1:N)
  #{
  #  cat(i,'\t')
  #}
  #cat('\n')
  #for(i in 1:N)
  #{
  #  cat('       ',i,'\t')
  #  for(j in 1:N)
  #  {
  #    if(i>j){cat('\t')}
  #    if(i<=j){
  #      if (1-C_SM[i,j]>=0) cat(sprintf("%.3f",abs(1-C_SM[i,j])),'\t')
  #      else cat(sprintf("%.3f#",0),'\t')
  #    }
  #  }
  #  cat('\n')
  #}
  #cat('\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #  
  #cat('(3)  NEARLY UNBIASED ESTIMATION OF MORISITA SIMILARITY IN ',N,'SUBPOPULATIONS:\n\n')
  #cat('    Estimator','     Estimate','      s.e.','      95% Confidence Interval\n\n')
  #  
  #if ( q==0 ) temp2n=x$overlap[1,]
  #if ( q==1 ) temp2n=x$overlap[3,]
  #if ( q==2 ) temp2n=x$overlap[4,]
  #if ( temp2n[1]>=1 ){
  #  cat('    C')
  #  cat(q)
  #  cat(N,sprintf("           %.3f#",1),'      ',sprintf("%.3f",temp2n[2]),'       (',
  #      sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')
  #}
  #if ( temp2n[1]<1 ){
  #  cat('    C')
  #  cat(q)
  #  cat(N,sprintf("           %.3f",temp2n[1]),'       ',sprintf("%.3f",temp2n[2]),'       (',
  #      sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')
  #}
  #cat('\n')
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('\n')
  #cat('    C')
  #cat(q)
  #cat(N,': A similarity measure of comparing 3 subpopulations based on allele shared information between any two subpopulations.\n')
  #cat('\n')
  #cat('    Pairwise Comparison:\n\n')
  #cat('    Estimator','     Estimate','      s.e.','      95% Confidence Interval\n\n')
  #Cqn_PC <- x$pairwise
  #no.temp=1
  #share_index = sapply(1:choose(N,2), function(i){x$info[2+2*N+i]})
  #no.temp2=1
  #for(i in 1:(N-1))
  #{
  #  for(j in (i+1):N)
  #  {
  #    temp=Cqn_PC[no.temp,]
  #    if(q==0){cat('    C02(')}
  #    if(q==1){cat('    C12(')}
  #    if(q==2){cat('    C22(')}
  #    cat(i)
  #    cat(',')
  #    cat(j)
  #    if ( share_index[no.temp2]!=0 ){
  #      if ( temp[1]<=1 ){
  #        cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
  #      }
  #      else {
  #        cat(')','     ',sprintf("%.3f #", 1),'     ',sprintf("%.3f",temp[2]),'        (',
  #            sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
  #      }
  #    }
  #    else
  #    {
  #      cat(')','     ',sprintf("%.3f",0),'       ',sprintf("%.3f",0),'        (',
  #          sprintf("%.3f",0),',',sprintf("%.3f",0),')\n')
  #    }
  #    # Cqn_PC[i,1]=0
  #    no.temp=no.temp+1
  #    no.temp2=no.temp2+1
  #  }
  #}
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('\n')
  #Cqn_PC <- x$pairwise
  #cat('    Average Pairwise =',sprintf("%.3f",mean(Cqn_PC[1:choose(N,2),1])),'\n')
  #cat('\n')
  #cat('    ## There are no shared species, thus estimated similarity is zero and should be used for caution.\n\n')
  #
  #cat('    Similarity Matrix: \n\n')
  #C_SM=x$similarity.matrix
  #if(q==0){cat('    C02(i,j)   \t')}
  #if(q==1){cat('    C12(i,j)   \t')}
  #if(q==2){cat('    C22(i,j)   \t')}
  #for(i in 1:N)
  #{
  #  cat(i,'\t')
  #}
  #cat('\n')
  #for(i in 1:N)
  #{
  #  cat('       ',i,'\t')
  #  for(j in 1:N)
  #  {
  #    if(i>j){cat('\t')}
  #    if(i<=j){
  #      if (C_SM[i,j]<=1) cat(sprintf("%.3f",abs(C_SM[i,j])),'\t')
  #      else cat(sprintf("%.3f#",1),'\t')
  #    }
  #  }
  #  cat('\n')
  #}
  #cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
  #    is greater than 1, replace it by 1.\n\n')
  #cat('\n')
  #cat('    References:\n')
  #cat('    Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-
  #    stage probabilistic approach to multiple-community similarity indices. 
  #    Biometrics, 64, 1178-1186.\n')
  #cat('    Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular
  #    Ecology, 17, 4015-4026.')
  
}