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
GST_se_equ <- function(X, nboot=50)
{
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
  cat('    (a) Classic richness-based dis-similarity\n\n')
  temp <- apply(as.matrix(x$Empirical_richness), 2, as.numeric)
  cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')
  cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')
  cat('    (b) Measures for comparing alleles relative abundances\n\n')
  temp <- apply(as.matrix(x$Empirical_relative), 2, as.numeric)
  cat('        1-C1');cat(N);cat('=1-U1');cat(N);cat(' (q=1, Horn)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('        1-C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
  cat('        1-U2');cat(N);cat(' (q=2, Regional diff.)      ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
  cat('        Gst                              ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
  cat('    (c) Measures for comparing size-weighted alleles relative abundances\n\n')
  temp <- x$Empirical_WtRelative
  cat('        Horn size-weighted(q=1)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('(3) ESTIMATED DIS-SIMILARITY INDICES: \n\n')
  cat('                                         Estimate       s.e.       95%Lower     95%Upper\n')
  cat('    (a) Classic richness-based dis-similarity\n\n')
  temp <- apply(as.matrix(x$estimated_richness), 2, as.numeric)
  if(temp[1,1]>1) {cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",1),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
  if(temp[1,1]<=1){cat('        1-C0');cat(N);cat(' (q=0, Sorensen)            ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n')}
  if(temp[2,1]>1) {cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",1) ,'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
  if(temp[2,1]<=1){cat('        1-U0');cat(N);cat(' (q=0, Jaccard)             ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n\n')}
  cat('    (b) Measures for comparing alleles relative abundances\n\n')
  temp <- apply(as.matrix(x$estimated_relative), 2, as.numeric)
  cat('        1-C1');cat(N);cat('=1-U1');cat(N);cat(' (q=1, Horn)          ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('        1-C2');cat(N);cat(' (q=2, Morisita-Horn)       ',sprintf("%.4f",temp[2,1]),'     ',sprintf("%.4f",temp[2,2]),'     ',sprintf("%.4f",temp[2,3]),'     ',sprintf("%.4f",temp[2,4]),'\n')
  cat('        1-U2');cat(N);cat(' (q=2, Regional diff.)      ',sprintf("%.4f",temp[3,1]),'     ',sprintf("%.4f",temp[3,2]),'     ',sprintf("%.4f",temp[3,3]),'     ',sprintf("%.4f",temp[3,4]),'\n\n')
  cat('        Gst                              ',sprintf("%.4f",temp[4,1]),'     ',sprintf("%.4f",temp[4,2]),'     ',sprintf("%.4f",temp[4,3]),'     ',sprintf("%.4f",temp[4,4]),'\n\n')
  cat('    (c) Measures for comparing size-weighted alleles relative abundances\n\n')
  temp <- x$estimated_WtRelative
  cat('        Horn size-weighted (q=1)         ',sprintf("%.4f",temp[1,1]),'     ',sprintf("%.4f",temp[1,2]),'     ',sprintf("%.4f",temp[1,3]),'     ',sprintf("%.4f",temp[1,4]),'\n\n')
  cat('(4) ESTIMATED PAIRWISE DIS-SIMILARITY:\n\n')
  if(q == 0){
    cat('    -----------------------Measure 1-C02------------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    ###################################################CqN_ Equal weight
    Cqn_PC <- x$pairwise$C02
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    1-C02(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','   ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','   ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise dis-similarity matrix: \n\n')
    C_SM=x$dissimilarity_matrix$C02
    cat('    1-C02(i,j)')
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
        if(i==j){cat(round(0,0),'  ')}
        if(i<j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    cat('    -----------------------Measure 1-U02------------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    Cqn_PC <- x$pairwise$U02
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    1-U02(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','   ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','   ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise dis-similarity matrix: \n\n')
    C_SM <- x$dissimilarity_matrix$U02
    
    cat('    1-U02(i,j)')
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
        if(i==j){cat(round(0,0),'  ')}
        if(i<j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
    cat('\n')
  }
  if( q == 1 ){
    cat('    --------------------Measure 1-C12 (=1-U12)----------------------\n\n')
    cat('    Estimator','     Estimate','     s.e.','       95% Confidence Interval\n\n')
    ###################################################CqN_ Equal weight
    Cqn_PC <- x$pairwise$C12
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    1-C12(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','   ',sprintf("%.3f",1)      ,'#      ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','   ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    
    cat('\n')
    cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise dis-similarity matrix: \n\n')
    C_SM=x$dissimilarity_matrix$C12
    cat('    1-C12(i,j)')
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
        if(i==j){cat(round(0,0),'  ')}
        if(i<j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    ###################################################UqN_ Equal weight
    cat('    ------------------Measure Horn size-weighted--------------------\n\n')
    cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
    Cqn_PC <- x$pairwise$Horn
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
        {cat(')','    ',sprintf("%.3f",1)      ,'#       ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','    ',sprintf("%.3f",temp[1]),'        ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    
    cat('\n')
    cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise dis-similarity matrix: \n\n')
    C_SM <- x$dissimilarity_matrix$Horn
    
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
        if(i==j){cat(round(0,0),'  ')}
        if(i<j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
    cat('\n')
  }
  if(q == 2){
    cat('    -----------------------Measure 1-C22------------------------\n\n')
    cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
    ###################################################CqN_ Equal weight
    Cqn_PC <- x$pairwise$C22
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
        if(temp[1]>1)
        {cat(')','   ',sprintf("%.3f",1)      ,'#        ',sprintf("%.3f",temp[2]),'         (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','   ',sprintf("%.3f",temp[1]),'         ',sprintf("%.3f",temp[2]),'         (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise dis-similarity matrix: \n\n')
    C_SM=x$dissimilarity_matrix$C22
    cat('    1-C22(i,j)')
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
        if(i==j){cat(round(0,0),'  ')}
        if(i<j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    cat('    -----------------------Measure 1-U22------------------------\n\n')
    cat('    Estimator','     Estimate','       s.e.','       95% Confidence Interval\n\n')
    Cqn_PC <- x$pairwise$U22
    no.temp=1
    for(i in 1:(N-1))
    {
      for(j in (i+1):N)
      {
        temp=Cqn_PC[no.temp,]
        cat('    1-U22(')
        cat(i)
        cat(',')
        cat(j)
        if(temp[1]>1)
        {cat(')','   ',sprintf("%.3f",1)      ,'#        ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
        }
        if(temp[1]<=1)
        {cat(')','   ',sprintf("%.3f",temp[1]),'         ',sprintf("%.3f",temp[2]),'        (',
             sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')}
        no.temp=no.temp+1
      }
    }
    cat('\n')
    cat('    Average pairwise dis-similarity =',sprintf("%.3f",mean(Cqn_PC[,1])),'\n\n')
    cat('    Pairwise dis-similarity matrix: \n\n')
    C_SM <- x$dissimilarity_matrix$U22
    
    cat('    1-U22(i,j)')
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
        if(i==j){cat(round(0,0),'  ')}
        if(i<j) {
          if(C_SM[i,j]<=1) cat(sprintf("%.3f",C_SM[i,j]),'  ')
          else cat(sprintf("%.3f#",1),'  ')
        }
      }
      cat('\n')
    }
    cat('\n')
    cat('    NOTE: Any estimate greater than 1 is replaced by 1; any estimate less than 0 is replaced by 0.')
    cat('\n')
  } 
}