print.spadeGenetic <- function(x, ...)
{
  cat('\n(1) BASIC DATA INFORMATION:\n\n')
  cat('    The loaded set includes abundance (or frequency) data from',x$info[1],'subpopulations\n')
  cat('    and a total of',x$info[2],'distinct alleles.\n\n')
  cat('   (Sample size in each subpopulation)                 n1 =', x$info[3],'\n')
  N <- x$info[1]
  q <- x$q
  for(j in 2:N){
    cat('                                                      ','n')
    cat(j,'=',x$info[2+j],'\n')   
  }
  cat('\n')
  cat('   (Number of alleles in one subpopulation)            D1 =', x$info[N+3],'\n')
  
  for(j in 2:N){
    cat('                                                      ','D')
    cat(j,'=',x$info[N+2+j],'\n')
  }
  cat('\n')
  cat('   (Number of shared alleles in two subpopulations)  ',' D12 =', x$info[3+2*N], '\n')
  
  if(N>2){
    k <- 1
    for(i in 1:(N-1)){     
      for(j in (i+1):N){
        if(i==1 & j==2) next
        cat('                                                      ','D')
        cat(i,j,' = ', x$info[3+2*N+k], '\n', sep="")
        k <- k + 1
      }
    }
  }
  cat('\n')
  if(N==3)
  {
    cat('   (Number of shared alleles in three subpopulations)',' D123','=',rev(x$info)[2],'\n\n') 
  }
  cat('   (Bootstrap replications for s.e. estimate)         ',rev(x$info)[1],'\n\n')
  cat('(2) NEARLY UNBIASED ESTIMATION OF ALLELIC DIFFERENTIATION OR MORISITA DISSIMILARITY IN ',N,'SUBPOPULATIONS:\n\n')
  cat('    Estimator','     Estimate','     Est_s.e.','     95% Confidence Interval\n\n')
  if ( q==0 ) temp0n=x$overlap[1,]
  if ( q==1 ) temp0n=x$overlap[3,]
  if ( q==2 ) temp0n=x$overlap[4,]
  if ( temp0n[1]<=1 ){
    cat('    1-C')
    cat(q)
    cat(N,sprintf("         %.3f",1-temp0n[1]),'       ',sprintf("%.3f",temp0n[2]),'        (',
        sprintf("%.3f",1-temp0n[4]),',',sprintf("%.3f",1-temp0n[3]),')\n')
  }
  if ( temp0n[1]>1 ){
    cat('    1-C')
    cat(q)
    cat(N,sprintf("         %.3f#",0),'      ',sprintf("%.3f",temp0n[2]),'        (',
        sprintf("%.3f",1-temp0n[4]),',',sprintf("%.3f",1-temp0n[3]),')\n')
  }
  
  cat('\n')
  cat('    1-C')
  cat(q)
  cat(N,':This is the genetic diversity measure defined in Jost (2008) for comparing subpopulations based 
      on allele shared information between any two subpopulations.\n')
  cat('    
      Confidence Interval: Based on an improved bootstrap percentile method. (recommend for use in the case when 
      similarity is close to 0 or 1 ) \n\n')
  cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
      is greater than 1, replace it by 1.\n\n')
  cat('    Pairwise Comparison:\n\n')
  cat('    Estimator','          Estimate','     Est_s.e.','     95% Confidence Interval\n\n')
  Cqn_PC <- x$pairwise  
  no.temp=1
  share_index = sapply(1:choose(N,2), function(i){x$info[2+2*N+i]})
  no.temp2=1
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      temp=Cqn_PC[no.temp,]
      cat('    1-C')
      cat(q)
      cat('2(')
      cat(i)
      cat(',')
      cat(j)
      if ( share_index[no.temp2]!=0 )
      {
        if ( temp[1]>1 ){
          cat(')','        ',sprintf("%.3f#",0),'      ',sprintf("%.3f",temp[2]),'        (',
              sprintf("%.3f",1-temp[4]),',',sprintf("%.3f",1-temp[3]),')\n')
        }
        if ( temp[1]<=1 ){
          cat(')','        ',sprintf("%.3f",1-temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
              sprintf("%.3f",1-temp[4]),',',sprintf("%.3f",1-temp[3]),')\n')
        }
      }
      else
      {
        cat(')','        ',sprintf("%.3f",1),'       ',sprintf("%.3f",0),'        (',
            sprintf("%.3f",1),',',sprintf("%.3f",1),')\n')
        #         cat(')','        ',sprintf("%.3f##",1),'     ',sprintf("%.3f##",0),'      (',
        #             sprintf("%.3f",1),',',sprintf("%.3f",1),')##\n')
      }
      no.temp=no.temp+1
      no.temp2=no.temp2+1
    }
  }
  cat('\n')
  cat('    Average Pairwise =',sprintf("%.3f",1-mean(Cqn_PC[1:choose(N,2),1])),'\n')
  cat('    ## There are no shared species, thus estimated similarity is zero and should be used for caution.\n\n')
  cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
      is greater than 1, replace it by 1.\n\n')
  cat('    1-C')
  cat(q)
  cat('2: This is the genetic diversity measure D defined in Jost (2008) for comparing 2 subpopulations.')
  cat('\n\n')
  cat('    Dissimilarity Matrix: \n\n')
  C_SM=x$similarity.matrix
  cat('    1-C')
  cat(q)
  cat('2(i,j)\t')
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
      if(i<=j){
        if (1-C_SM[i,j]>=0) cat(sprintf("%.3f",abs(1-C_SM[i,j])),'\t')
        else cat(sprintf("%.3f#",0),'\t')
      }
    }
    cat('\n')
  }
  cat('\n')
  cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
      is greater than 1, replace it by 1.\n\n')
  
  cat('(3)  NEARLY UNBIASED ESTIMATION OF MORISITA SIMILARITY IN ',N,'SUBPOPULATIONS:\n\n')
  cat('    Estimator','     Estimate','     Est_s.e.','     95% Confidence Interval\n\n')
  
  if ( q==0 ) temp2n=x$overlap[1,]
  if ( q==1 ) temp2n=x$overlap[3,]
  if ( q==2 ) temp2n=x$overlap[4,]
  if ( temp2n[1]>=1 ){
    cat('    C')
    cat(q)
    cat(N,sprintf("           %.3f#",1),'      ',sprintf("%.3f",temp2n[2]),'       (',
        sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')
  }
  if ( temp2n[1]<1 ){
    cat('    C')
    cat(q)
    cat(N,sprintf("           %.3f",temp2n[1]),'       ',sprintf("%.3f",temp2n[2]),'       (',
        sprintf("%.3f",temp2n[3]),',',sprintf("%.3f",temp2n[4]),')\n')
  }
  cat('\n')
  cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
      is greater than 1, replace it by 1.\n\n')
  cat('\n')
  cat('    C')
  cat(q)
  cat(N,': A similarity measure of comparing 3 subpopulations based on allele shared information between any two subpopulations.\n')
  cat('\n')
  cat('    Pairwise Comparison:\n\n')
  cat('    Estimator','     Estimate','     Est_s.e.','     95% Confidence Interval\n\n')
  Cqn_PC <- x$pairwise
  no.temp=1
  share_index = sapply(1:choose(N,2), function(i){x$info[2+2*N+i]})
  no.temp2=1
  for(i in 1:(N-1))
  {
    for(j in (i+1):N)
    {
      temp=Cqn_PC[no.temp,]
      if(q==0){cat('    C02(')}
      if(q==1){cat('    C12(')}
      if(q==2){cat('    C22(')}
      cat(i)
      cat(',')
      cat(j)
      if ( share_index[no.temp2]!=0 ){
        cat(')','     ',sprintf("%.3f",temp[1]),'       ',sprintf("%.3f",temp[2]),'        (',
            sprintf("%.3f",temp[3]),',',sprintf("%.3f",temp[4]),')\n')
      }
      else
      {
        cat(')','     ',sprintf("%.3f",0),'       ',sprintf("%.3f",0),'        (',
            sprintf("%.3f",0),',',sprintf("%.3f",0),')\n')
      }
      Cqn_PC[i,1]=0
      no.temp=no.temp+1
      no.temp2=no.temp2+1
    }
  }
  cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
      is greater than 1, replace it by 1.\n\n')
  cat('\n')
  cat('    Average Pairwise =',sprintf("%.3f",mean(Cqn_PC[1:choose(N,2),1])),'\n')
  cat('\n')
  cat('    ## There are no shared species, thus estimated similarity is zero and should be used for caution.\n\n')
  
  cat('    Similarity Matrix: \n\n')
  C_SM=x$similarity.matrix
  if(q==0){cat('    C02(i,j)   \t')}
  if(q==1){cat('    C12(i,j)   \t')}
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
      if(i<=j){
        if (C_SM[i,j]<=1) cat(sprintf("%.3f",abs(C_SM[i,j])),'\t')
        else cat(sprintf("%.3f#",1),'\t')
      }
    }
    cat('\n')
  }
  cat('    # If an estimator is less than 0, replace it by 0; if an estimator 
      is greater than 1, replace it by 1.\n\n')
  cat('\n')
  cat('    References:\n')
  cat('    Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H. and Chazdon, R. (2008). A Two-
      stage probabilistic approach to multiple-community similarity indices. 
      Biometrics, 64, 1178-1186.\n')
  cat('    Jost, L. (2008). GST and its relatives do not measure differentiation. Molecular
      Ecology, 17, 4015-4026.')
}
