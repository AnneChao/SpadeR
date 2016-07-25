#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
double shared_entropy(NumericVector x, NumericVector y,int z,int n1,int n2,double w1)
{
     double sum=0;
     int i,t1,t2,k;
     double temp3=0;
     double w2=1-w1;
     for(i=0;i<z;i++)
     {
         for(t1=0;t1<=(n1-1);t1++) 
         {
             for(t2=0;t2<=n2;t2++)
             {
                 k=t1+t2;
                 if(k>0 && (n1-x[i]>=t1) &&  (n2-y[i]>=t2) )
                 {
                   {temp3=exp((t1+1)*log(w1)+(t2)*log(w2)+lgamma(k+1)-lgamma(t1+1)-lgamma(t2+1)+lgamma(n1-x[i]+1)+lgamma(n1-t1)-lgamma(n1-x[i]-t1+1)-lgamma(n1)+lgamma(n2-y[i]+1)+lgamma(n2-t2+1)-lgamma(n2-y[i]-t2+1)-lgamma(n2+1));}         
                   sum=sum+1.0/k*x[i]/n1*temp3;
                 }      
             }
         }       
     }
     return(sum);
}


//[[Rcpp::export]]
double shared_entropy_sub(NumericVector X1, NumericVector X2,NumericVector p1, NumericVector p2,int I,int n1,int n2,double w1)
{
     double sum1=0.0;
     double sum2=0.0;
     double temp1=0.0;
     int i,t1,t2,k;
     double w2=1-w1;
     for(i=0;i<I;i++)
     {
         for(t1=n1;t1<=n1+n2-2;t1++)
         {
             for(t2=0;t2<=n1+n2-2-t1;t2++)
             {
                 if( (n2-X2[i])>=t2    )
                 {
                     k=t1+t2;
                     temp1=exp((t1+1)*log(w1)+(t2)*log(w2)+lgamma(k+1)-lgamma(t1+1)-lgamma(t2+1));
                     sum1=sum1+1.0*temp1/(t1+t2)*p1[i]*pow(1-p1[i],t1)*exp(lgamma(n2-X2[i]+1)+lgamma(n2-t2+1)-lgamma(n2-X2[i]-t2+1)-lgamma(n2+1));
                 }
             }
         }
         for(t1=0;t1<=n1-3;t1++)
         {
             for(t2=n2+1;t2<=n1+n2-2-t1;t2++)
             {
                 if( (n1-X1[i])>=t1 )
                 {
                      k=t1+t2;
                      temp1=exp((t1+1)*log(w1)+(t2)*log(w2)+lgamma(k+1)-lgamma(t1+1)-lgamma(t2+1));
                      sum2=sum2+1.0*temp1/(t1+t2)*X1[i]/n1*pow(1-p2[i],t2)*exp(lgamma(n1-X1[i]+1)+lgamma(n1-t1)-lgamma(n1-X1[i]-t1+1)-lgamma(n1));
                 }
             }
         }
     }
     return(sum1+sum2);     
}