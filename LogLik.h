#ifndef RCPP_BASIC
#define RCPP_BASIC
#include "RcppBasic.h"
#endif

using namespace Rcpp ;

double F_dpm(double y, NumericVector mm, NumericVector ss, NumericVector Pi){
  
  double res = 0.0;
  int L = mm.size();
  
  for(int k=0; k<L; k++)
    res += (double) R::pnorm(y, mm[k], sqrt(ss[k]), true, false) * Pi[k];
  
  return(res);
}

NumericVector calc_loglik(int n, int m, NumericVector mm, NumericVector ss, NumericVector Pi,
                          NumericMatrix R,
                          NumericVector b, NumericVector gSI, NumericVector logC, IntegerVector State){
  
  int N = n*m, i, K = max(State), L = mm.size();
  double p, y0, a0, b0, R_sum ;
  NumericVector logl(n) ;
  
  for(int kk=0; kk<n; kk++){
    
    logl[kk] = 0.0;
    
    for(int jj=0; jj<m; jj++){
      
      i = (kk*m) + jj ;
      p = 0.0;
      
      if(State[i] == 0){
        y0 = logC[i] - log(R(i, 0)) - b[i] -gSI[i] ;
        
        p = 1 - F_dpm(y0, mm, ss, Pi) ;
      }
      else if(State[i] <K){
        
        R_sum = 0.0 ;
        for(int j=0; j<State[i]; j++)
          R_sum += R(i, j) ;
        
        b0 = logC[i] - log(R_sum) - b[i] -gSI[i] ;
        a0 = logC[i] - log(R_sum + R(i, State[i])) - b[i] -gSI[i] ;
        
        p = F_dpm(b0, mm, ss, Pi)  - F_dpm(a0, mm, ss, Pi) ;
      }
      
      else{
        y0 = logC[i] - b[i] -gSI[i] ;
        
        p = F_dpm(y0, mm, ss, Pi) ;
      }
      
      logl[kk] += log(p) ;
      
    }
  }
  
  return(logl);
  
}