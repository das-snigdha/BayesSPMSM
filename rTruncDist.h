#ifndef RCPP_BASIC
#define RCPP_BASIC
#include "RcppBasic.h"
#endif

using namespace Rcpp ;

NumericMatrix rtmvnormHMC
(int n, NumericVector mu, NumericMatrix Sigma, NumericVector init,
 NumericMatrix ff, NumericVector gg, int n_burn=0)
{
  //// constraints: innderProduct(ff(j,_), x) + gg[j] \geq 0, j=1,...,m
  
  int d=mu.size(), m=ff.nrow(), h ;
  double u, phi, tmp, tmp2, T_end, T_h, alpha_h ;
  NumericVector x(d), s(d), a(d), b(d), T(m), x_dot(d), g(m) ;
  NumericMatrix f(m,d), res(n,d), Sigma_chol_L(d,d), Sigma_chol_U(d,d) ;
  
  chol(Sigma, Sigma_chol_U) ;
  transpose(Sigma_chol_U, Sigma_chol_L) ;
  for(int j=0; j<m; j++){
    g[j] = innerProduct(ff(j,_), mu) + gg[j] ;
    for(int i=0; i<d; i++)
      f(j,i) = innerProduct(Sigma_chol_U(i,_), ff(j,_)) ;
  }
  
  product(solve(Sigma_chol_L), init-mu, x) ;
  
  for(int nn=0; nn<n+n_burn; nn++){
    
    s = rnorm(d) ;
    vec_mem_cpy(s, a) ;
    vec_mem_cpy(x, b) ;
    
    T_end = M_PI/2.0 ;
    
    while(1){
      for(int j=0; j<m; j++){
        tmp = innerProduct(f(j,_), a) ;
        tmp2 = innerProduct(f(j,_), b) ;
        u = sqrt(tmp*tmp + tmp2*tmp2) ;
        phi = atan2(-1.0*tmp, tmp2) ;
        
        if((u < g[j]) || (u < -1.0*g[j]))
          T[j] = T_end ;
        else{
          T[j] = acos(-1.0*g[j]/u) - phi ;
          if(T[j] < 0.0){
            T[j] += 2.0*M_PI ;
            Rprintf("1\n") ;
          }
        }
      }
      
      h = 0 ;
      T_h = T[0] ;
      for(int j=1; j<m; j++){
        if(T[j] < T_h){
          T_h = T[j] ;
          h = j ;
        }
      }
      
      if(T_h < T_end){
        for(int i=0; i<d; i++){
          tmp = sin(T_h - 1e-10) ;
          tmp2 = cos(T_h - 1e-10) ;
          x[i] = a[i]*tmp + b[i]*tmp2 ;
          x_dot[i] = a[i]*tmp2 - b[i]*tmp ;
        }
        alpha_h = innerProduct(f(h,_), x_dot) / norm(f(h,_),2,true) ;
        for(int i=0; i<d; i++)
          a[i] = x_dot[i] - 2.0*alpha_h*f(h,i) ;
        vec_mem_cpy(x, b) ;
        T_end -= T_h ;
      }
      else{
        for(int i=0; i<d; i++){
          tmp = sin(T_end) ;
          tmp2 = cos(T_end) ;
          x[i] = a[i]*tmp + b[i]*tmp2 ;
        }
        break ;
      }
    }
    
    if(nn >= n_burn)
      for(int i=0; i<d; i++)
        res(nn-n_burn,i) = innerProduct(Sigma_chol_L(i,_), x) + mu[i] ;
  }
  
  return(res) ;
}


double rtruncnorm_c(double m, double s, double lower, double upper){
  
  // obtain environment containing function
  Environment truncnorm = Rcpp::Environment::namespace_env("truncnorm");
  // make function callable from Cpp
  Function rtruncnorm = truncnorm["rtruncnorm"];
  
  // return object
  // code is interpreted as dbvnorm(x = a, y = b, cor = p)
  return as<double>(rtruncnorm(Named("n") = 1, _["a"] = lower, _["b"] = upper, _["mean"] = m, _["sd"] = s));
}


double rtruncbeta(double shape1, double shape2, double lower, double upper){
  double res, u, a, b;
  
  a = R::pbeta(lower, shape1, shape2, true, false) ;
  b = R::pbeta(upper, shape1, shape2, true, false) ;
  u = runif(1, a, b)[0];
  
  res = R::qbeta(u, shape1, shape2, true, false) ;
  
  return(res) ;
}