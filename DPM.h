#ifndef RCPP_BASIC
#define RCPP_BASIC
#include "RcppBasic.h"
#endif

using namespace Rcpp;

List update_pi(IntegerVector S, int K, double gam){
  int n = S.size();
  NumericVector n_group(K), nu(K), Pi(K);
  
  for(int k = 0; k<K; k++){
    n_group[k] = sum(S==(k+1));
  }
  NumericVector n_cum = cumsum(n_group);
  
  for(int k = 0; k<(K-1); k++){
    nu[k] = rbeta(1, n_group[k]+1, gam + n - n_cum[k])[0];
  }
  NumericVector nu_prod = cumprod(1-nu);
  nu[K-1] = 1;
  
  Pi[0] = nu[0];
  for(int i = 1; i<K; i++){
    Pi[i] = nu[i]*nu_prod[i-1];
  }
  Pi [Pi < pow(10, -6)] = pow(10, -6) ;
  Pi = Pi/sum(Pi);
  return Rcpp::List::create(Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("Pi_tilde") = nu);
}


IntegerVector update_S(NumericVector Pi, NumericVector y, NumericVector m, NumericVector s2){
  int n = y.size(), K = m.size();
  
  IntegerVector S(n);
  NumericVector prob_S(K);;
  for(int i=0; i<n; i++){
    
    for(int k=0; k<K; k++){
      
      prob_S[k] = log(Pi[k]) + R::dnorm(y[i], m[k], sqrt(s2[k]), true);
    }
    
    prob_S = prob_S - max(prob_S) ;
    prob_S = exp_vec(prob_S);
    prob_S = prob_S/sum(prob_S);
    
    S[i] = sample(seq_len(K), prob_S, K);
  }
  return(S);
}


NumericMatrix update_phi(NumericVector y, IntegerVector S, double mu, 
                         double nu, double alpha, double lambda, int K){
  int n = S.size();
  NumericVector n_group(K), mean_y(K), ssq_y(K), mu_group(K), nu_group(K), alpha_group(K), lambda_group(K);
  NumericMatrix phi(K, 2);
  
  
  for(int k = 0; k<K; k++){
    
    n_group[k] = sum(S==(k+1));
    mean_y[k] = sum(y[S == (k+1)])/n_group[k];
    ssq_y[k] = sum(y[S == (k+1)]*y[S == (k+1)]) - (mean_y[k]*mean_y[k]*n_group[k]);
    
    if(n_group[k] == 0){
      mean_y[k] = 0.0;
      ssq_y[k] = 0.0;
      
      mu_group[k] = mu;
      nu_group[k] = nu;
      alpha_group[k] = alpha;
      lambda_group[k] = lambda;
    }
    else {
      
      nu_group[k] = nu + n_group[k];
      mu_group[k] = (mu*nu + n_group[k]*mean_y[k])/(nu_group[k]);
      alpha_group[k] = alpha + (0.5 *n_group[k]);
      lambda_group[k] = lambda + 0.5*ssq_y[k] + (0.5*nu*n_group[k]*((mean_y[k] - mu)*(mean_y[k] - mu)))/nu_group[k];
    }
    
    phi(k,1) = 1.0/rgamma(1, alpha_group[k], 1.0/lambda_group[k])[0] ;
    phi(k,0) = rnorm(1, mu_group[k], sqrt(phi(k,1)/nu_group[k]))[0] ;
  }
  
  return(phi);
}


double update_gam(NumericVector Pi_tilde, int K, double aa, double bb){
  double res, y = 1.0;
  
  for(int i=0; i<(K-2); i++)
    y = y*(1-Pi_tilde[i]);
  
  res = rgamma(1, (K-1+aa), 1.0/(bb - log(y)))[0];
  
  return(res);
}
