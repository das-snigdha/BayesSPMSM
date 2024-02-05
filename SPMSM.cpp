// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

#define RCPP_BASIC
#include "RcppBasic.h"
#define RTRUNCDIST
#include "rTruncDist.h"
#include "singleIndex.h"
#include "DPM.h"
#include "LogLik.h"

using namespace Rcpp ;

// [[Rcpp::export]]
double g(double x, NumericVector xi, NumericVector u){
  // g(x) = \sum_{l in 0:L} xi_(l+1)*psi_l(x)
  if(x > 1.0)
    return(g(1.0, xi, u)) ; 
  else if(x < -1.0)
    return(0.0) ;
  else{
    double res = 0.0 ;
    for(int l=0; l<xi.size();l++)
      res += xi[l]*psi_l(x,l,u) ;
    return(res) ;
  }
}

double log_post_alpha(NumericVector alpha, NumericMatrix X, double a, double lambda){
  
  int n=X.nrow(), K=X.ncol() ;
  double lgamma_sum, sum_lgamma_alpha, sum_alpha, sum_alpha_logX, sum_log_alpha ;
  double res ;
  NumericVector colSum_logX(K) ;
  
  for(int k=0; k<K; k++){
    for(int i=0; i<n; i++){
      colSum_logX[k] += log(X(i,k)) ;
    }
    
    sum_lgamma_alpha += lgamma(alpha[k]) ;
    sum_alpha_logX += (alpha[k]-1) * colSum_logX[k];
    sum_log_alpha += log(alpha[k]); 
  }
  
  sum_alpha = sum(alpha) ;
  lgamma_sum = lgamma(sum_alpha) ;
  
  res = n*(lgamma_sum - sum_lgamma_alpha) + sum_alpha_logX - (lambda*sum_alpha) + (a-1)*sum_log_alpha ;
  return(res);
  
}

List Bayesian_gamma_Dirichlet
(NumericMatrix X, NumericVector init_alpha, double step_size, double a=1.0, double lambda=0.01){
  
  // Sampling from the posterior distribution
  // with a random walk in a Gibbs sampler
  // X(i,_) ~ Dirichlet(alhpa_k), k=1,...,K
  // alpha_k ~ Gamma(a,lambda) (indep)
  
  int n=X.nrow(), K=X.ncol() ;
  double log_u, log_r, dec ;
  NumericVector new_alpha(K), res_alpha(K), eps(K) ;
  
  eps = rnorm(K, 0, step_size);
  new_alpha = init_alpha + eps;
  
  
  if(is_true(any(new_alpha < 0))){
    vec_mem_cpy(init_alpha, res_alpha); 
    dec = 0.0;
  }
  else{
    log_r = log_post_alpha(new_alpha, X, a, lambda) - log_post_alpha(init_alpha, X, a, lambda) ;
    if(log_r > 0){
      log_r = 0.0;
    }
    log_u = log(runif(1)[0]);
    if(R_IsNA(log_r)){
      vec_mem_cpy(init_alpha, res_alpha);
      dec = 0.0;
    } 
    else if(log_u < log_r){
      vec_mem_cpy(new_alpha, res_alpha) ;
      dec = 1.0;
    }
    else {
      vec_mem_cpy(init_alpha, res_alpha);
      dec = 0.0;
    }
  }
  return Rcpp::List::create(Rcpp::Named("alpha") = res_alpha,
                            Rcpp::Named("dec") = dec);
}

void update_b(NumericMatrix mu_mtx, NumericMatrix Sigma_b, NumericVector sigma_sq_eps, NumericVector b_mtx){
  int n=mu_mtx.nrow(), m=mu_mtx.ncol() ;
  NumericVector mean_b(m/2), tmp_m(m/2) ;
  NumericMatrix tmp_mm(m/2,m/2), Var_b(m/2,m/2), Sigma_eps(n,m) ;
  
  asMatrix(sigma_sq_eps, Sigma_eps, true) ;
  
  for(int i=0; i<n; i++){
    solve(Sigma_b, tmp_mm) ;
    for(int j=0; j<m/2; j++)
      tmp_mm(j,j) += 1.0/Sigma_eps(i,j) ;
    for(int j=m/2; j<m; j++)
      tmp_mm(m-1-j,m-1-j) += 1.0/Sigma_eps(i,j) ;
    solve(tmp_mm, Var_b) ;
    
    for(int j=0; j<m/2; j++)
      tmp_m[j] = mu_mtx(i,j)/Sigma_eps(i,j) ;
    for(int j=m/2; j<m; j++)
      tmp_m[m-1-j] += mu_mtx(i,j)/Sigma_eps(i,j) ;
    product(Var_b, tmp_m, mean_b) ;
    rmvnorm(mean_b, Var_b, tmp_m) ;
    for(int j=0; j<m/2; j++)
      b_mtx(i,j) = tmp_m[j] ;
  }
}

// [[Rcpp::export]]
List CS_MCMC
(IntegerVector id, IntegerVector tooth, NumericVector C, IntegerVector State, NumericMatrix X,
 NumericVector u, NumericMatrix S_Omega, int c_Omega, 
 NumericVector init_xi, NumericVector init_beta, NumericVector init_alpha,
 NumericMatrix init_Sigma_b, NumericVector init_logT,
 NumericMatrix Kinv, double l_matern, double nu_matern = 2.5, double eta_xi = 50, 
 double a_xi = 3.0, double b_xi = 100, double sigma_sq_beta = 10.0, 
 double alpha_a = 1.0, double alpha_lambda = 0.1, double alpha_stepsize = 0.2,
 double f_mu = 0.0, double f_nu = 1.0, double f_lambda = 1.0, double f_alpha = 3.0,
 double f_aa = 0.1, double f_bb = 0.1, int Lmax = 5,
 int M_mcmc=1000, int M_burn=100, int thin=1)
{
  // Obs must be sorted by (id, tooth) in advance
  int n=max(id), n_tooth=max(tooth), N=n*n_tooth, p=X.ncol(), K=max(State), Kc ;
  int L=u.size()-1, M = M_mcmc/thin ;
  double alpha_dec, tmp1, tmp2, tmp3;
  
  List alpha_out, pi_out;
  
  NumericVector mu(N), eps(N), b(N), alpha(K), beta(p), beta_tilde(p), logT(N), logT_rep(N), logT_rep_sq(N) ;
  NumericVector logC(N), Z(N), SI(N), gSI(N), xi(L+1) ;
  NumericVector sigma_sq_eps(N,0.01), tmp_n(n), tmp_K(K) ;
  NumericMatrix mu_mtx(n,n_tooth), b_mtx(n,n_tooth/2), Sigma_b(n_tooth/2,n_tooth/2) ;
  NumericMatrix D_Omega(n_tooth/2,n_tooth/2), R(N,K) ;
  
  IntegerVector K_vec(M), labels(N) ;
  NumericVector gam_vec(M), loglik(n), alpha_dec_vec(M);
  NumericMatrix eps_mtx(M, N), xi_mtx(M,L+1), gSI_mtx(M, N),  beta_tilde_mtx(M,p), beta_mtx(M,p), phi(Lmax, 2) ;
  NumericMatrix alpha_mtx(M,K), loglik_mtx(M, n), Sigma_b_mtx(M, n_tooth/2), hat_Sigma_b(n_tooth/2,n_tooth/2);
  NumericMatrix pi_mtx(M, Lmax), mm_mtx(M, Lmax), ss_mtx(M, Lmax);
  
  // initialization
  
  vec_mem_cpy(init_xi, xi) ;
  vec_mem_cpy(init_alpha, alpha) ;
  
  vec_mem_cpy(init_beta, beta) ;
  product(X, beta, SI) ;
  for(int i=0; i<N; i++)
    gSI[i] = g(SI[i], xi, u) ;
  for(int i=0; i<N; i++)
    logC[i] = log(C[i]) ;
  
  mtx_mem_cpy(init_Sigma_b, Sigma_b) ;
  for(int i=0; i<N; i++)
    b[i] = 1e-5 ;
  vec_mem_cpy(init_logT, logT) ;
  
  singleIndex_varyingSigma SI_obj(mu, X, u, nu_matern, l_matern, Kinv, eta_xi, beta, xi, 
                                  rep(1.0,N), sigma_sq_beta, a_xi, b_xi) ;
  
  NumericVector Pi = rep(1.0/Lmax, Lmax);
  NumericVector ss = 1.0/rgamma(Lmax, f_alpha, 1.0/f_lambda);
  NumericVector mm = rnorm(Lmax);
  
  double gam = rgamma(1, f_aa, 1.0/f_bb)[0];
  
  // run MCMC
  for(int m=0; m<M_mcmc+M_burn; m++){
    
    if(m%200==0){
      Rprintf("Number of iterations = %d\n", m) ;
      Rprintf("beta: ") ;
      print(beta) ;
    }
    
    // update R
    for(int i=0; i<N; i++){
      if(State[i]==0){
        tmp1 = alpha[0] ;
        tmp2 = sum(alpha) - alpha[0] ;
        R(i,0) = rtruncbeta(tmp1, tmp2, C[i]/exp(logT[i]), 1.0) ; 
        vec_mem_cpy(rDirichlet(1, alpha[Range(1,K-1)])(0,_), tmp_K, K-1) ;
        for(int k=1; k<K; k++)
          R(i,k) = (1.0-R(i,0))*tmp_K[k-1] ;
      }
      else if(State[i] < K-1){
        double V1, V2 ;
        tmp1 = 0.0 ;
        for(int k=0; k<State[i]; k++)
          tmp1 += alpha[k] ;
        tmp2 = sum(alpha) - tmp1 ;
        V1 = rtruncbeta(tmp1, tmp2, 0.0, C[i]/exp(logT[i])) ;
        tmp2 -= alpha[State[i]] ;
        tmp1 = rtruncbeta(alpha[State[i]], tmp2, (C[i]/exp( logT[i])-V1)/(1.0-V1), 1.0) ;
        R(i,State[i]) = (1.0-V1)*tmp1 ;
        V2 = 1.0 - V1 - R(i,State[i]) ;
        vec_mem_cpy(rDirichlet(1, alpha[Range(0,State[i]-1)])(0,_), tmp_K, State[i]) ;
        for(int k=0; k<State[i]; k++)
          R(i,k) = V1*tmp_K[k] ;
        vec_mem_cpy(rDirichlet(1, alpha[Range(State[i]+1,K-1)])(0,_), tmp_K, K-State[i]-1) ;
        for(int k=State[i]+1; k<K; k++)
          R(i,k) = V2*tmp_K[k-State[i]-1] ;
      }
      else if(State[i] < K){
        double V ;
        tmp1 = 0.0 ;
        for(int k=0; k<State[i]; k++)
          tmp1 += alpha[k] ;
        tmp2 = sum(alpha) - tmp1 ;
        V = rtruncbeta(tmp1, tmp2, 0.0, C[i]/exp(logT[i])) ;
        R(i,K-1) = 1.0 - V ;
        vec_mem_cpy(rDirichlet(1, alpha[Range(0,K-2)])(0,_), tmp_K, K-1) ;
        for(int k=0; k<K-1; k++)
          R(i,k) = V*tmp_K[k] ;
      }
      else{
        vec_mem_cpy(rDirichlet(1, alpha)(0,_), tmp_K) ;
        for(int k=0; k<K; k++)
          R(i,k) = tmp_K[k] ;
      }
    }
    
    // update logT
    for(int i=0; i<N; i++){
      mu[i] = b[i] + gSI[i] + Z[i] ;
      //mu[i] = b[i] + init_gSI[i] + Z[i] ;
      
      if(State[i] == 0){
        tmp1 = logC[i] - log(R(i,0)) ;
        tmp2 = R_PosInf ;
        
        logT[i] = rtruncnorm_c(mu[i], sqrt(sigma_sq_eps[i]), tmp1, tmp2) ;
        
      }
      else if(State[i] < K){
        tmp3 = 0.0 ;
        for(int k=0; k<State[i]; k++)
          tmp3 += R(i,k) ;
        tmp2 = logC[i] - log(tmp3) ;
        tmp1 = logC[i] - log(tmp3 + R(i,State[i])) ;
        
        logT[i] = rtruncnorm_c(mu[i], sqrt(sigma_sq_eps[i]), tmp1, tmp2) ;
      }
      else{
        tmp1 = R_NegInf ;
        tmp2 = logC[i] ;
        
        logT[i] = rtruncnorm_c(mu[i], sqrt(sigma_sq_eps[i]), tmp1, tmp2) ;
      }
    }
    
    // update b
    mu = logT - Z - gSI ;
    asMatrix(mu, mu_mtx, true) ;
    update_b(mu_mtx, Sigma_b, sigma_sq_eps, b_mtx) ;
    for(int i=0, ind=0; i<n; i++){
      for(int j=0; j<n_tooth/2; j++, ind++)
        b[ind] = b_mtx(i,j) ;
      for(int j=n_tooth/2; j<n_tooth; j++, ind++)
        b[ind] = b_mtx(i,n_tooth-1-j) ;
    }
    
    // update DP
    mu = logT - b - gSI ;
    labels = update_S(Pi, mu, mm, ss) ;
    pi_out = update_pi(labels, Lmax, gam) ;
    Pi = pi_out["Pi"];
    phi = update_phi(mu, labels, f_mu, f_nu, f_alpha, f_lambda, Lmax);
    
    mm = phi(_,0);
    ss = phi(_,1);
    gam = update_gam(pi_out["Pi_tilde"], Lmax, f_aa, f_bb);
    
    for(int i=0; i<N; i++){
      Z[i] = phi(labels[i]-1,0) ;
      sigma_sq_eps[i] = phi(labels[i]-1,1) ;
    }
    
    eps = mu - Z ;
    SI_obj.change_sigma_sq_eps(sigma_sq_eps) ;
    
    // update g, beta
    mu = logT - b - Z ;
    SI_obj.change_y(mu) ;
    SI_obj.update_xi() ;
    
    vec_mem_cpy(SI_obj.xi, xi) ;
    SI_obj.update_beta() ;
    vec_mem_cpy(SI_obj.beta, beta) ;
    vec_mem_cpy(SI_obj.beta_tilde, beta_tilde) ;
    product(X, beta, SI) ;
    for(int i=0; i<N; i++)
      gSI[i] = g(SI[i], xi, u) ;
    
    // update Omega_b/Sigma_b
    AtA(b_mtx, D_Omega) ;
    for(int j=0; j<n_tooth/2; j++)
      for(int jj=0; jj<n_tooth/2; jj++)
        D_Omega(j,jj) += S_Omega(j,jj) ;
    rInvWishart(c_Omega+n, solve(D_Omega), Sigma_b) ;
    
    // update alpha
    alpha_out = Bayesian_gamma_Dirichlet(R, alpha, alpha_stepsize, alpha_a, alpha_lambda) ;
    alpha = alpha_out["alpha"] ;
    alpha_dec = alpha_out["dec"];
    
    loglik = calc_loglik(n, n_tooth, mm, ss, Pi, 
                         R, b, gSI, logC, State) ;
    
    // store sample
    if(m >= M_burn && m%thin == 0){  
      loglik_mtx((m-M_burn)/thin,_) = loglik ;
      gam_vec[(m-M_burn)/thin] = gam ;
      
      eps_mtx((m-M_burn)/thin,_) = eps;
      mm_mtx((m-M_burn)/thin,_) = mm;
      ss_mtx((m-M_burn)/thin,_) = ss;
      pi_mtx((m-M_burn)/thin,_) = Pi;
      
      xi_mtx((m-M_burn)/thin,_) = xi ;
      gSI_mtx((m-M_burn)/thin,_) = gSI ;
      beta_mtx((m-M_burn)/thin,_) = beta ;
      beta_tilde_mtx((m-M_burn)/thin,_) = beta_tilde ;
      alpha_mtx((m-M_burn)/thin,_) = alpha ;
      alpha_dec_vec[(m-M_burn)/thin] = alpha_dec ;
      
      Sigma_b_mtx((m-M_burn)/thin,_) = diag(Sigma_b);
      
      for(int j1=0; j1<n_tooth/2; j1++)
        for(int j2=0; j2<n_tooth/2; j2++)
          hat_Sigma_b(j1,j2) += Sigma_b(j1,j2)/(double)M ;
    }
  }
  
  List res ;
  res["loglik"] = loglik_mtx ;
  res["beta"] = beta_mtx ;
  res["beta_tilde"] = beta_tilde_mtx ;
  res["alpha"] = alpha_mtx ;
  res["alpha_dec"] = alpha_dec_vec ;
  
  res["xi"] = xi_mtx ;
  res["eps"] = eps_mtx ;
  res["gam"] = gam_vec ;
  res["pi"] = pi_mtx ;
  res["mm"] = mm_mtx ;
  res["ss"] = ss_mtx ;
  
  res["Sigma_b"] = Sigma_b_mtx ;
  res["hat_Sigma_b"] = hat_Sigma_b ;

  return(res) ;
}


// [[Rcpp::export]]
List CS_MCMC_BP
(IntegerVector id, IntegerVector tooth, NumericVector C, IntegerVector State, NumericMatrix X,
 NumericMatrix S_Omega, int c_Omega,
 NumericVector init_psi, NumericVector init_beta, NumericVector init_alpha,
 NumericMatrix init_Sigma_b, NumericVector init_logT, 
 double sigma_sq_beta = 10.0, double a_psi = 3.0, double b_psi = 1000,
 double alpha_a = 1.0, double alpha_lambda = 0.1, double alpha_stepsize = 0.2,
 double f_mu = 0.0, double f_nu = 1.0, double f_lambda = 1.0, double f_alpha = 0.01,
 double f_aa = 0.1, double f_bb = 0.1, int Lmax = 5,
 int M_mcmc=1000, int M_burn=100, int thin=1)
{
  // Obs must be sorted by (id, tooth) in advance
  int n=max(id), n_tooth=max(tooth), N=n*n_tooth, p=X.ncol(), K=max(State) ;
  int L=init_psi.size(), M = M_mcmc/thin ; //ngrid=grid_x.size() ;
  double alpha_dec, tmp1, tmp2, tmp3 ;
  
  List alpha_out, pi_out ;
  
  NumericVector mu(N), b(N), eps(N), alpha(K), beta(p), beta_tilde(p), logT(N), logT_rep(N), logT_rep_sq(N) ;
  NumericVector logC(N), Z(N), SI(N), gSI(N), psi(L) ;
  NumericVector sigma_sq_eps(N,0.01), tmp_n(n), tmp_K(K) ;
  NumericMatrix mu_mtx(n,n_tooth), b_mtx(n,n_tooth/2), Sigma_b(n_tooth/2,n_tooth/2) , hat_Sigma_b(n_tooth/2,n_tooth/2);
  NumericMatrix D_Omega(n_tooth/2,n_tooth/2), R(N,K), Db(N, L) ;
  
  IntegerVector K_vec(M), labels(N) ;
  NumericVector gam_vec(M), loglik(n), alpha_dec_vec(M);
  NumericMatrix eps_mtx(M, N), psi_mtx(M,L), gSI_mtx(M, N), beta_tilde_mtx(M,p), beta_mtx(M,p);//, eps_density(M,ngrid) ;
  NumericMatrix alpha_mtx(M,K), Sigma_b_mtx(M, n_tooth/2), loglik_mtx(M, n), phi(Lmax, 2) ;
  NumericMatrix pi_mtx(M, Lmax), mm_mtx(M, Lmax), ss_mtx(M, Lmax);
  
  vec_mem_cpy(init_psi, psi) ;
  vec_mem_cpy(init_alpha, alpha) ;
  
  vec_mem_cpy(init_beta, beta) ;
  product(X, beta, SI) ;
  
  Db = Dbeta(L, SI);
  product(Db, psi, gSI) ;
  
  for(int i=0; i<N; i++)
    logC[i] = log(C[i]) ;
  
  mtx_mem_cpy(init_Sigma_b, Sigma_b) ;
  for(int i=0; i<N; i++)
    b[i] = 1e-5 ;
  vec_mem_cpy(init_logT, logT) ;
  
  singleIndex_BP SI_obj(mu, X, beta, psi, rep(1.0, N), sigma_sq_beta, a_psi, b_psi) ;
  
  NumericVector Pi = rep(1.0/Lmax, Lmax);
  NumericVector ss = 1.0/rgamma(Lmax, f_alpha, 1.0/f_lambda);
  NumericVector mm = rnorm(Lmax);
  
  double gam = rgamma(1, f_aa, 1.0/f_bb)[0];
  
  
  // run MCMC
  for(int m=0; m<M_mcmc+M_burn; m++){
    
    if(m%200==0){
      Rprintf("Number of iterations = %d\n", m) ;
      Rprintf("beta: ") ;
      print(beta) ;
    }
    
    // update R
    for(int i=0; i<N; i++){
      if(State[i]==0){
        tmp1 = alpha[0] ;
        tmp2 = sum(alpha) - alpha[0] ;
        R(i,0) = rtruncbeta(tmp1, tmp2, C[i]/exp(logT[i]), 1.0) ;
        vec_mem_cpy(rDirichlet(1, alpha[Range(1,K-1)])(0,_), tmp_K, K-1) ;
        for(int k=1; k<K; k++)
          R(i,k) = (1.0-R(i,0))*tmp_K[k-1] ;
      }
      else if(State[i] < K-1){
        double V1, V2 ;
        tmp1 = 0.0 ;
        for(int k=0; k<State[i]; k++)
          tmp1 += alpha[k] ;
        tmp2 = sum(alpha) - tmp1 ;
        V1 = rtruncbeta(tmp1, tmp2, 0.0, C[i]/exp(logT[i])) ;
        tmp2 -= alpha[State[i]] ;
        tmp1 = rtruncbeta(alpha[State[i]], tmp2, (C[i]/exp( logT[i])-V1)/(1.0-V1), 1.0) ;
        R(i,State[i]) = (1.0-V1)*tmp1 ;
        V2 = 1.0 - V1 - R(i,State[i]) ;
        vec_mem_cpy(rDirichlet(1, alpha[Range(0,State[i]-1)])(0,_), tmp_K, State[i]) ;
        for(int k=0; k<State[i]; k++)
          R(i,k) = V1*tmp_K[k] ;
        vec_mem_cpy(rDirichlet(1, alpha[Range(State[i]+1,K-1)])(0,_), tmp_K, K-State[i]-1) ;
        for(int k=State[i]+1; k<K; k++)
          R(i,k) = V2*tmp_K[k-State[i]-1] ;
      }
      else if(State[i] < K){
        double V ;
        tmp1 = 0.0 ;
        for(int k=0; k<State[i]; k++)
          tmp1 += alpha[k] ;
        tmp2 = sum(alpha) - tmp1 ;
        V = rtruncbeta(tmp1, tmp2, 0.0, C[i]/exp(logT[i])) ;
        R(i,K-1) = 1.0 - V ;
        vec_mem_cpy(rDirichlet(1, alpha[Range(0,K-2)])(0,_), tmp_K, K-1) ;
        for(int k=0; k<K-1; k++)
          R(i,k) = V*tmp_K[k] ;
      }
      else{
        vec_mem_cpy(rDirichlet(1, alpha)(0,_), tmp_K) ;
        for(int k=0; k<K; k++)
          R(i,k) = tmp_K[k] ;
      }
    }
    
    // update logT
    for(int i=0; i<N; i++){
      mu[i] = b[i] + gSI[i] + Z[i] ;
      
      if(State[i] == 0){
        tmp1 = logC[i] - log(R(i,0)) ;
        tmp2 = R_PosInf ;
        
        logT[i] = rtruncnorm_c(mu[i], sqrt(sigma_sq_eps[i]), tmp1, tmp2) ;
        
      }
      else if(State[i] < K){
        tmp3 = 0.0 ;
        for(int k=0; k<State[i]; k++)
          tmp3 += R(i,k) ;
        tmp2 = logC[i] - log(tmp3) ;
        tmp1 = logC[i] - log(tmp3 + R(i,State[i])) ;
        
        logT[i] = rtruncnorm_c(mu[i], sqrt(sigma_sq_eps[i]), tmp1, tmp2) ;
      }
      else{
        tmp1 = R_NegInf ;
        tmp2 = logC[i] ;
        
        logT[i] = rtruncnorm_c(mu[i], sqrt(sigma_sq_eps[i]), tmp1, tmp2) ;
      }
    }
    
    // update b
    for(int i=0; i<N; i++)
      mu[i] = logT[i] - Z[i] - gSI[i] ;
    
    asMatrix(mu, mu_mtx, true) ;
    update_b(mu_mtx, Sigma_b, sigma_sq_eps, b_mtx) ;
    for(int i=0, ind=0; i<n; i++){
      for(int j=0; j<n_tooth/2; j++, ind++)
        b[ind] = b_mtx(i,j) ;
      for(int j=n_tooth/2; j<n_tooth; j++, ind++)
        b[ind] = b_mtx(i,n_tooth-1-j) ;
    }
    
    // update DP
    for(int i=0; i<N; i++)
      mu[i] = logT[i] - b[i] - gSI[i] ;
    
    labels = update_S(Pi, mu, mm, ss) ;
    pi_out = update_pi(labels, Lmax, gam) ;
    Pi = pi_out["Pi"];
    phi = update_phi(mu, labels, f_mu, f_nu, f_alpha, f_lambda, Lmax);
    
    mm = phi(_,0);
    ss = phi(_,1);
    gam = update_gam(pi_out["Pi_tilde"], Lmax, f_aa, f_bb);
    
    for(int i=0; i<N; i++){
      Z[i] = phi(labels[i]-1,0) ;
      sigma_sq_eps[i] = phi(labels[i]-1,1) ;
    }
    
    eps = mu - Z ;
    SI_obj.change_sigma_sq_eps(sigma_sq_eps) ;
    
    // update g, beta
    for(int i=0; i<N; i++)
      mu[i] = logT[i] - b[i] - Z[i] ;
    SI_obj.change_y(mu) ;
    SI_obj.update_psi() ;
    vec_mem_cpy(SI_obj.psi, psi) ;
    
    SI_obj.update_beta() ;
    vec_mem_cpy(SI_obj.beta, beta) ;
    vec_mem_cpy(SI_obj.beta_tilde, beta_tilde) ;
    product(X, beta, SI) ;
    Db = Dbeta(L, SI);
    product(Db, psi, gSI) ;
    
    
    // update Omega_b/Sigma_b
    AtA(b_mtx, D_Omega) ;
    for(int j=0; j<n_tooth/2; j++)
      for(int jj=0; jj<n_tooth/2; jj++)
        D_Omega(j,jj) += S_Omega(j,jj) ;
    rInvWishart(c_Omega+n, solve(D_Omega), Sigma_b) ;
    
    // update alpha
    alpha_out = Bayesian_gamma_Dirichlet(R, alpha, alpha_stepsize, alpha_a, alpha_lambda) ;
    alpha = alpha_out["alpha"] ;
    alpha_dec = alpha_out["dec"];
    
    
    loglik = calc_loglik(n, n_tooth, mm, ss, Pi, 
                         R, b, gSI, logC, State) ;
    
    // store sample
    if(m >= M_burn && m%thin == 0){
      loglik_mtx((m-M_burn)/thin,_) = loglik ;
      gam_vec[(m-M_burn)/thin] = gam ;
      eps_mtx((m-M_burn)/thin,_) = eps;
      mm_mtx((m-M_burn)/thin,_) = mm;
      ss_mtx((m-M_burn)/thin,_) = ss;
      pi_mtx((m-M_burn)/thin,_) = Pi;
      
      psi_mtx((m-M_burn)/thin,_) = psi ;
      gSI_mtx((m-M_burn)/thin,_) = gSI ;
      beta_mtx((m-M_burn)/thin,_) = beta ;
      beta_tilde_mtx((m-M_burn)/thin,_) = beta_tilde ;
      alpha_mtx((m-M_burn)/thin,_) = alpha ;
      alpha_dec_vec[(m-M_burn)/thin] = alpha_dec ;
      
      Sigma_b_mtx((m-M_burn)/thin,_) = diag(Sigma_b);
      
      for(int j1=0; j1<n_tooth/2; j1++)
        for(int j2=0; j2<n_tooth/2; j2++)
          hat_Sigma_b(j1,j2) += Sigma_b(j1,j2)/(double)M ;
      
      // for(int i=0; i<N; i++){
      //   tmp1 = b[i] + gSI[i] + Z[i] + sqrt(sigma_sq_eps[i])*(rnorm(1)[0]) ;
      //   logT_rep[i] = logT_rep[i] + tmp1/M ;
      //   logT_rep_sq[i] = logT_rep_sq[i] + tmp1*tmp1/M ;
      // }
    }
  }
  
  List res ;
  res["loglik"] = loglik_mtx ;
  res["beta"] = beta_mtx ;
  res["beta_tilde"] = beta_tilde_mtx ;
  res["alpha"] = alpha_mtx ;
  res["alpha_dec"] = alpha_dec_vec ;
  res["psi"] = psi_mtx ;
  
  res["eps"] = eps_mtx ;
  res["gam"] = gam_vec ;
  res["pi"] = pi_mtx ;
  res["mm"] = mm_mtx ;
  res["ss"] = ss_mtx ;
  
  res["Sigma_b"] = Sigma_b_mtx ;
  res["hat_Sigma_b"] = hat_Sigma_b ;
  
  return(res) ;
}