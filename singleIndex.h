#ifndef RCPP_BASIC
#define RCPP_BASIC
#include "RcppBasic.h"
#endif

#ifndef RTRUNCDIST
#define RTRUNCDIST
#include "rTruncDist.h"
#endif

using namespace Rcpp ;

double h(double x){
  // h(x) = (1-|x|)_+
  if(x < -1.0)
    return(0.0) ;
  else if( x < 0.0 )
    return(1.0+x) ;
  else if( x < 1.0 )
    return(1.0-x) ;
  else
    return(0.0) ;
}

double h_l(double x, int l, NumericVector u){
  if(x < u[l]){
    double u1 ;
    if(l==0)
      u1 = 2.0*u[0] - u[1] ;
    else
      u1 = u[l-1] ;
    return( h((x - u[l])/(u[l]-u1)) ) ;
  }
  else{
    double u3 ;
    int L=u.size()-1 ;
    if(l==L)
      u3 = 2.0*u[L] - u[L-1] ;
    else
      u3 = u[l+1] ;
    return( h((x - u[l])/(u3-u[l])) ) ;
  }
}

double psi(double x){
  // \int h(x) dx
  if( x < -1.0 )
    return(0.0) ;
  else if( x < 0.0 )
    return( 0.5*(x+1.0)*(x+1.0) ) ;
  else if(x < 1.0)
    return( 0.5 + 0.5*(2.0-x)*x ) ;
  else
    return(1.0) ;
}

double psi_l(double x, int l, NumericVector u){
  // \int_{-1}^x h_l(t) dt
  
  if(l==0){
    if(x < u[1]){
      return( ( u[1]*(1+x) + 0.5 - 0.5*x*x )/(u[1]+1) ) ;
    }
    else
      return(0.5*(u[1]+1)) ;
  }
  else{
    double u1, tmp ;
    u1 = u[l-1] ;
    
    if(x < u[l]){
      tmp = u[l] - u1 ;
      return( tmp*psi((x-u[l])/tmp) ) ;
    }
    else{
      double u3 ;
      int L=u.size()-1 ;
      if(l==L)
        u3 = 2.0*u[L] - u[L-1] ;
      else
        u3 = u[l+1] ; 
      tmp = u3 - u[l] ;
      return( 0.5*(u[l]-u1) + tmp*(psi((x-u[l])/tmp)-0.5) ) ;
    }
  }
}

double gg(double x, NumericVector xi, NumericVector u){
  // g(x) = \sum_{l in 0:L} xi_(l+1)*psi_l(x)
  if(x > 1.0)
    return(gg(1.0, xi, u)) ; 
  else if(x < -1.0)
    return(0.0) ;
  else{
    double res = 0.0 ;
    for(int l=0; l<xi.size();l++)
      res += xi[l]*psi_l(x,l,u) ;
    return(res) ;
  }
}


//[[Rcpp::export]]
NumericVector psi_vec(double x, NumericVector u){
  int L = u.size();
  NumericVector res(L);
  
  for(int l=0; l<L; l++)
    res[l]= psi_l(x, l, u);
  
  return(res);
}

double log_L(NumericVector y, NumericVector Xbeta, NumericVector xi, 
             double sigma_sq_eps, NumericVector u){
  double res=0.0, tmp ; 
  for(int i=0; i<y.size(); i++){
    tmp = y[i] - gg(Xbeta[i], xi, u) ;
    res += tmp*tmp ;
  }
  res = -0.5*res/sigma_sq_eps ;
  return(res) ;
}

double log_L_varyingSigma(NumericVector y, NumericVector Xbeta, NumericVector xi, 
                          NumericVector sigma_sq_eps, NumericVector u){
  double res=0.0, tmp ;
  for(int i=0; i<y.size(); i++){
    tmp = y[i] - gg(Xbeta[i], xi, u) ;
    res += tmp*tmp/sigma_sq_eps[i] ;
  }
  res = -0.5*res ;
  return (res) ;
}

double log_L_linear(NumericVector y, NumericVector Xbeta,
                    NumericVector sigma_sq_eps){
  double res=0.0, tmp ;
  for(int i=0; i<y.size(); i++){
    tmp = y[i] - Xbeta[i] ;
    res += tmp*tmp/sigma_sq_eps[i] ;
  }
  res = -0.5*res ;
  return (res) ;
}

// calculates B_{L,j}(t)
// L - degree of bernstein polynomials, number of basis functions: L+1
// j - indicates which of the L+1 bernstein polynomials are being evaluated
// t - the value the polynomial is being evaluated at, in [-1,1]
double Btilde(int L, int j, double t) {
  
  // if(t >= 1.0 || t <= -1.0)
  //   return(0.0) ;
  // else
  return((0.5/(L+1)) * R::dbeta((0.5*(t+1)), (j+1), (L-j+1), false)) ;
}


// calculates D_beta matrix
// t is a vector of length n, t in [-1,1]
// L is the degree of the bernstein polynomials
//[[Rcpp::export]]
NumericMatrix Dbeta(int L, NumericVector t) {
  
  int n = t.size();
  NumericMatrix BigB(n, L), A(L, L), D_mat(n, L);
  
  for(int i=0; i<n; i++)
    for(int j=0; j<L; j++)
      BigB(i,j) = Btilde(L, j+1, t[i]) ;
  // rows are over t, columns are over 0:L, values[i,j] is Btilde(L,j,t[i])
  
  for(int i=0; i<L; i++)
    for(int j=0; j<L; j++){
      if(i >= j) {
        A(i, j) = 1;
      }
      else {
        A(i, j) = 0;
      }
    }
    product(BigB, A, D_mat);
  return(D_mat);
}


double log_L_BP(NumericVector y, NumericVector Xbeta, 
                NumericVector psi, NumericVector sigma_sq_eps){
  int L, N;
  L = psi.size() ;
  N = y.size();
  double res=0.0, tmp ;
  NumericVector Db_psi(N) ;
  NumericMatrix Db;
  
  Db = Dbeta(L, Xbeta);
  product(Db, psi, Db_psi) ;
  
  for(int i=0; i<N; i++){
    tmp = y[i] - Db_psi[i] ;
    res += tmp*tmp/sigma_sq_eps[i] ;
  }
  res = -0.5*res ;
  return (res) ;
}


class singleIndex_varyingSigma{
  
  //// y[i] = g( sum(X(i,_)*beta) ) + sigma_sq_eps[i]*eps[i]
  //// eps[i] ~ N(0, sigma_sq_eps[i])
  //// n_HMC: # of HMC iteration in truncated normal sampling
  
private:
  int n, p, L ;
  double sigma_sq_beta, a_xi, b_xi, nu_matern, l_matern, eta_xi ;
  // NumericMatrix Psi, PsiPsi, Omega_xi ;
  
  // Following variables should be defined in sampling methods
  // They are defined here for speed-up by avoiding repetitive memory allocation
  NumericVector nu_xi, tmp_L, nu, beta_new, beta_tilde_new, Xbeta_new ;
  NumericMatrix Kinv; //B_xi, B_xi_inv, ff ;
  
public:
  
  NumericVector beta_tilde;
  NumericVector y, u, Xbeta, eps, xi, beta, sigma_sq_eps;
  NumericMatrix X;
  // List xi_out ;
  singleIndex_varyingSigma(NumericVector y_obs, NumericMatrix X_obs, NumericVector _u, 
                           double _nu_matern, double _l_matern, NumericMatrix _Kinv, double _eta_xi,
                           NumericVector init_beta, NumericVector init_xi,
                           NumericVector init_sigma_sq_eps, double _sigma_sq_beta=10.0, 
                           double _a_xi=3, double _b_xi=10){
    
    n = y_obs.size() ;
    p = X_obs.ncol() ;
    L = init_xi.size() - 1 ;
    
    NumericVector _y(n), _beta(p), _beta_tilde(p), _xi(L+1), _Xbeta(n), _eps(n), _sigma_sq_eps(n) ;
    NumericMatrix _X(n,p) ; //_Psi(n,L+1), _PsiPsi(L+1,L+1), _Omega_xi(L+1, L+1) ;
    
    NumericVector _nu_xi(L+1), _tmp_L(L+1), _nu(p), _beta_new(p), _beta_tilde_new(p), _Xbeta_new(n) ;
    // NumericMatrix _B_xi(L+1,L+1), _B_xi_inv(L+1,L+1), _ff(L+1,L+1) ;
    
    
    y = _y ;
    X = _X ;
    u = _u ;
    beta = _beta ;
    beta_tilde = _beta_tilde ;
    xi = _xi ;
    Xbeta = _Xbeta ;
    eps = _eps ;
    // Psi = _Psi ;
    // PsiPsi = _PsiPsi ;
    // Omega_xi = _Omega_xi ;
    a_xi = _a_xi ;
    b_xi = _b_xi;
    nu_matern = _nu_matern ;
    l_matern = _l_matern ;
    eta_xi = _eta_xi ;
    Kinv = _Kinv ;
    
    sigma_sq_eps = _sigma_sq_eps ;
    sigma_sq_beta = _sigma_sq_beta ;
    
    nu_xi = _nu_xi ;
    tmp_L = _tmp_L ;
    nu = _nu ;
    beta_new = _beta_new ;
    beta_tilde_new = _beta_tilde_new ;
    Xbeta_new = _Xbeta_new ; 
    
    // B_xi = _B_xi ;
    // B_xi_inv = _B_xi_inv ;
    // ff = _ff ;
    
    // initial setting
    vec_mem_cpy(y_obs, y) ;
    mtx_mem_cpy(X_obs, X) ;
    vec_mem_cpy(init_beta, beta) ;
    vec_mem_cpy(init_sigma_sq_eps, sigma_sq_eps) ;
    beta_tilde = sqrt(p*sigma_sq_beta)*beta ;
    // beta_tilde = rnorm(p, 0, sqrt(sigma_sq_beta)) ;
    vec_mem_cpy(init_xi, xi) ;
    product(X, beta, Xbeta) ;
    
  }
  
  // NumericVector y, xi, beta, sigma_sq_eps ;
  // NumericMatrix X ;
  void update_xi(){
    Function update_basis_coeff("update_basis_coeff"); 
    
    xi = update_basis_coeff(y, Xbeta, u, sigma_sq_eps, nu_matern, l_matern, Kinv, eta_xi, xi,
                            a_xi, b_xi);
    // xi = xi_out["xi"] ;
  }
  void update_beta(){
    double uu, log_y, theta, theta_min, theta_max, log_l_new ;
    uu = runif(1)[0] ;
    nu = rnorm(p, 0.0, sqrt(sigma_sq_beta)) ;
    log_y = log_L_varyingSigma(y, Xbeta, xi, sigma_sq_eps, u) + log(uu) ;
    // Rprintf("log L beta evaluated \n");
    theta = 2.0*M_PI*runif(1)[0] ;
    theta_min = theta - 2.0*M_PI ;
    theta_max = theta ;
    beta_tilde_new = cos(theta)*beta_tilde + sin(theta)*nu ;
    beta_new = beta_tilde_new / norm(beta_tilde_new) ;
    product(X, beta_new, Xbeta_new) ;
    // Rprintf("Beta \n");
    // Rcpp::print(beta);
    
    // Rprintf("Begin ESS \n") ;
    for(int i=0; ; i++){
      // Rprintf("Iteration :%d\n", i+1) ;
      log_l_new = log_L_varyingSigma(y, Xbeta_new, xi, sigma_sq_eps, u) ; 
      // Rprintf("Beta \n");
      // Rcpp::print(beta_new);
      // Rprintf("log L beta evaluated \n");
      
      double diff = log_l_new - log_y ;
      // Rprintf("Diff %d\n", diff);
      
      if(diff> 0){
        // Rprintf("Number of iterations in ESS = %d\n", i+1) ;
        break ;
      }
      
      if(theta < 0)
        theta_min = theta ;
      else
        theta_max = theta ;
      
      theta = (theta_max-theta_min)*runif(1)[0] + theta_min ;
      beta_tilde_new = cos(theta)*beta_tilde + sin(theta)*nu ;
      beta_new = beta_tilde_new / norm(beta_tilde_new) ;
      product(X, beta_new, Xbeta_new) ;
    }
    // Rprintf("End ESS \n") ;
    
    vec_mem_cpy(beta_new, beta) ;
    vec_mem_cpy(beta_tilde_new, beta_tilde) ;
    vec_mem_cpy(Xbeta_new, Xbeta) ;
  } 
  void change_sigma_sq_eps(NumericVector sigma_sq_eps_new){
    vec_mem_cpy(sigma_sq_eps_new, sigma_sq_eps) ;
  }
  void change_y(NumericVector y_new){
    vec_mem_cpy(y_new, y) ;
  }
  void print(){
    Rprintf("Results of Single Index::print()\n") ;
    Rprintf("beta: ") ;
    Rcpp::print(beta) ;
  }
  
} ;

class singleIndex_BP{
  
  //// y[i] = g( sum(X(i,_)*beta) ) + sigma_sq_eps[i]*eps[i]
  //// eps[i] ~ N(0, sigma_sq_eps[i])
  //// n_HMC: # of HMC iteration in truncated normal sampling
  
private:
  int n, p, L ;
  double sigma_sq_beta, sigma_sq_psi, a_psi, b_psi ;
  NumericVector Xbeta, eps ;
  NumericMatrix D, DD, Omega_psi, Db ;
  
  // Following variables should be defined in sampling methods
  // They are defined here for speed-up by avoiding repetitive memory allocation
  NumericVector nu_psi, tmp_L, nu, beta_new, beta_tilde_new, Xbeta_new ;
  NumericMatrix B_psi, B_psi_inv, ff ;
  
public:
  NumericVector beta_tilde;
  singleIndex_BP(NumericVector y_obs, NumericMatrix X_obs,
                 NumericVector init_beta, NumericVector init_psi, NumericVector init_sigma_sq_eps, 
                 double _sigma_sq_beta=1.0, double _a_psi=0.1, double _b_psi=1000){
    
    n = y_obs.size() ;
    p = X_obs.ncol() ;
    L = init_psi.size() ;
    
    NumericVector _y(n), _beta(p), _beta_tilde(p), _psi(L), _Xbeta(n), _eps(n), _sigma_sq_eps(n) ;
    NumericMatrix _X(n,p), _D(n,L), _DD(L,L), _Omega_psi(L, L) ;
    
    NumericVector _nu_psi(L), _tmp_L(L), _nu(p), _beta_new(p), _beta_tilde_new(p), _Xbeta_new(n) ;
    NumericMatrix _B_psi(L, L), _B_psi_inv(L, L), _ff(L,L) ;
    
    
    y = _y ;
    X = _X ;
    // u = _u ;
    beta = _beta ;
    beta_tilde = _beta_tilde ;
    psi = _psi ;
    Xbeta = _Xbeta ;
    eps = _eps ;
    D = _D ;
    DD = _DD ;
    Omega_psi = _Omega_psi ;
    a_psi = _a_psi ;
    b_psi = _b_psi;
    
    sigma_sq_eps = _sigma_sq_eps ;
    sigma_sq_beta = _sigma_sq_beta ;
    
    nu_psi = _nu_psi ;
    tmp_L = _tmp_L ;
    nu = _nu ;
    beta_new = _beta_new ;
    beta_tilde_new = _beta_tilde_new ;
    Xbeta_new = _Xbeta_new ; 
    B_psi = _B_psi ;
    B_psi_inv = _B_psi_inv ;
    ff = _ff ;
    
    // initial setting
    vec_mem_cpy(y_obs, y) ;
    mtx_mem_cpy(X_obs, X) ;
    vec_mem_cpy(init_beta, beta) ;
    vec_mem_cpy(init_sigma_sq_eps, sigma_sq_eps) ;
    beta_tilde = sqrt(p*sigma_sq_beta)*beta ;
    vec_mem_cpy(init_psi, psi) ;
    product(X, beta, Xbeta) ;
    Db = Dbeta(L, Xbeta) ;
    
    for(int i=0; i<n; i++)
      for(int l=0; l<L; l++)
        D(i,l) = Db(i, l)/sqrt(sigma_sq_eps[i]) ;
    AtA(D, DD) ;
    // Omega_psi = diag_c((L+1), (1.0/sigma_sq_psi)) ;
    // Omega_psi = (1/sigma_sq_psi)* I_((L+1))
    // solve(Sigma_psi, Omega_psi) ;
    
    for(int l=0; l<L; l++)
      ff(l,l) = 1.0 ;
  }
  
  NumericVector y, psi, beta, sigma_sq_eps ;
  NumericMatrix X ;
  void update_psi(int n_HMC=1){
    
    sigma_sq_psi = 1.0/rgamma(1, a_psi + 0.5*(L+1), 1.0/(b_psi + 0.5*innerProduct(psi,psi)))[0] ;
    
    Omega_psi = diag_c((L+1), (1.0/sigma_sq_psi)) ;
    sum(DD, Omega_psi, B_psi_inv) ;
    solve(B_psi_inv, B_psi) ;
    Atx(D, y/sqrt(sigma_sq_eps), tmp_L) ;
    product(B_psi, tmp_L, nu_psi) ;
    psi = rtmvnormHMC(n_HMC, nu_psi, B_psi, psi, ff, rep(0.0, L), 0)(n_HMC-1,_) ;
  }
  void update_beta(){
    double uu, log_y, theta, theta_min, theta_max ;
    uu = runif(1)[0] ;
    nu = rnorm(p, 0.0, sqrt(sigma_sq_beta)) ;
    log_y = log_L_BP(y, Xbeta, psi, sigma_sq_eps) + log(uu) ;
    theta = 2.0*M_PI*runif(1)[0] ;
    theta_min = theta - 2.0*M_PI ;
    theta_max = theta ;
    beta_tilde_new = cos(theta)*beta_tilde + sin(theta)*nu ;
    beta_new = beta_tilde_new / norm(beta_tilde_new) ;
    product(X, beta_new, Xbeta_new) ;
    for(int i=0; ; i++){
      if(log_L_BP(y, Xbeta_new, psi, sigma_sq_eps) > log_y){
        //Rprintf("Number of iterations in ESS = %d\n", i+1) ;
        break ;
      }
      if(theta < 0)
        theta_min = theta ;
      else
        theta_max = theta ;
      theta = (theta_max-theta_min)*runif(1)[0] + theta_min ;
      beta_tilde_new = cos(theta)*beta_tilde + sin(theta)*nu ;
      beta_new = beta_tilde_new / norm(beta_tilde_new) ;
      product(X, beta_new, Xbeta_new) ;
    }
    
    vec_mem_cpy(beta_new, beta) ;
    vec_mem_cpy(beta_tilde_new, beta_tilde) ;
    vec_mem_cpy(Xbeta_new, Xbeta) ;
    
    Db = Dbeta(L, Xbeta) ;
    
    for(int i=0; i<n; i++)
      for(int l=0; l<L; l++)
        D(i,l) = Db(i, l)/sqrt(sigma_sq_eps[i]) ;
    AtA(D, DD) ;
  } 
  void change_sigma_sq_eps(NumericVector sigma_sq_eps_new){
    vec_mem_cpy(sigma_sq_eps_new, sigma_sq_eps) ;
  }
  void change_y(NumericVector y_new){
    vec_mem_cpy(y_new, y) ;
  }
} ;
