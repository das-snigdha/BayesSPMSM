# R wrapper for the Cpp function to implement the MCMC sampler
library(gtools); library(MASS); library(Rcpp);

source("data_generate.R")
source("circulantembedding.R")
sourceCpp("SPMSM.cpp")

GP_MSM = function(id, tooth, C, State, X, sigma.sq.beta, S.Omega, c.Omega,
                  u = u, l.matern, nu.matern, eta.xi, a.xi, b.xi, L, 
                  alpha.a, alpha.lambda, alpha.tune, 
                  f.mu, f.nu, f.lambda, f.alpha, L.max, M = 1000, M.burn = 500, M.thin = 1){
  
  n = max(id); p = ncol(X); K = max(State) ; n_tooth = max(tooth); N = n*n_tooth
  
  # initialize parameters
  
  # basis coefficients
  K_mat = covmat(u, nu.matern, l.matern); K_inv = tinv(K_mat)
  init.xi = pmax(0, mvrnorm(1, rep(0, L+1), K_mat))
  
  # Regression parameter
  tmp = rnorm(p); init.beta = tmp/norm(tmp,"2")
  
  # Dirichlet parameters for the relative disease increment times
  init.alpha = rep(1,K)
  
  # Covariance matrix of the random effects
  init.Sigma.b = matrix(0.99, nrow=m/2, ncol=m/2); 
  diag(init.Sigma.b) = 1; init.Sigma.b = 0.1*init.Sigma.b
  
  # Time to a missing tooth
  init.logT = rep(0,N)
  for(i in 1:N){
    if(obs.full[i, "State"] == K){
      init.logT[i] = log(obs.full[i, "C"]*runif(1))
    }
    else init.logT[i] = log(obs.full[i, "C"] + 0.1*runif(1))
  }
  
  # Call the Cpp function to run the MCMC chain
  out = GP_MSM_c(id = id, tooth = tooth, C = C, State = State, X = X, 
                 u = u, S_Omega = S.Omega, c_Omega = c.Omega, init_xi = init.xi, 
                 init_beta = init.beta, init_alpha = init.alpha, 
                 init_Sigma_b = init.Sigma.b, init_logT = init.logT, 
                 Kinv = K_inv, l_matern = l.matern, nu_matern = nu.matern, eta_xi = eta.xi,
                 a_xi = a.xi, b_xi = b.xi, sigma_sq_beta = sigma.sq.beta,  
                 alpha_a = alpha.a, alpha_lambda = alpha.lambda, alpha_stepsize = alpha.tune,
                 f_mu = f.mu, f_nu = f.nu, f_lambda = f.lambda, f_alpha = f.alpha,
                 Lmax = L.max, M_mcmc = M, M_burn = M.burn, thin = M.thin)
  
  return(out)
}