rm(list=ls())
library(gtools); library(MASS); library(Rcpp) ; library(truncnorm)
source("data_generate.R")
source("circulantembedding.R")
sourceCpp("SPMSM.cpp")

# poly link function
true.g = function(x, a = 1){
  n = length(x)
  res = rep(0, n)
  for(i in 1:n){
    if(x[i]<=0){
      res[i] = 0
    }
    else res[i] = x[i]^2
  }
  return(a*res)
}
true.alpha = c(3, 4, 2)
tmp = c(-1, 1, -1); true.beta=tmp/norm(tmp,"2")
mu_eps = c(-0.5, 0, 0.5)

n = 500; m =10; p = 3 ; K=3;  N=n*m
xnames = paste0("x", 1:p)

L=30
u = seq(-1, 1, length=L+1)
nu_matern = 0.75
l_matern=l_est(nu_matern,c(u[1],u[length(u)]),0.05)
# prior covariance K:
K_mat = covmat(u,nu_matern,l_matern)
# prior precision:
K_inv = tinv(K_mat)


grid.g = seq(-1,1,length=100)

M = 2000; M_burn = 5000 ; M_thin = 1

#############################################################################

library(doParallel)
nReplicate = 100

mycombine = function(...){
  
  list.all = list(...)
  list.out = vector(mode = 'list', length = length(list.all[[1]]))
  for(r in 1:length(list.all)){
    
    Z = 10*4+1
    for(k in 1:Z){
      list.out[[k]] = rbind(list.out[[k]], list.all[[r]][[k]])
    }
  }
  
  list.out
}

doParallel::registerDoParallel(cores = detectCores()-1)

####################################################################################

out = foreach::foreach(r = 1:nReplicate, .combine = 'mycombine', .multicombine = T, 
                       .errorhandling = 'remove') %dopar% {
  
  data_all = data_misspecification(n = n, m=m, p=p, true.alpha = true.alpha, 
                           true.beta = true.beta, C.param = c(0.5, 0.2))
  obs.full = data_all$obs.full
  const.g = data_all$const.g
  
  grid.x = seq(min(obs.full[,"eps"])-1, max(obs.full[,"eps"])+1, length=10) # grid pts evaluating f(x)
  
  init.xi = pmax(0,mvrnorm(1,rep(0,L+1),K_mat))
  init.psi = rep(1,L)
  
  tmp=rnorm(p); init.beta=tmp/norm(tmp,"2")
  init.alpha = rep(1,K)
  
  init.Sigma.b=matrix(0.99,nrow=m/2,ncol=m/2); 
  diag(init.Sigma.b)=1; init.Sigma.b=0.1*init.Sigma.b
  
  init.logT = rep(0,N)
  for(i in 1:N){
    if(obs.full[i,"State"] == K){
      init.logT[i] = log(obs.full[i,"C"]*runif(1))
    }
    else init.logT[i] = log(obs.full[i,"C"] + 0.1*runif(1))
  }
  
  rho = 0.9
  tmp=diag(rep(2,m/2)); tmp[1,1]=tmp[m/2,m/2]=1
  for(j in 1:(m/2-1))
    tmp[j,j+1] = tmp[j+1,j] = -rho ;
  E.Sigma.b = tmp
  
  E.Sigma.b = solve(tmp)
  S.Omega = E.Sigma.b
  c.Omega = m/2 + 2
  
  res = GP_MSM_c(id = obs.full[,"id"], tooth = obs.full[,"tooth"], C  = obs.full[,"C"], 
                State = obs.full[,"State"], X = obs.full[, xnames], 
                u = u, S_Omega = S.Omega, c_Omega = c.Omega, init_xi = init.xi, 
                init_beta = init.beta, init_alpha = init.alpha, init_Sigma_b = init.Sigma.b, 
                init_logT = init.logT, 
                Kinv = K_inv, l_matern = l_matern, nu_matern = nu_matern, eta_xi = 100,
                a_xi = 0.01, b_xi = 0.01, sigma_sq_beta = 100,  
                alpha_a = 1, alpha_lambda = 0.1, alpha_stepsize = 0.08,
                f_mu = 0, f_nu = 1, f_lambda = 0.1, f_alpha = 2.001,
                Lmax = 10, M_mcmc = M, M_burn = M_burn, thin = M_thin) 
  
  resBP = BP_MSM_c(id = obs.full[,"id"], tooth = obs.full[,"tooth"], C  = obs.full[,"C"], 
                     State = obs.full[,"State"], X = obs.full[, xnames], 
                     S_Omega = S.Omega, c_Omega = c.Omega, init_psi = init.psi, 
                     init_beta = init.beta, init_alpha = init.alpha, 
                     init_Sigma_b = init.Sigma.b, init_logT = init.logT, 
                     sigma_sq_beta = 100, a_psi = 0.01, b_psi = 0.01, 
                     alpha_a = 1, alpha_lambda = 0.1, alpha_stepsize = 0.08,
                     f_mu = 0, f_nu = 1, f_lambda = 0.1, f_alpha = 2.001, 
                     Lmax = 10, M_mcmc = M, M_burn = M_burn, thin = M_thin) 
  
  res_norm = GP_MSM_c(id = obs.full[,"id"], tooth = obs.full[,"tooth"], C  = obs.full[,"C"], 
                     State = obs.full[,"State"], X = obs.full[, xnames], 
                     u = u, S_Omega = S.Omega, c_Omega = c.Omega, init_xi = init.xi, 
                     init_beta = init.beta, init_alpha = init.alpha, init_Sigma_b = init.Sigma.b, 
                     init_logT = init.logT, 
                     Kinv = K_inv, l_matern = l_matern, nu_matern = nu_matern, eta_xi = 100,
                     a_xi = 0.01, b_xi = 0.01, sigma_sq_beta = 100,  
                     alpha_a = 1, alpha_lambda = 0.1, alpha_stepsize = 0.08,
                     f_mu = 0, f_nu = 1, f_lambda = 0.1, f_alpha = 2.001,
                     Lmax = 1, M_mcmc = M, M_burn = M_burn, thin = M_thin) 
  
  resBP_norm = BP_MSM_c(id = obs.full[,"id"], tooth = obs.full[,"tooth"], C  = obs.full[,"C"], 
                          State = obs.full[,"State"], X = obs.full[, xnames], 
                          S_Omega = S.Omega, c_Omega = c.Omega, init_psi = init.psi, 
                          init_beta = init.beta, init_alpha = init.alpha, 
                          init_Sigma_b = init.Sigma.b, init_logT = init.logT, 
                          sigma_sq_beta = 100, a_psi = 0.01, b_psi = 0.01, 
                          alpha_a = 1, alpha_lambda = 0.1, alpha_stepsize = 0.08,
                          f_mu = 0, f_nu = 1, f_lambda = 0.1, f_alpha = 2.001, 
                          Lmax = 1, M_mcmc = M, M_burn = M_burn, thin = M_thin) 
  
  ####################################################################
  # store samples
  # GP_DP
  sam.xi = res$xi
  sam.beta = res$beta
  
  # BP_DP
  sam.beta.BP = resBP$beta
  sam.psi.BP = resBP$psi
  
  # GP_norm
  sam.xi.norm = res_norm$xi
  sam.beta.norm = res_norm$beta
  
  # BP_norm
  sam.beta.BP.norm = resBP_norm$beta
  sam.psi.BP.norm = resBP_norm$psi
  
  #####################################################################
  # beta
  # GP_DP
  hat.beta = colMeans(sam.beta)
  hat.beta = hat.beta / norm(hat.beta, "2")
  CI.beta = apply(sam.beta, 2, quantile, probs=c(0.025,0.975))
  cov.beta = sapply(1:p, function(k) 
    1*(true.beta[k] >= CI.beta[1, k] && true.beta[k] <= CI.beta[2, k]))
  MSE.beta = sapply(1:p, function(k) (true.beta[k] - hat.beta[k])^2)
  RB.beta = sapply(1:p, function(k) (hat.beta[k] - true.beta[k])/true.beta[k])
  
  # BP_DP
  hat.beta.BP = colMeans(sam.beta.BP)
  hat.beta.BP = hat.beta.BP / norm(hat.beta.BP, "2")
  CI.beta.BP = apply(sam.beta.BP, 2, quantile, probs=c(0.025,0.975))
  cov.beta.BP = sapply(1:p, function(k) 
    1*(true.beta[k] >= CI.beta.BP[1, k] && true.beta[k] <= CI.beta.BP[2, k]))
  MSE.beta.BP = sapply(1:p, function(k) (true.beta[k] - hat.beta.BP[k])^2)
  RB.beta.BP = sapply(1:p, function(k) (hat.beta.BP[k] - true.beta[k])/true.beta[k])
  
  # GP_norm
  hat.beta.norm = colMeans(sam.beta.norm)
  hat.beta.norm = hat.beta.norm / norm(hat.beta.norm, "2")
  CI.beta.norm = apply(sam.beta.norm, 2, quantile, probs=c(0.025,0.975))
  cov.beta.norm = sapply(1:p, function(k) 
    1*(true.beta[k] >= CI.beta.norm[1, k] && true.beta[k] <= CI.beta.norm[2, k]))
  MSE.beta.norm = sapply(1:p, function(k) (true.beta[k] - hat.beta.norm[k])^2)
  RB.beta.norm = sapply(1:p, function(k) (hat.beta.norm[k] - true.beta[k])/true.beta[k])
  
  # BP_norm
  hat.beta.BP.norm = colMeans(sam.beta.BP.norm)
  hat.beta.BP.norm = hat.beta.BP.norm / norm(hat.beta.BP.norm, "2")
  CI.beta.BP.norm = apply(sam.beta.BP.norm, 2, quantile, probs=c(0.025,0.975))
  cov.beta.BP.norm = sapply(1:p, function(k) 
    1*(true.beta[k] >= CI.beta.BP.norm[1, k] && true.beta[k] <= CI.beta.BP.norm[2, k]))
  MSE.beta.BP.norm = sapply(1:p, function(k) (true.beta[k] - hat.beta.BP.norm[k])^2)
  RB.beta.BP.norm = sapply(1:p, function(k) (hat.beta.BP.norm[k] - true.beta[k])/true.beta[k])
  
  #####################################################################
  # g
  g.true = rep(0,length(grid.g))
  for(i in 1:length(grid.g))
    g.true[i] = true.g(grid.g[i], a = const.g)
  D_mat = Dbeta(L, grid.g)
  
  Xbeta.true = obs.full[,"SI"]
  g.Xbeta.true = obs.full[, "gSI"]
  
  # GP_DP
  hat.xi = colMeans(sam.xi)
  hat.g = rep(NA, length(grid.g))
  for(i in 1:length(grid.g)){
    hat.g[i] = g(grid.g[i], hat.xi, u)
  }
  MISE.g = mean((g.true - hat.g)^2)
  
  hat.Xbeta = c(obs.full[, xnames] %*% hat.beta)
  hat.g.Xbeta = rep(0,N)
  for(i in 1:N){
    hat.g.Xbeta[i] = g(hat.Xbeta[i], hat.xi, u)
  }
  MISE.g.Xbeta = mean((g.Xbeta.true - hat.g.Xbeta)^2)
  
  
  # BP_DP
  g.BP = t(sapply(1:M, function(m) D_mat%*%sam.psi.BP[m, ]))
  hat.g.BP = colMeans(g.BP)
  MISE.g.BP = mean((g.true - hat.g.BP)^2)
  
  hat.Xbeta.BP = c(obs.full[, xnames] %*% hat.beta.BP)
  D_mat_Xbeta.BP = Dbeta(L, hat.Xbeta.BP)
  g.BP.Xbeta = t(sapply(1:M, function(m) D_mat_Xbeta.BP%*%sam.psi.BP[m, ]))
  hat.g.BP.Xbeta = colMeans(g.BP.Xbeta)
  MISE.g.BP.Xbeta = mean((g.Xbeta.true - hat.g.BP.Xbeta)^2)
  
  
  # GP_norm
  hat.xi.norm = colMeans(sam.xi.norm)
  hat.g.norm = rep(NA, length(grid.g))
  for(i in 1:length(grid.g)){
    hat.g.norm[i] = g(grid.g[i], hat.xi.norm, u)
  }
  MISE.g.norm = mean((g.true - hat.g.norm)^2)
  
  hat.Xbeta.norm = c(obs.full[, xnames] %*% hat.beta.norm)
  hat.g.norm.Xbeta = rep(0,N)
  for(i in 1:N){
    hat.g.norm.Xbeta[i] = g(hat.Xbeta.norm[i], hat.xi, u)
  }
  MISE.g.norm.Xbeta = mean((g.Xbeta.true - hat.g.norm.Xbeta)^2)
  
  
  # BP_norm
  g.BP.norm = t(sapply(1:M, function(m) D_mat%*%sam.psi.BP.norm[m, ]))
  hat.g.BP.norm = colMeans(g.BP.norm)
  MISE.g.BP.norm = mean((g.true - hat.g.BP.norm)^2)
  
  hat.Xbeta.BP.norm = c(obs.full[, xnames] %*% hat.beta.BP.norm)
  D_mat_Xbeta.BP.norm = Dbeta(L, hat.Xbeta.BP.norm)
  g.BP.norm.Xbeta = t(sapply(1:M, function(m) D_mat_Xbeta.BP.norm%*%sam.psi.BP.norm[m, ]))
  hat.g.BP.norm.Xbeta = colMeans(g.BP.norm.Xbeta)
  MISE.g.BP.norm.Xbeta = mean((g.Xbeta.true - hat.g.BP.norm.Xbeta)^2)
  
  
  save.list = list(hat.beta, hat.beta.BP, hat.beta.norm, hat.beta.BP.norm,
                   CI.beta[1,], CI.beta.BP[1,], CI.beta.norm[1,], CI.beta.BP.norm[1,],
                   CI.beta[2,], CI.beta.BP[2,], CI.beta.norm[2,], CI.beta.BP.norm[2,],
                   cov.beta, cov.beta.BP, cov.beta.norm, cov.beta.BP.norm,
                   MSE.beta, MSE.beta.BP, MSE.beta.norm, MSE.beta.BP.norm,
                   RB.beta, RB.beta.BP, RB.beta.norm, RB.beta.BP.norm,
                   hat.g, hat.g.BP, hat.g.norm, hat.g.BP.norm,
                   MISE.g, MISE.g.BP, MISE.g.norm, MISE.g.BP.norm, 
                   hat.g.Xbeta, hat.g.BP.Xbeta, hat.g.norm.Xbeta, hat.g.BP.norm.Xbeta,
                   MISE.g.Xbeta, MISE.g.BP.Xbeta, MISE.g.norm.Xbeta, MISE.g.BP.norm.Xbeta,
                   const.g)
  
  ##########################################################################
  
  c(save.list)
}

names(out) = c("hat.beta", "hat.beta.BP", "hat.beta.norm", "hat.beta.BP.norm",
               "CI.beta.1", "CI.beta.BP.1", "CI.beta.norm.1", "CI.beta.BP.norm.1",
               "CI.beta.2", "CI.beta.BP.2", "CI.beta.norm.2", "CI.beta.BP.norm.2",
               "cov.beta", "cov.beta.BP", "cov.beta.norm", "cov.beta.BP.norm",
               "MSE.beta", "MSE.beta.BP", "MSE.beta.norm", "MSE.beta.BP.norm",
               "RB.beta", "RB.beta.BP", "RB.beta.norm", "RB.beta.BP.norm",
               "hat.g", "hat.g.BP", "hat.g.norm", "hat.g.BP.norm",
               "MISE.g", "MISE.g.BP", "MISE.g.norm", "MISE.g.BP.norm", 
               "hat.g.Xbeta", "hat.g.BP.Xbeta", "hat.g.norm.Xbeta", "hat.g.BP.norm.Xbeta",
               "MISE.g.Xbeta", "MISE.g.BP.Xbeta", "MISE.g.norm.Xbeta", "MISE.g.BP.norm.Xbeta",
               "const.g")

out4 = c(list(n=n, m=m, p=p, true.beta = true.beta), out)
filename.rdata = paste0("out_poly_misspec_n", n, ".Rdata") 
save(out4, file = filename.rdata)
