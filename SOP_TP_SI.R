library(gtools); library(Rcpp)
sourceCpp("SPMSM.cpp")

gSI = function(xi, x){
  
  L = length(xi) - 1;
  u = seq(-1, 1, length=L+1)
  
  res = g(x = x, xi = xi, u = u)
  return(res)
  
}

get_sop_SI = function(j, x, State, n_MC, t_vec, alpha, mm, ss, Pi, Sigma_b, xi){
  
  L = length(mm); K = length(alpha); 
  
  s_b = c(Sigma_b, rev(Sigma_b))[j]
  g_x = gSI(xi = xi, x = x)
  
  Z = sample(L, n_MC, TRUE, Pi)
  m_vec = mm[Z]; s_vec = ss[Z]
  
  samR = rdirichlet(n_MC, alpha)
  samT = sapply(1:n_MC, function(i) 
    exp(rnorm(1, mean = g_x + m_vec[i], sd = sqrt(s_b + s_vec[i]))))
  
  calc.SOP = function(k, t){ 
    if(k == 1){
      ind_S = (samT * samR[,1] > t)
    }
    else if(k == (K+1)){
      ind_S = (samT <= t)
    }
    else {
      ind_S = (samT * rowSums(samR[,1:(k-1), drop = FALSE]) <= t) & 
        (t < samT * rowSums(samR[,1:k, drop = FALSE]))
    }
    return(sum(ind_S)/n_MC)
  }
  
  res = sapply(t_vec, function(t) sapply(1:(K+1), calc.SOP, t))
  
  return(res)
}


get_tp_SI = function(j, x, State, n_MC, t0, t_vec, alpha, mm, ss, Pi, Sigma_b, xi){
  
  L = length(mm); K = length(alpha); 
  
  s_b = rep(c(Sigma_b, rev(Sigma_b)))[j]
  g_x = gSI(xi = xi, x = x)
  
  Z = sample(L, n_MC, TRUE, Pi)
  m_vec = mm[Z]; s_vec = ss[Z]
  
  samR = rdirichlet(n_MC, alpha)
  samT = exp((g_x + m_vec) + (rnorm(n_MC) * sqrt(s_b + s_vec)))
  
  
  calc.TP = function(k, t){ 
    if(k == 1){
      ind_S = (samT * samR[,1] > t0)
    }
    else {
      ind_S = (samT * rowSums(samR[,1:(k-1), drop = FALSE]) <= t0) & 
        (t0 < samT * rowSums(samR[,1:k, drop = FALSE]))
    }
    
    if(k == K){
      ind_T = (samT <= t) 
    }
    else{
      ind_T = (samT * rowSums(samR[,1:k, drop = FALSE]) <= t) & 
        (t < samT * rowSums(samR[,1:(k+1), drop = FALSE]))
    }
    
    if(sum(ind_S) == 0){
      x = 0
    }
    else x = sum(ind_T & ind_S)/sum(ind_S)
    return(x)
  }
  
  res = sapply(t_vec, function(t) sapply(1:K, calc.TP, t))
  
  return(res)
  
}