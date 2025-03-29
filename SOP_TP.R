library(gtools); library(Rcpp)
sourceCpp("SPMSM.cpp")

scale.new.covariate = function(x.new, x.data, wts){
  
  center.x = rbind(attr(x.data, "scaled:center"),attr(x.data, "scaled:center"))
  scale.x = rbind(attr(x.data, "scaled:scale"),attr(x.data, "scaled:scale"))
  
  x.new.scaled = x.new - center.x
  x.new.scaled = x.new.scaled / scale.x
  x.new.scaled = x.new.scaled/wts
  
  return(x.new.scaled) 
}

gSI_beta = function(xi, x, Beta){
  
  L = length(xi) - 1;
  u = seq(-1, 1, length=L+1)
  
  xb = sum(x*Beta)
  
  res = g(x = xb, xi = xi, u = u)
  return(res)
  
}

get_sop = function(j, x, State, n_MC, t_vec, alpha, mm, ss, Pi, Sigma_b, xi, Beta){
  
  L = length(mm); K = length(alpha); 
  
  s_b = c(Sigma_b, rev(Sigma_b))[j]
  g_x = gSI_beta(xi = xi, x = x, Beta = Beta)
  
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


get_tp = function(j, x, State, n_MC, t0, t_vec, alpha, mm, ss, Pi, Sigma_b, xi, Beta){
  
  L = length(mm); K = length(alpha); 
  
  s_b = rep(c(Sigma_b, rev(Sigma_b)))[j]
  g_x = gSI_beta(xi = xi, x = x, Beta = Beta)
  
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

est_SOP_TP = function(M_mc, out, x.new.scaled, tooth_num, t0, t_vec, t_vec_tp){

  M_mcmc = length(out$alpha_dec)
  SOP = vector(mode = "list", length = nrow(x.new))
  TP = vector(mode = "list", length = nrow(x.new))
  
  for(k in 1:nrow(x.new)){
    print(paste("k :", k))
    
    SOP[[k]] = lapply(1:M_mcmc, function(i) 
      get_sop(j = tooth_num[k], x =  x.new.scaled[k, ], State = obs.full[, "State"],
              n_MC = M_mc, t_vec = t_vec, alpha = out$alpha[i, ], 
              mm = out$mm[i, ], ss = out$ss[i, ], Pi = out$pi[i, ],
              Sigma_b = out$Sigma_b[i, ], xi = out$xi[i, ],
              Beta = out$beta[i, ]) )
    
    TP[[k]] = lapply(1:M_mcmc, function(i) 
      get_tp(j = tooth_num[k], x =  x.new.scaled[k, ], State = obs.full[, "State"],
             n_MC = M_mc, t0 = t0, t_vec = t_vec_tp, alpha = out$alpha[i, ], 
             mm = out$mm[i, ], ss = out$ss[i, ], Pi = out$pi[i, ],
             Sigma_b = out$Sigma_b[i, ], xi = out$xi[i, ],
             Beta = out$beta[i, ]) )
  }
  
  SOP1.state1 = t(sapply(SOP[[1]], function(x) x[1,]))
  SOP1.state2 = t(sapply(SOP[[1]], function(x) x[2,]))
  SOP1.state3 = t(sapply(SOP[[1]], function(x) x[3,]))
  SOP1.state4 = t(sapply(SOP[[1]], function(x) x[4,]))
  
  hat.SOP1.state1 = colMeans(SOP1.state1)
  band.SOP1.state1 = apply(SOP1.state1, 2, quantile, probs=c(0.025,0.975))
  hat.SOP1.state2 = colMeans(SOP1.state2)
  band.SOP1.state2 = apply(SOP1.state2, 2, quantile, probs=c(0.025,0.975))
  hat.SOP1.state3 = colMeans(SOP1.state3)
  band.SOP1.state3 = apply(SOP1.state3, 2, quantile, probs=c(0.025,0.975))
  hat.SOP1.state4 = colMeans(SOP1.state4)
  band.SOP1.state4 = apply(SOP1.state4, 2, quantile, probs=c(0.025,0.975))
  
  SOP2.state1 = t(sapply(SOP[[2]], function(x) x[1,]))
  SOP2.state2 = t(sapply(SOP[[2]], function(x) x[2,]))
  SOP2.state3 = t(sapply(SOP[[2]], function(x) x[3,]))
  SOP2.state4 = t(sapply(SOP[[2]], function(x) x[4,]))
  
  hat.SOP2.state1 = colMeans(SOP2.state1)
  band.SOP2.state1 = apply(SOP2.state1, 2, quantile, probs=c(0.025,0.975))
  hat.SOP2.state2 = colMeans(SOP2.state2)
  band.SOP2.state2 = apply(SOP2.state2, 2, quantile, probs=c(0.025,0.975))
  hat.SOP2.state3 = colMeans(SOP2.state3)
  band.SOP2.state3 = apply(SOP2.state3, 2, quantile, probs=c(0.025,0.975))
  hat.SOP2.state4 = colMeans(SOP2.state4)
  band.SOP2.state4 = apply(SOP2.state4, 2, quantile, probs=c(0.025,0.975))
  
  TP1.state1 = t(sapply(TP[[1]], function(x) x[1,]))
  TP1.state2 = t(sapply(TP[[1]], function(x) x[2,]))
  TP1.state3 = t(sapply(TP[[1]], function(x) x[3,]))
  
  hat.TP1.state1 = colMeans(TP1.state1)
  band.TP1.state1 = apply(TP1.state1, 2, quantile, probs=c(0.025,0.975))
  hat.TP1.state2 = colMeans(TP1.state2)
  band.TP1.state2 = apply(TP1.state2, 2, quantile, probs=c(0.025,0.975))
  hat.TP1.state3 = colMeans(TP1.state3)
  band.TP1.state3 = apply(TP1.state3, 2, quantile, probs=c(0.025,0.975))
  
  TP2.state1 = t(sapply(TP[[2]], function(x) x[1,]))
  TP2.state2 = t(sapply(TP[[2]], function(x) x[2,]))
  TP2.state3 = t(sapply(TP[[2]], function(x) x[3,]))
  
  hat.TP2.state1 = colMeans(TP2.state1)
  band.TP2.state1 = apply(TP2.state1, 2, quantile, probs=c(0.025,0.975))
  hat.TP2.state2 = colMeans(TP2.state2)
  band.TP2.state2 = apply(TP2.state2, 2, quantile, probs=c(0.025,0.975))
  hat.TP2.state3 = colMeans(TP2.state3)
  band.TP2.state3 = apply(TP2.state3, 2, quantile, probs=c(0.025,0.975))
  
  est.SOP = list(hat.SOP1.state1 = hat.SOP1.state1, hat.SOP1.state2 = hat.SOP1.state2,
                 hat.SOP1.state3 = hat.SOP1.state3, hat.SOP1.state4 = hat.SOP1.state4,
                 hat.SOP2.state1 = hat.SOP2.state1, hat.SOP2.state2 = hat.SOP2.state2, 
                 hat.SOP2.state3 = hat.SOP2.state3, hat.SOP2.state4 = hat.SOP2.state4,
                 band.SOP1.state1 = band.SOP1.state1, band.SOP1.state2 = band.SOP1.state2,
                 band.SOP1.state3 = band.SOP1.state3, band.SOP1.state4 = band.SOP1.state4,
                 band.SOP2.state1 = band.SOP2.state1, band.SOP2.state2 = band.SOP2.state2, 
                 band.SOP2.state3 = band.SOP2.state3, band.SOP2.state4 = band.SOP2.state4
  )
  
  est.TP = list(hat.TP1.state1 = hat.TP1.state1, hat.TP1.state2 = hat.TP1.state2,
                hat.TP1.state3 = hat.TP1.state3, hat.TP2.state1 = hat.TP2.state1,
                hat.TP2.state2 = hat.TP2.state2, hat.TP2.state3 = hat.TP2.state3,
                band.TP1.state1 = band.TP1.state1, band.TP1.state2 = band.TP1.state2,
                band.TP1.state3 = band.TP1.state3, band.TP2.state1 = band.TP2.state1,
                band.TP2.state2 = band.TP2.state2, band.TP2.state3 = band.TP2.state3
  )
  
  # plot friendly versions
  SOP = c(hat.SOP1.state1, hat.SOP1.state2, hat.SOP1.state3, hat.SOP1.state4,
         hat.SOP2.state1, hat.SOP2.state2, hat.SOP2.state3, hat.SOP2.state4)
  
  CB.SOP.1 = c(band.SOP1.state1[1,], band.SOP1.state2[1,],
               band.SOP1.state3[1,], band.SOP1.state4[1,],
               band.SOP2.state1[1,], band.SOP2.state2[1,],
               band.SOP2.state3[1,], band.SOP2.state4[1,])
  
  CB.SOP.2 = c(band.SOP1.state1[2,], band.SOP1.state2[2,],
               band.SOP1.state3[2,], band.SOP1.state4[2,],
               band.SOP2.state1[2,], band.SOP2.state2[2,],
               band.SOP2.state3[2,], band.SOP2.state4[2,])
  
  TP = c(hat.TP1.state1, hat.TP1.state2, hat.TP1.state3,
        hat.TP2.state1, hat.TP2.state2, hat.TP2.state3)
  
  CB.TP.1 = c(band.TP1.state1[1,], band.TP1.state2[1,], band.TP1.state3[1,],
              band.TP2.state1[1,], band.TP2.state2[1,], band.TP2.state3[1,])
  
  CB.TP.2 = c(band.TP1.state1[2,], band.TP1.state2[2,], band.TP1.state3[2,],
              band.TP2.state1[2,], band.TP2.state2[2,], band.TP2.state3[2,])
  
  return(est.prob = list(est.SOP = est.SOP, est.TP = est.TP, 
                         plot.SOP = SOP, CB.SOP.1 = CB.SOP.1, CB.SOP.2 = CB.SOP.2,
                         plot.TP = TP, CB.TP.1 = CB.TP.1, CB.TP.2 = CB.TP.2))
  
}
