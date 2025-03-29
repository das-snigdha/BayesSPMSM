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


est_SOP_TP_SI = function(M_mc, out, SI.new, tooth_num, t0, t_vec, t_vec_tp){
  
  SOP = vector(mode = "list", length = length(SI.new))
  TP = vector(mode = "list", length = length(SI.new))
  M_mcmc = length(out$alpha_dec)
  
  for(k in 1:length(SI.new)){
    print(paste("k :", k))
    
    SOP[[k]] = lapply(1:M_mcmc, function(i) 
      get_sop_SI(j = tooth.num[k], x =  SI.new[k], State = obs[, "State"],
                 n_MC = M_mc, t_vec = t_vec, alpha = out$alpha[i, ], 
                 mm = out$mm[i, ], ss = out$ss[i, ], Pi = out$pi[i, ],
                 Sigma_b = out$Sigma_b[i, ], xi = out$xi[i, ]) )
    
    TP[[k]] = lapply(1:M_mcmc, function(i) 
      get_tp_SI(j = tooth.num[k], x =  SI.new[k], State = obs[, "State"],
                n_MC = M_mc, t0 = t0, t_vec = t_vec_tp, alpha = out$alpha[i, ], 
                mm = out$mm[i, ], ss = out$ss[i, ], Pi = out$pi[i, ],
                Sigma_b = out$Sigma_b[i, ], xi = out$xi[i, ]) )
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
  
  SOP3.state1 = t(sapply(SOP[[3]], function(x) x[1,]))
  SOP3.state2 = t(sapply(SOP[[3]], function(x) x[2,]))
  SOP3.state3 = t(sapply(SOP[[3]], function(x) x[3,]))
  SOP3.state4 = t(sapply(SOP[[3]], function(x) x[4,]))
  
  hat.SOP3.state1 = colMeans(SOP3.state1)
  band.SOP3.state1 = apply(SOP3.state1, 2, quantile, probs=c(0.025,0.975))
  hat.SOP3.state2 = colMeans(SOP3.state2)
  band.SOP3.state2 = apply(SOP3.state2, 2, quantile, probs=c(0.025,0.975))
  hat.SOP3.state3 = colMeans(SOP3.state3)
  band.SOP3.state3 = apply(SOP3.state3, 2, quantile, probs=c(0.025,0.975))
  hat.SOP3.state4 = colMeans(SOP3.state4)
  band.SOP3.state4 = apply(SOP3.state4, 2, quantile, probs=c(0.025,0.975))
  
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
  
  TP3.state1 = t(sapply(TP[[3]], function(x) x[1,]))
  TP3.state2 = t(sapply(TP[[3]], function(x) x[2,]))
  TP3.state3 = t(sapply(TP[[3]], function(x) x[3,]))
  
  hat.TP3.state1 = colMeans(TP3.state1)
  band.TP3.state1 = apply(TP3.state1, 2, quantile, probs=c(0.025,0.975))
  hat.TP3.state2 = colMeans(TP3.state2)
  band.TP3.state2 = apply(TP3.state2, 2, quantile, probs=c(0.025,0.975))
  hat.TP3.state3 = colMeans(TP3.state3)
  band.TP3.state3 = apply(TP3.state3, 2, quantile, probs=c(0.025,0.975))
  
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