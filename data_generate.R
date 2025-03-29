# function to generate data from the correctly specified model
data_model = function(n=100, m=14, p=3, snr = 5, mu_eps = c(-0.5, 0, 0.5),
                         true.alpha, true.beta, C.param){
  
  N=n*m; p1 = floor((p-1)/2); p2 = (p-1)-p1
  
  tmp=diag(rep(2,m/2)); tmp[1,1]=tmp[m/2,m/2]=1
  for(j in 1:(m/2-1))
    tmp[j,j+1] = tmp[j+1,j] = -0.9
  true.Sigma.b = solve(tmp)/100
  
  tmp = cbind(1:n, matrix(runif(p1*n, min = -5, max = 5), ncol = p1), 
              matrix(rbinom(p2*n, 1, 0.5), ncol=p2))
  obs.full = cbind(tmp, 1)
  for(i in 2:m)
    obs.full = rbind(obs.full, cbind(tmp, i))
  xnames = paste0("x", 1:(p-1))
  colnames(obs.full) = c("id", xnames, "tooth")
  
  obs.full = obs.full[order(obs.full[,"id"], obs.full[,"tooth"]), ]
  
  obs.full = cbind(obs.full, upperjaw=0)
  obs.full[obs.full[,"tooth"] > m/2, "upperjaw"] = 1
  
  colnames(obs.full) = c("id", xnames, "tooth", paste0("x",p))
  
  xnames = c(xnames, paste0("x",p))
  obs.full = obs.full[,c("id", "tooth", xnames)]
  
  x.mtx = as.matrix(obs.full[, xnames])
  
  x.scaled = scale(x.mtx, center = TRUE, scale = TRUE)
  tmp = max(sqrt(rowSums(x.scaled*x.scaled)))
  x.scaled = x.scaled / tmp
  obs.full[, xnames] = x.scaled
  x.scaled.wts = tmp
  
  # b
  tmp = mvrnorm(n, rep(0,m/2), true.Sigma.b)
  tmp = cbind(tmp, tmp[,(m/2):1])
  obs.full = cbind(obs.full, b=as.vector(t(tmp)))
  
  n.mu = length(mu_eps)
  ind = sample(1:n.mu, N, T, prob = rep(1/n.mu, each = n.mu))
  eps = rep(0, N)
  
  mu = mu_eps
  for(i in 1:N)
    eps[i] = rnorm(1, mean = mu[ind[i]], sd = 0.1)
  
  obs.full = cbind(obs.full, eps=eps)
  
  # SI
  obs.full = cbind(obs.full, SI=(obs.full[, xnames] %*% true.beta)[,1])
  
  # determine constant for g
  var.error = var(obs.full[, "b"]) + var(obs.full[, "eps"])
  
  var.g = var(true.g(obs.full[,"SI"], a=1))
  const.g = sqrt(snr*var.error / var.g)
  
  # g(SI)
  gSI = true.g(obs.full[,"SI"], a = const.g)
  
  # logT
  obs.full = cbind(obs.full, logT=obs.full[,"b"] + gSI + obs.full[,"eps"], gSI = gSI)
  obs.full = cbind(obs.full, C=rgamma(N, C.param[1], C.param[2]))
  
  
  # R, T_k, S
  tmp = rdirichlet(N, true.alpha)
  colnames(tmp) = c("R1", "R2", "R3")
  obs.full = cbind(obs.full, tmp)
  
  obs.full = cbind(obs.full, T1=0, T2=0, T3=exp(obs.full[,"logT"]), State=-1)
  obs.full[,"T1"] = obs.full[,"T3"] * obs.full[,"R1"]
  obs.full[,"T2"] = obs.full[,"T3"] * (obs.full[,"R1"]+obs.full[,"R2"])
  
  ind0 = (obs.full[,"C"] < obs.full[,"T1"])
  ind1 = (obs.full[,"C"] >= obs.full[,"T1"]) & (obs.full[,"C"] < obs.full[,"T2"])
  ind2 = (obs.full[,"C"] >= obs.full[,"T2"]) & (obs.full[,"C"] < obs.full[,"T3"])
  ind3 = (obs.full[,"C"] >= obs.full[,"T3"])
  obs.full[ind0, "State"] = 0
  obs.full[ind1, "State"] = 1
  obs.full[ind2, "State"] = 2
  obs.full[ind3, "State"] = 3
  
  return(list(obs.full = obs.full, const.g = const.g,
              x.scaled = x.scaled, x.scaled.wts = x.scaled.wts))
  
}

# function to generate data from the misspecified model
data_misspecification = function(n=100, m=10, p=4, snr = 5, nu = 3, error.percent = 0.9, 
                                 true.alpha, true.beta, C.param){
  
  N=n*m; p1 = floor((p-1)/2); p2 = (p-1)-p1
  
  tmp=diag(rep(2,m/2)); tmp[1,1]=tmp[m/2,m/2]=1
  for(j in 1:(m/2-1))
    tmp[j,j+1] = tmp[j+1,j] = -0.9
  true.Sigma.b = solve(tmp)/100
  
  tmp = cbind(1:n, matrix(runif(p1*n, min = -5, max = 5), ncol = p1), 
              matrix(rbinom(p2*n, 1, 0.5), ncol=p2))
  obs.full = cbind(tmp, 1)
  for(i in 2:m)
    obs.full = rbind(obs.full, cbind(tmp, i))
  xnames = paste0("x", 1:(p-1))
  colnames(obs.full) = c("id", xnames, "tooth")
  
  obs.full = obs.full[order(obs.full[,"id"], obs.full[,"tooth"]), ]
  
  obs.full = cbind(obs.full, upperjaw=0)
  obs.full[obs.full[,"tooth"] > m/2, "upperjaw"] = 1
  
  colnames(obs.full) = c("id", xnames, "tooth", paste0("x",p))
  
  xnames = c(xnames, paste0("x",p))
  obs.full = obs.full[,c("id", "tooth", xnames)]
  
  x.mtx = as.matrix(obs.full[, xnames])
  
  x.scaled = scale(x.mtx, center = TRUE, scale = TRUE)
  tmp = max(sqrt(rowSums(x.scaled*x.scaled)))
  x.scaled = x.scaled / tmp
  obs.full[, xnames] = x.scaled
  x.scaled.wts = tmp
  
  ## b
  tmp = mvrnorm(n, rep(0,m/2), true.Sigma.b)
  tmp = cbind(tmp, tmp[,(m/2):1])
  obs.full = cbind(obs.full, b=as.vector(t(tmp)))
  
  # eps
  ind = sample(1:2, N, T, prob = c(error.percent, 1-error.percent))
  eps = rep(0, N)
  
  for(i in 1:N){
    if(ind[i]==1)
      eps[i] = rnorm(1, mean = 0, sd = 0.1)
    else eps[i] = rt(1, nu)
  }
  
  obs.full = cbind(obs.full, eps=eps)
  
  # SI
  obs.full = cbind(obs.full, SI=(obs.full[, xnames] %*% true.beta)[,1])
  
  # determine constant for g
  var.error = var(obs.full[, "b"]) + var(obs.full[, "eps"])
  
  var.g = var(true.g(obs.full[,"SI"], a=1))
  const.g = sqrt(snr*var.error / var.g)
  
  # g(SI)
  gSI = true.g(obs.full[,"SI"], a = const.g)
  
  # logT
  obs.full = cbind(obs.full, logT=obs.full[,"b"] + gSI + obs.full[,"eps"], gSI = gSI)
  obs.full = cbind(obs.full, C=rgamma(N, C.param[1], C.param[2]))
  
  
  # R, T_k, S
  tmp = rdirichlet(N, true.alpha)
  colnames(tmp) = c("R1", "R2", "R3")
  obs.full = cbind(obs.full, tmp)
  
  obs.full = cbind(obs.full, T1=0, T2=0, T3=exp(obs.full[,"logT"]), State=-1)
  obs.full[,"T1"] = obs.full[,"T3"] * obs.full[,"R1"]
  obs.full[,"T2"] = obs.full[,"T3"] * (obs.full[,"R1"]+obs.full[,"R2"])
  
  ind0 = (obs.full[,"C"] < obs.full[,"T1"])
  ind1 = (obs.full[,"C"] >= obs.full[,"T1"]) & (obs.full[,"C"] < obs.full[,"T2"])
  ind2 = (obs.full[,"C"] >= obs.full[,"T2"]) & (obs.full[,"C"] < obs.full[,"T3"])
  ind3 = (obs.full[,"C"] >= obs.full[,"T3"])
  obs.full[ind0, "State"] = 0
  obs.full[ind1, "State"] = 1
  obs.full[ind2, "State"] = 2
  obs.full[ind3, "State"] = 3
  
  return(list(obs.full = obs.full, const.g = const.g,
              x.scaled = x.scaled, x.scaled.wts = x.scaled.wts))
  
}
