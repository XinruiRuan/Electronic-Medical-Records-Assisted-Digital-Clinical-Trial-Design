# R functions for digital clinical trial design                       
# input: parameters (adaptive experiment setups)                                           
# output: the design strategies and design-based estimators     

##### Design strategy in Problem 2 (perfect transportability) #####
opt_design<-function(obs.df,rct.df,tau,tau.hat,sigma1.hat,sigma0.hat1,r.hat,n,emr=1){
  df <- rbind(obs.df, rct.df) # collected dataset
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(cut(obs.df$X1,quantiles, labels = FALSE)==k)
  }
  N<-nrow(obs.df)
  # objective function
  eval_f0 <- function(x){
    p <- x[1:K]  # subgroup proportion
    e <- x[(K+1):(2*K)]  # propensity score
    
    V <- numeric(K)
    # EMR data assisted
    for (k in 1:K) {
      if (emr==1){
        V[k] <- ((tau[k] - tau.hat)^2 / (N.k[k] + p[k] * n)) +
          ((sigma1.hat[k]^2) / (p[k] * e[k] * n)) +
          ((sigma0.hat1[k]^2) / (p[k] * (1 - e[k]) * n + N.k[k] * r.hat[k]))
      }
      if (emr==0){
        V[k] <- ((tau[k] - tau.hat)^2 / (p[k] * n)) +
          ((sigma1.hat[k]^2) / (p[k] * e[k] * n)) +
          ((sigma0.hat1[k]^2) / (p[k] * (1 - e[k]) * n))
      }
    }
    
    return(sum(N.k^2 / N^2 * V))
  }
  
  # Define the constraint function
  eval_g0 <- function(x) {
    p <- x[1:K]
    e <- x[(K+1):(2*K)]
    return(c(sum(p)-1, -p,c1-e, -1 + e + c1))
  }
  
  # Solve using NLOPT_LN_COBYLA without gradient information
  res1 <- nloptr(x0=c(rep(1/K,K),rep(1/2,K)),
                 eval_f=eval_f0,
                 lb = rep(0,2*K),
                 ub = rep(1,2*K),
                 eval_g_ineq = eval_g0,
                 opts = list("algorithm" = "NLOPT_LN_COBYLA",
                             "xtol_rel"=1.0e-6,
                             "print_level" = 0,
                             "check_derivatives" = TRUE,
                             "warn_every_iter"=FALSE,
                             "check_derivatives_print" = "all"))
  
  # optimal p* and e*
  return(res1$solution)
}

adaptive_design <- function(N,n,emr=1){   
  # N denotes the number of units in external data
  # n is a vector denoting the number of units in each stage of the experiments
  
  T<-length(n)
  p<-matrix(NA, nrow = T, ncol = K)
  e<-matrix(NA, nrow = T, ncol = K)
  tau_list<-c()
  std_list<-c()
  
  # Stage 1
  p[1,]<-rep(1/K,K)
  e[1,]<-rep(1/2,K)
  obs.df<-generate_data_0(N)
  rct.df<-generate_data_1(n[1],p[1,],e[1,])
  df<-rbind(obs.df,rct.df)
  
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(cut(obs.df$X1,quantiles, labels = FALSE)==k)
  }
  
  # initialize parameters
  tau<-tau_hat(obs.df,rct.df,emr)[1:K]
  tau.hat<-tau_hat(obs.df,rct.df,emr)[K+1]
  tau_list<-c(tau_list,tau.hat)
  sigma1.hat<-numeric(K)
  sigma0.hat1<-numeric(K)
  sigma0.hat0<-numeric(K)
  for (k in 1:K) {
    df.k <- df[df$X1 <= quantiles[k + 1] & df$X1 > quantiles[k], ]
    sigma1.hat[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==1,]$Y))
    sigma0.hat1[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==0,]$Y))
    sigma0.hat0[k]<-sqrt(var(df.k[df.k$Z==0,]$Y))
  }
  r.hat<-sigma0.hat1^2/sigma0.hat0^2
  std_list<-c(std_list,std_hat(obs.df,rct.df,emr))
  
  # multi-stage design
  for(t in 1:(T-1)){ # start of sequential assignments
    
    p_opt <- opt_design(obs.df,rct.df,tau,tau.hat,sigma1.hat,sigma0.hat1,r.hat,sum(n[1:T]),emr)[1:K]
    e_opt <- opt_design(obs.df,rct.df,tau,tau.hat,sigma1.hat,sigma0.hat1,r.hat,sum(n[1:T]),emr)[(K+1):(2*K)]
    
    # rescale
    pt<-(p_opt*sum(n[1:(t+1)])-p[t,]*sum(n[1:t]))/n[t+1]
    et<-(e_opt*p_opt*sum(n[1:(t+1)])-e[t,]*p[t,]*sum(n[1:t]))/(n[t+1]*pt)
    pt<-ifelse(pt < 0.01, 0.01, pt)
    pt<-pt/sum(pt)
    et<-ifelse(et < c1, c1, et)
    et<-ifelse(et > 1-c1, 1-c1, et)
    
    
    # sample
    rct.t<-generate_data_1(n[t+1],pt,et)
    rct.df<-rbind(rct.df,rct.t)
    df<-rbind(obs.df,rct.df)
    
    # update parameters
    tau<-tau_hat(obs.df,rct.df,emr)[1:K]
    tau.hat<-tau_hat(obs.df,rct.df,emr)[K+1]
    tau_list<-c(tau_list,tau.hat)
    std_list<-c(std_list,std_hat(obs.df,rct.df,emr))
    p[t+1,]<-p_opt
    e[t+1,]<-e_opt
    sigma1.hat<-numeric(K)
    sigma0.hat1<-numeric(K)
    sigma0.hat0<-numeric(K)
    for (k in 1:K) {
      df.k <- df[df$X1 <= quantiles[k + 1] & df$X1 > quantiles[k], ]
      sigma1.hat[k]<-sqrt(var(df[df$Z==1 & df$D==1,]$Y))
      sigma0.hat1[k]<-sqrt(var(df[df$Z==1 & df$D==0,]$Y))
      sigma0.hat0[k]<-sqrt(var(df[df$Z==0,]$Y))
    }
    r.hat<-sigma0.hat1^2/sigma0.hat0^2
  }
  
  result_list <-
    list(
      obs.df = obs.df,
      rct.df = rct.df,
      p = p,
      e = e,
      tau_list = tau_list
    )
  return(result_list)
}


##### Design strategy in Problem 3 (imperfect transportability) #####

opt_design_2<-function(obs.df,rct.df,tau,tau.hat,sigma1.hat,sigma0.hat1,r.hat,n,Gamma0,Gamma1){
  df <- rbind(obs.df, rct.df) # collected dataset
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(cut(obs.df$X1,quantiles, labels = FALSE)==k)
  }
  N<-nrow(obs.df)
  # objective function
  eval_f0 <- function(x){
    p <- x[1:K]  # subgroup proportion
    e <- x[(K+1):(2*K)]  # propensity score
    
    V <- numeric(K)
    # external data assisted
    for (k in 1:K) {
      V[k] <- ((tau[k] - tau.hat) ^ 2 / (N.k[k] + p[k] * n)) +
        ((sigma1.hat[k] ^ 2) / (p[k] * e[k] * n)) +
        ((sigma0.hat1[k] ^ 2) / (p[k] * (1 - e[k]) * n + N.k[k] * r.hat[k])) +
        (Gamma1 + (N.k[k] * r.hat[k]) / (n * p[k] * (1 - e[k]) + N.k[k] * r.hat[k]) *
           Gamma0) ^ 2
    }
    
    return(sum(N.k^2 / N^2 * V))
  }
  
  # Define the constraint function
  eval_g0 <- function(x) {
    p <- x[1:K]
    e <- x[(K+1):(2*K)]
    return(c(sum(p)-1, -p,c1-e, -1 + e + c1))
  }
  
  # Solve using NLOPT_LN_COBYLA without gradient information
  res1 <- nloptr(x0=c(rep(1/K,K),rep(1/2,K)),
                 eval_f=eval_f0,
                 lb = rep(0,2*K),
                 ub = rep(1,2*K),
                 eval_g_ineq = eval_g0,
                 opts = list("algorithm" = "NLOPT_LN_COBYLA",
                             "xtol_rel"=1.0e-6,
                             "print_level" = 0,
                             "check_derivatives" = TRUE,
                             "warn_every_iter"=FALSE,
                             "check_derivatives_print" = "all"))
  
  # optimal p* and e*
  return(res1$solution)
}


adaptive_design_2 <- function(N,n,Gamma0,Gamma1){   
  # N denotes the number of units in external data
  # n is a vector denoting the number of units in each stage of the experiments
  
  T<-length(n)
  p<-matrix(NA, nrow = T, ncol = K)
  e<-matrix(NA, nrow = T, ncol = K)
  tau_list<-c()
  std_list<-c()
  
  # Stage 1
  p[1,]<-rep(1/K,K)
  e[1,]<-rep(1/2,K)
  obs.df<-generate_data_0(N)
  rct.df<-generate_data_1(n[1],p[1,],e[1,])
  df<-rbind(obs.df,rct.df)
  
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(cut(obs.df$X1,quantiles, labels = FALSE)==k)
  }
  
  # estimate parameters
  tau<-tau_hat(obs.df,rct.df)[1:K]
  tau.hat<-tau_hat(obs.df,rct.df)[K+1]
  tau_list<-c(tau_list,tau.hat)
  sigma1.hat<-numeric(K)
  sigma0.hat1<-numeric(K)
  sigma0.hat0<-numeric(K)
  for (k in 1:K) {
    df.k <- df[df$X1 <= quantiles[k + 1] & df$X1 > quantiles[k], ]
    sigma1.hat[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==1,]$Y))
    sigma0.hat1[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==0,]$Y))
    sigma0.hat0[k]<-sqrt(var(df.k[df.k$Z==0,]$Y)) # sigma and rhat should be in the update stage, too
  }
  r.hat<-sigma0.hat1^2/sigma0.hat0^2
  std_list<-c(std_list,std_hat(obs.df,rct.df))
  
  # multi-stage design
  for(t in 1:(T-1)){ # start of sequential assignments
    
    p_opt <- opt_design_2(obs.df,rct.df,tau,tau.hat,sigma1.hat,sigma0.hat1,r.hat,sum(n[1:T]),Gamma0,Gamma1)[1:K]
    e_opt <- opt_design_2(obs.df,rct.df,tau,tau.hat,sigma1.hat,sigma0.hat1,r.hat,sum(n[1:T]),Gamma0,Gamma1)[(K+1):(2*K)]
    
    # rescale
    pt<-(p_opt*sum(n[1:(t+1)])-p[t,]*sum(n[1:t]))/n[t+1]
    et<-(e_opt*p_opt*sum(n[1:(t+1)])-e[t,]*p[t,]*sum(n[1:t]))/(n[t+1]*pt)
    pt<-ifelse(pt < 0.01, 0.01, pt)
    pt<-pt/sum(pt)
    et<-ifelse(et < c1, c1, et)
    et<-ifelse(et > 1-c1, 1-c1, et)
    
    
    # sample
    rct.t<-generate_data_1(n[t+1],pt,et)
    rct.df<-rbind(rct.df,rct.t)
    df<-rbind(obs.df,rct.df)
    
    # update
    tau<-tau_hat(obs.df,rct.df)[1:K]
    tau.hat<-tau_hat(obs.df,rct.df)[K+1]
    tau_list<-c(tau_list,tau.hat)
    std_list<-c(std_list,std_hat(obs.df,rct.df))
    p[t+1,]<-p_opt
    e[t+1,]<-e_opt
    sigma1.hat<-numeric(K)
    sigma0.hat1<-numeric(K)
    sigma0.hat0<-numeric(K)
    for (k in 1:K) {
      df.k <- df[df$X1 <= quantiles[k + 1] & df$X1 > quantiles[k], ]
      sigma1.hat[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==1,]$Y))
      sigma0.hat1[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==0,]$Y))
      sigma0.hat0[k]<-sqrt(var(df.k[df.k$Z==0,]$Y)) # sigma and rhat should be in the update stage, too
    }
    r.hat<-sigma0.hat1^2/sigma0.hat0^2
  }
  result_list <-
    list(
      obs.df = obs.df,
      rct.df = rct.df,
      p = p,
      e = e,
      tau_list = tau_list
    )
  return(result_list)
}



