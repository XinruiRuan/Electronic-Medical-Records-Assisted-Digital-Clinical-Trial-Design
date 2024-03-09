# R functions for data integration estimators                       
# input: data, parameters                                             
# output: the estimates                                 

# estimate working models (assume correct specification)
estimate_method<-function(df){
  fit.1 <- lm(Y~X1+X2,data=df) # linear expectation
  residuals_squared <- (residuals(fit.1))^2
  log_residuals_squared <- log(residuals_squared)
  fit.2<-lm(log_residuals_squared ~ X1 + X2, data = df) # log-linear variance
  
  return(list(fit.1=fit.1, fit.2=fit.2))
}

# estimate ATE
tau_hat <- function(obs.df, rct.df, emr = 1) {
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(cut(obs.df$X1,quantiles, labels = FALSE)==k)
  }
  N<-nrow(obs.df)
  n<-nrow(rct.df)
  df <- rbind(obs.df, rct.df) # collected dataset
  
  # linear regression
  fit.m1 <- estimate_method(rct.df[rct.df$D==1,])$fit.1
  # fit m0(X) with all control units
  fit.m0<-estimate_method(df[df$D==0,])$fit.1
  fit.V1<-estimate_method(rct.df[rct.df$D==0,])$fit.2
  fit.V0<-estimate_method(obs.df)$fit.2
  
  tau <- numeric(K)
  for (k in 1:K) {
    df.k <- df[df$X1 <= quantiles[k + 1] & df$X1 > quantiles[k], ]
    
    # estimates of some quantities
    m1.hat <- predict(fit.m1, df.k)
    m0.hat <- predict(fit.m0, df.k)
    V1.hat <- exp(predict(fit.V1, df.k))
    V0.hat <- exp(predict(fit.V0, df.k))
    r.hat <- V1.hat/V0.hat
    e.hat <- mean(df.k[df.k$Z == 1, ]$D)
    pi.hat <- mean(df.k$Z)
    
    # estimate of tau_k
    if (emr == 1) {
      # with EMR data
      tau[k] <-
        mean((1 - df.k$Z) / (1 - pi.hat) * (m1.hat - m0.hat) + df.k$D * df.k$Z /
               (pi.hat * e.hat) * (df.k$Y - m1.hat) - (df.k$Z * (1 - df.k$D) + (1 - df.k$Z) *
                                                         r.hat) / (pi.hat * (1 - e.hat) + (1 - pi.hat) * r.hat) * (df.k$Y - m0.hat)
        )
    }
    
    if (emr == 0) {
      # without EMR data
      tau[k] <-
        mean((1 - df.k$Z) / (1 - pi.hat) * (m1.hat - m0.hat) + df.k$D * df.k$Z /
               (pi.hat * e.hat) * (df.k$Y - m1.hat) - (df.k$Z * (1 - df.k$D)) / (pi.hat * (1 - e.hat)) * (df.k$Y - m0.hat)
        )
    }
  }
  tau.hat <- sum(tau * N.k / N) # weighted average
  return(c(tau, tau.hat))
}

# estimate standard error
std_hat<-function(obs.df,rct.df,emr=1){
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(cut(obs.df$X1,quantiles, labels = FALSE)==k)
  }
  N<-nrow(obs.df)
  n<-nrow(rct.df)
  df <- rbind(obs.df, rct.df) # collected dataset
  tau<-tau_hat(obs.df,rct.df)[1:K]
  tau.hat<-tau_hat(obs.df,rct.df)[K+1]
  sigma1.hat<-numeric(K)
  sigma0.hat1<-numeric(K)
  sigma0.hat0<-numeric(K)
  p.hat<-numeric(K)
  e.hat<-numeric(K)
  r.hat<-numeric(K)
  V<-numeric(K)
  for (k in 1:K) {
    df.k <- df[df$X1 <= quantiles[k + 1] & df$X1 > quantiles[k], ]
    sigma1.hat[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==1,]$Y))
    sigma0.hat1[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==0,]$Y))
    sigma0.hat0[k]<-sqrt(var(df.k[df.k$Z==0,]$Y))
    p.hat[k]<-mean(df.k$Z)
    e.hat[k]<-mean(df.k[df.k$Z == 1, ]$D)
    r.hat[k]<-sigma0.hat1[k]^2/sigma0.hat0[k]^2
    if (emr==1){ # with EMR data
      V[k]<-((tau[k] - tau.hat)^2 / (N.k[k] + p.hat[k] * n)) +
        ((sigma1.hat[k]^2) / (p.hat[k] * e.hat[k] * n)) +
        ((sigma0.hat1[k]^2) / (p.hat[k] * (1 - e.hat[k]) * n + N.k[k] * r.hat[k]))
    }
    if (emr==0){ # without EMR data
      V[k] <- ((tau[k] - tau.hat)^2 / (p.hat[k] * n)) +
        ((sigma1.hat[k]^2) / (p.hat[k] * e.hat[k] * n)) +
        ((sigma0.hat1[k]^2) / (p.hat[k] * (1 - e.hat[k]) * n))
    }
  }
  return(sqrt(sum(N.k^2 / N^2 * V)))
}





