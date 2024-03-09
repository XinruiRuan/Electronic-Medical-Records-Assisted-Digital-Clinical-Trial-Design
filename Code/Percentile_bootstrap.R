# Percentile bootstrap confidence intervals
# input: data and parameters                                           
# output: confidence interval 

##### estimation based on bootstrapped samples #####
boot.k <- function(df,Gamma0,Gamma1) {
  obs.df<-df[df$Z==0,]
  rct.df<-df[df$Z==1,]
  
  fit.m1 <- estimate_method(rct.df[rct.df$D==1,])$fit.1
  fit.m0<-estimate_method(df)$fit.1
  fit.V1<-estimate_method(rct.df[rct.df$D==0,])$fit.2
  fit.V0<-estimate_method(obs.df)$fit.2
  
  m1.hat <- predict(fit.m1, df)
  m0.hat <- predict(fit.m0, df)
  V1.hat <- exp(predict(fit.V1, df))
  V0.hat <- exp(predict(fit.V0, df))
  r.hat <- V1.hat/V0.hat
  e.hat <- mean(df[df$Z == 1, ]$D)
  pi.hat <- mean(df$Z)
  
  # estimate of tau_k
  tau.k<-mean((1 - df$Z) / (1 - pi.hat) * (m1.hat - m0.hat) + df$D * df$Z /
                (pi.hat * e.hat) * (df$Y - m1.hat) - (df$Z * (1 - df$D) + (1 - df$Z) *
                                                        r.hat) / (pi.hat * (1 - e.hat) + (1 - pi.hat) * r.hat) * (df$Y - m0.hat)
  )
  bias.k<-Gamma1 + mean(((1-df$Z) * r.hat) / (pi.hat * (1 - e.hat) + (1-pi.hat) * r.hat)) *
    Gamma0
  return(c(tau.k,bias.k))
}

##### percentile bootstrap ci #####

boot.ci<-function(obs.df,rct.df,Gamma0,Gamma1,B){
  lower_bounds <- numeric(B)
  upper_bounds <- numeric(B)
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(cut(obs.df$X1,quantiles, labels = FALSE)==k)
  }
  N<-nrow(obs.df)
  df <- rbind(obs.df, rct.df)
  for (b in 1:B) {
    lb<-numeric(K)
    ub<-numeric(K)
    for (k in 1:K) {
      df.k <- df[df$X1 <= quantiles[k + 1] & df$X1 > quantiles[k],]
      sample_size <- nrow(df.k)
      sample_indices<-sample(1:sample_size, replace = TRUE)
      bootstrap_sample <- df.k[sample_indices,]
      tau.k <- boot.k(bootstrap_sample, Gamma0, Gamma1)[1]
      bias.k <- boot.k(bootstrap_sample, Gamma0, Gamma1)[2]
      lb[k] <- tau.k - bias.k
      ub[k] <- tau.k + bias.k
    }
    
    lower_bound <- sum(N.k / N * lb)
    upper_bound <- sum(N.k / N * ub)
    
    lower_bounds[b] <- lower_bound
    upper_bounds[b] <- upper_bound
  }
  
  alpha <- 0.05
  lower_percentile <- quantile(lower_bounds, alpha / 2)
  upper_percentile <- quantile(upper_bounds, 1 - alpha / 2)
  
  return(c(lower_percentile, upper_percentile))
}
