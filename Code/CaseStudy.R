# Case Study (HIV cash transfer)

library(ggplot2)
library(MASS)
library(mnormt)
library(locfit)
library(nloptr)

###### generate synthetic data #####

K<-2 # 0/1 gender
c1<-0.1

generate_data_0 <- function(n) {
  p0 <- c(0.38, 0.62)
  X <- factor(sample(c(1, 2), n, prob = p0, replace = TRUE))
  # discrete covariate representing WHO stages
  # X <- cbind(1, X)
  
  alpha.0 <- c(0.90, 0.91)  # alpha is the parameter for mean, beta is for std
  # we don't need beta, since the outcome is binary
  Y <- rbinom(n, size = 1, prob = alpha.0[X])
  Z <- rep(0, n)
  D <- rep(0, n)
  obs.df <- data.frame(Y = Y,
                       Z = Z,
                       D = D,
                       X = X)
  return(obs.df)
}

generate_data_1 <- function(n, p, e) {
  #Standardize p to prevent sum(p) from deviating from 1 due to slight computational errors
  p <- p / sum(p)
  
  alpha.1 <- c(0.96, 0.95)
  alpha.0 <- c(0.80, 0.92)
  
  X <- factor(sample(c(1,2), n, prob = p, replace = TRUE)) 
  # discrete covariate representing WHO stages
  # X <- cbind(1, X)
  K <- 2  # number of strata
  X <- c()
  D <- c()
  for (k in 1:K) {
    nk <- allocate_units(n,p)[k]
    Xk <- rep(k,nk)
    Dk <- rbinom(nk, 1, e[k]) # assign treatment indicator D
    X <- c(X, Xk)
    D <- c(D, Dk)
  }
  X<-factor(X)
  
  Y1 <- rbinom(n, size = 1, prob = alpha.1[X])
  Y0 <- rbinom(n, size = 1, prob = alpha.0[X])
  Y <- D * Y1 + (1 - D) * Y0  # observed Y
  Z <- rep(1, n)
  rct.df <- data.frame(
    Y = Y,
    Z = Z,
    D = D,
    X = X
  )
  return(rct.df)
}

# assess the transportability assumption (t test)
test<-function(obs.df,rct.df){
  p.value<-numeric(K)
  for (k in 1:K) {
    obs.df.k<-obs.df[(obs.df$X==k & obs.df$D==0),]
    rct.df.k<-rct.df[(rct.df$X==k & rct.df$D==0),]
    p.value[k]<-t.test(obs.df.k$Y, rct.df.k$Y, var.equal = FALSE)$p.value
  }
  return(p.value)
}

# p-values
obs_data <- generate_data_0(1000)
rct_data <- generate_data_1(1000, rep(1/K, K), rep(1/2, K))
p_values <- test(obs_data, rct_data)

# estimate working models
estimate_method<-function(df){
  fit.1 <- lm(Y~X-1,data=df)
  residuals_squared <- (residuals(fit.1))^2
  fit.2<-lm(residuals_squared ~ X-1, data = df) # log-linear variance
  
  return(list(fit.1=fit.1, fit.2=fit.2))
}

# estimate ATE
tau_hat <- function(obs.df, rct.df) {
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(obs.df$X==k)
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
    df.k <- df[df$X==k, ]
    
    # estimates of some quantities
    m1.hat <- predict(fit.m1, df.k)
    m0.hat <- predict(fit.m0, df.k)
    V1.hat <- predict(fit.V1, df.k)
    V0.hat <- predict(fit.V0, df.k)
    r.hat <- V1.hat/V0.hat
    e.hat <- mean(df.k[df.k$Z == 1, ]$D)
    pi.hat <- mean(df.k$Z)
    
    # estimate of tau_k
    tau[k] <- mean((1 - df.k$Z) / (1 - pi.hat) * (m1.hat - m0.hat) + df.k$D * df.k$Z /
                     (pi.hat * e.hat) * (df.k$Y - m1.hat) - (df.k$Z * (1 - df.k$D) + (1 - df.k$Z) *
                                                               r.hat) / (pi.hat * (1 - e.hat) + (1 - pi.hat) * r.hat) * (df.k$Y - m0.hat)
    )
    
    # estimate simply by experimental data

    # tau[k] <- mean((1 - df.k$Z) / (1 - pi.hat) * (m1.hat - m0.hat) + df.k$D * df.k$Z /
    #                  (pi.hat * e.hat) * (df.k$Y - m1.hat) - (df.k$Z * (1 - df.k$D)) / (pi.hat * (1 - e.hat)) * (df.k$Y - m0.hat)
    # )
  }
  tau.hat <- sum(tau * N.k / N) # weighted average
  return(c(tau, tau.hat))
}

# estimate standard error
std_hat<-function(obs.df,rct.df){
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(obs.df$X==k)
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
    df.k <- df[df$X==k, ]
    sigma1.hat[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==1,]$Y))
    sigma0.hat1[k]<-sqrt(var(df.k[df.k$Z==1 & df.k$D==0,]$Y))
    sigma0.hat0[k]<-sqrt(var(df.k[df.k$Z==0,]$Y))
    p.hat[k]<-mean(df.k$Z)
    e.hat[k]<-mean(df.k[df.k$Z == 1, ]$D)
    r.hat[k]<-sigma0.hat1[k]^2/sigma0.hat0[k]^2
    V[k]<-((tau[k] - tau.hat)^2 / (N.k[k] + p.hat[k] * n)) +
      ((sigma1.hat[k]^2) / (p.hat[k] * e.hat[k] * n)) +
      ((sigma0.hat1[k]^2) / (p.hat[k] * (1 - e.hat[k]) * n + N.k[k] * r.hat[k]))
  }
  return(sqrt(sum(N.k^2 / N^2 * V)))
}

# design strategy (Problem 3)
opt_design_2<-function(obs.df,rct.df,tau,tau.hat,sigma1.hat,sigma0.hat1,r.hat,n,Gamma0,Gamma1){
  df <- rbind(obs.df, rct.df) # collected dataset
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(obs.df$X==k)
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
        (Gamma1[k] + (N.k[k] * r.hat[k]) / (n * p[k] * (1 - e[k]) + N.k[k] * r.hat[k]) *
           Gamma0[k]) ^ 2
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
    N.k[k]<-sum(obs.df==k)
  }
  
  # estimate parameters
  tau<-tau_hat(obs.df,rct.df)[1:K]
  tau.hat<-tau_hat(obs.df,rct.df)[K+1]
  tau_list<-c(tau_list,tau.hat)
  sigma1.hat<-numeric(K)
  sigma0.hat1<-numeric(K)
  sigma0.hat0<-numeric(K)
  for (k in 1:K) {
    df.k <- df[df$X==k, ]
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
      df.k <- df[df$X==k, ]
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

### percentile bootstrap CI

boot.k <- function(df,Gamma0,Gamma1) {
  obs.df<-df[df$Z==0,]
  rct.df<-df[df$Z==1,]
  
  m1.hat <- mean(rct.df[rct.df$D==1,]$Y)
  m0.hat<-mean(df[df$D==0,]$Y)
  V1.hat<-var(rct.df[rct.df$D==0,]$Y)
  V0.hat<-var(obs.df$Y)
  
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

boot.ci<-function(obs.df,rct.df,Gamma0,Gamma1,B){
  N.k<-numeric(K)
  for (k in 1:K){
    N.k[k]<-sum(obs.df$X==k)
  }
  N<-nrow(obs.df)
  n<-nrow(rct.df)
  lower_bounds <- numeric(B)
  upper_bounds <- numeric(B)
  df <- rbind(obs.df, rct.df)
  for (b in 1:B) {
    lb<-numeric(K)
    ub<-numeric(K)
    for (k in 1:K) {
      df.k <- df[df$X==k,]
      sample_size <- nrow(df.k)
      sample_indices<-sample(1:sample_size, replace = TRUE)
      bootstrap_sample <- df.k[sample_indices,]
      tau.k <- boot.k(bootstrap_sample, Gamma0[k], Gamma1[k])[1]
      bias.k <- boot.k(bootstrap_sample, Gamma0[k], Gamma1[k])[2]
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

##### Percentile bootstrap confidence interval with respect to various choices of Gamma1 #####
gamma1_list <- c(0, 0.02, 0.04, 0.06)
lb_list <- numeric(length(gamma1_list))
ub_list <- numeric(length(gamma1_list))
df <-
  suppressWarnings(adaptive_design_2(500, c(250, 500, 750, 1000), c(0.1, 0), c(0, 0)))
obs.df <- df$obs.df
rct.df <- df$rct.df
for (i in 1:length(gamma1_list)) {
  ci <- boot.ci(obs.df, rct.df, c(0.1, 0), rep(gamma1_list[i], 2), 1000)
  lb_list[i] <- ci[1]
  ub_list[i] <- ci[2]
}

ci.df <- data.frame(gamma1 = factor(c(0, 0.02, 0.04, 0.06)),
                    lb <- lb_list,
                    # Mean values for each category
                    ub <- ub_list,
                    # Standard deviations for each category
                    type = factor(c("Type 1", "Type 1", "Type 1", 'Type 2')))

# plot
custom_colors <- c("Type 1" = "steelblue", "Type 2" = "gray")

p5 <-
  ggplot(ci.df, aes(
    x = gamma1,
    y = (lb + ub) / 2,
    color = type
  )) +
  geom_point(position = position_dodge(0.9), size = 2) +
  geom_errorbar(
    aes(ymin = lb, ymax = ub),
    position = position_dodge(0.9),
    size = 1,
    width = 0.5
  ) +
  labs(x = "Gamma1", y = "Treatment Effect Estimates") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "black") +
  scale_color_manual(values = custom_colors) +
  ggtitle("Percentile Bootstrap Confidence Interval") +
  guides(color = FALSE) +
  theme_minimal()

ggsave("case1.pdf",
       p5,
       width = 4,
       height = 4,
       units = "in")
