# R code for simulation studies

library(MASS)
library(mnormt)
library(locfit)
library(nloptr)

source('Estimators.R')
source('Adaptive_allocation.R')
source('Percentile_bootstrap.R')

##### Under transportability assumption #####

##### DGP #####

# ensure that the sum adds up to n while considering potential computation errors
allocate_units <- function(n, p) {
  subgroup_units <- floor(p * n)
  remaining_units <- n - sum(subgroup_units)
  rounding_errors <- p * n - subgroup_units
  sorted_errors <- sort(rounding_errors, decreasing = TRUE)
  for (i in 1:remaining_units) {
    if (remaining_units > 0) {
      subgroup_units[which(rounding_errors == sorted_errors[i])[1]] <- subgroup_units[which(rounding_errors == sorted_errors[i])[1]] + 1
      remaining_units <- remaining_units - 1
    } else {
      break
    }
  }
  return(subgroup_units)
}


# generate EMR data from general population
generate_data_0 <- function(n) {
  X <- mvrnorm(n, c(0, 0), matrix(c(1, 0.1, 0.1, 1), nrow = 2)) # 2-dimensional covariates
  X <- cbind(1, X)
  
  alpha.0 <- c(1, 1.5, 1)  # alpha is the parameter for mean, beta is for std 
  beta.0 <- c(0.5, 0.5, 0.5)
  Y <- X %*% alpha.0 + rnorm(n, sd = sqrt(exp(X %*% beta.0))) # only control data in the external dataset
  Z <- rep(0, n)
  D <- rep(0, n)
  obs.df <- data.frame(
    Y = Y,
    Z = Z,
    D = D,
    X1 = X[, 2],
    X2 = X[, 3]
  )
  return(obs.df)
}

# generate RCT data from trial population
generate_data_1 <- function(n, p, e) {
  # standardize p to prevent sum(p) from deviating from 1 due to slight computational errors
  p <- p / sum(p)
  
  alpha.1 <- c(3, 1, 0.5)
  beta.1 <- c(0.5, 0.5,0.5)
  alpha.0 <- c(1, 1.5, 1)
  beta.0 <- c(2, 0.5,0.5)

  X <- c()
  D <- c()
  for (k in 1:K) {
    nk <- allocate_units(n,p)[k]
    lb <- c(quantiles[k], -Inf)
    ub <- c(quantiles[k + 1], Inf)
    Xk <-
      rmtruncnorm(nk, c(0.5, 0.5), matrix(c(1, 0.1, 0.1, 1), nrow = 2), lb, ub) # sample covariates from k'th stratum
    Dk <- rbinom(nk, 1, e[k]) # assign treatment indicator D
    X <- rbind(X, Xk)
    D <- c(D, Dk)
  }
  
  X <- cbind(1, X)
  Y1 <- X %*% alpha.1 + rnorm(n, sd = sqrt(exp(X %*% beta.1)))
  Y0 <- X %*% alpha.0 + rnorm(n, sd = sqrt(exp(X %*% beta.0)))
  Y <- D * Y1 + (1 - D) * Y0  # observed Y
  Z <- rep(1, n)
  rct.df <- data.frame(
    Y = Y,
    Z = Z,
    D = D,
    X1 = X[, 2],
    X2 = X[, 3]
  )
  return(rct.df)
}

##### parameters #####
K<-4 # number of strata
T<-4
c1<-0.1
# X1 quantiles for stratification
quantiles <- qnorm(seq(0, 1, 0.25), mean = 0.5, sd = 1)


# Monte-Carlo experiment for the efficiency comparison
M <- 100 
# number of randomizations (You can alter this number as needed)

tau_mc_1 <- matrix(data = NA, nrow = M, ncol = T)
tau_mc_0 <- matrix(data = NA, nrow = M, ncol = T)
for (m in 1:M) {
  tau_mc_1[m, ] <-
    suppressWarnings(adaptive_design(500, c(250, 500, 750, 1000), 1)$tau_list)
  tau_mc_0[m, ] <-
    suppressWarnings(adaptive_design(500, c(250, 500, 750, 1000), 0)$tau_list)
}

# save the results
write.csv(tau_mc_1, file = 'tau_mc_1.csv', row.names = FALSE)
write.csv(tau_mc_0, file = 'tau_mc_0.csv', row.names = FALSE)


# Monte-Carlo experiment for design strategies
n_list <- seq(100, 2500, by = 200)
p_mc<-matrix(data=NA,nrow=M,ncol=length(n_list))
e_mc<-matrix(data=NA,nrow=M,ncol=length(n_list))
for (i in 1:length(n_list)) {
  for (m in 1:M) {
    design<-suppressWarnings(adaptive_design(n_list[i],c(n_list[i],n_list[i]*2)))
    p_mc[m,i]<-design$p[2,1]
    e_mc[m,i]<-design$e[2,1]
  }
}

# save the results
bias <- c(apply(p_mc,2,mean),apply(e_mc,2,mean))
lower_ci <- bias-1.96*sqrt(c(apply(p_mc,2,var),apply(e_mc,2,var)))
upper_ci <- bias+1.96*sqrt(c(apply(p_mc,2,var),apply(e_mc,2,var)))

df_max <- data.frame(
  m = rep(n_list,2),
  bias = bias,
  Design = c(rep('p1',length(n_list)),rep('e1',length(n_list))),
  lower_ci = lower_ci,
  upper_ci = upper_ci
)

write.csv(df_max,file='df_max.csv',row.names = FALSE)




###### violation of transportability assumption #######

# bias parameters
u0 <- 0.1
delta <- 0

# DGP

generate_data_1 <- function(n, p, e) {
  #Standardize p to prevent sum(p) from deviating from 1 due to slight computational errors
  p <- p / sum(p)
  
  alpha.1 <- c(3 + u0 + delta, 1, 0.5)
  beta.1 <- c(0.5, 0.5, 0.5)
  alpha.0 <- c(1 + u0, 1.5, 1)
  beta.0 <- c(2, 0.5, 0.5)
  
  quantiles <-
    qnorm(seq(0, 1, 0.25), mean = 0.5, sd = 1) # quantiles for stratification
  K <- 4 # number of strata
  X <- c()
  D <- c()
  for (k in 1:K) {
    nk <- allocate_units(n,p)[k]
    lb <- c(quantiles[k], -Inf)
    ub <- c(quantiles[k + 1], Inf)
    Xk <-
      rmtruncnorm(nk, c(0.5, 0.5), matrix(c(1, 0.1, 0.1, 1), nrow = 2), lb, ub) # sample covariates from k'th stratum
    Dk <- rbinom(nk, 1, e[k]) # assign treatment indicator D
    X <- rbind(X, Xk)
    D <- c(D, Dk)
  }
  
  X <- cbind(1, X)
  Y1 <- X %*% alpha.1 + rnorm(n, sd = sqrt(exp(X %*% beta.1)))
  Y0 <- X %*% alpha.0 + rnorm(n, sd = sqrt(exp(X %*% beta.1)))
  Y <- D * Y1 + (1 - D) * Y0  # observed Y
  Z <- rep(1, n)
  rct.df <- data.frame(
    Y = Y,
    Z = Z,
    D = D,
    X1 = X[, 2],
    X2 = X[, 3]
  )
  return(rct.df)
}

# Design strategies comparison between Problems 2 and 3
design2_star <-
  suppressWarnings(adaptive_design_2(500, c(250, 500, 750, 1000), 0.5, 0.5))
design2_star <- c(design2_star$p[T, ], design2_star$e[T, ])
write.csv(design2_star, "design2_star.csv", row.names = FALSE)
design1_star <-
  suppressWarnings(adaptive_design(500, c(250, 500, 750, 1000)))
design1_star <- c(design1_star$p[T, ], design1_star$e[T, ])

write.csv(design1_star, "design1_star.csv", row.names = FALSE)


# Monte-Carlo experiment for testing the validity of percentile bootstrap
m <- 10 
# You can alter this number of randomizations, but it's time-costly; in the ultimate simulation, m = 100
ci_matrix <- matrix(data = NA, nrow = m, ncol = 2)
indicator <- numeric(m)
for (j in 1:m) {
  result_list <-
    suppressWarnings(adaptive_design_2(500, c(250, 500, 750, 1000), 0.5, 0.5))
  obs.df <- result_list$obs.df
  rct.df <- result_list$rct.df
  ci <-
    boot.ci(obs.df, rct.df, Gamma0 = 0, Gamma1 = 0.5, 100) 
  # here you can alter Gamma0 and Gamma1 (if Gamma0>=0.5, Gamma1>=0.5, then coverage rate > 0.95)
  ci_matrix[j, ] <- ci
  if (2 >= ci[1] && 2 <= ci[2]) {
    indicator[j] <- 1
  }
}

# empirical coverage rate
empirical_coverage_rate <- sum(indicator) / m
cat("Empirical Coverage Rate:", empirical_coverage_rate, "\n")


# Additional simulations (set delta=0 and revise u0)
M <- 1000 # time costly
T<-2
tau_mc_b <- matrix(data = NA, nrow = M, ncol = T) # for design B
tau_mc_c <- matrix(data = NA, nrow = M, ncol = T) # for design C
for (m in 1:M) {
  tau_mc_c[m, ] <-
    suppressWarnings(adaptive_design_2(500, c(250, 250),0.1,0)$tau_list)
  tau_mc_b[m, ] <-
    suppressWarnings(adaptive_design(500, c(250, 250), 0)$tau_list)
}

# compute MSE
mse_b<-colMeans((tau_mc_b - 2)^2)
mse_c<-colMeans((tau_mc_c - 2)^2)
# efficiency gain in terms of MSE
(mse_b-mse_c)/mse_b


