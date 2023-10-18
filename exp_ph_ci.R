#!/bin/Rscript
## sandwich covariance matrix of $\mu_1, \mu_2$
library("sandwich")
library(R2jags)
library(parallel)
library(survival)
set.seed(4)
NSIM <- 1000
param_grid <- expand.grid(n = seq(10, 50, by = 5),
                          alpha = c(0.5, 1, 2),
                          beta = c(0, 1, 2))

job_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

n <- param_grid[job_id, "n"]
alpha <- param_grid[job_id, "alpha"]
beta <- param_grid[job_id, "beta"]
## read current parameters

sim_res <- mclapply(1:NSIM, function(ii) {
  ## generate covariates
  zz <- rnorm(n)
  
  ## parametrization according to Weibull PH model
  bb <- exp(-zz * beta / alpha)
  
  ## generate survival time
  ss <- rweibull(n, shape = alpha, scale = bb)
  ss <- ceiling(ss * 100) / 100 ## make sure ss >= 0.01
  
  ## censoring time -- administrative censoring
  cc <- 10
  
  ## event indicator
  delta <- as.numeric(ss <= cc)

  
  ## observed failure time
  tt <- ss * delta + rep(cc, n) * (1 - delta) 
  
  df_surv <- data.frame(tt = tt, delta = delta, zz = zz)
  
  ## exponential proportional hazards model
  fit_exp <- survreg(Surv(time = tt, event = delta) ~ 1 + zz, dist = "exponential",
                     data = df_surv)
  
  mle <- -coef(fit_exp)[2]
  model_se <- sqrt(vcov(fit_exp)[2, 2])
  hw_se <- sqrt(vcovHC(fit_exp, "HC0")[2, 2])
  
  
  ## re-format the data for bayesian analysis with rjags
  ## censoring indicator
  IsCensored <- 1 - delta
  surv_time <- tt
  surv_time[IsCensored == 1] <- NA
  
  jags_model <- textConnection("model{ 
for (i in 1:n) {
  IsCensored[i] ~ dinterval(t[i], 10)
  t[i] ~ dexp(lambda[i])
  lambda[i] <- exp(beta[1] + z[i] * beta[2])
  
}

for (j in 1:2) {
  beta[j] ~ dnorm(0, 0.001)
}
	
} 
")
  
  exp_model <- jags.model(file = jags_model, 
                         data = list(IsCensored = IsCensored, t = surv_time, n = n, z = zz),
                         n.chains = 1, n.adapt = 1000)
  
  beta_samples <- as.matrix(coda.samples(exp_model, c('beta'), n.iter = 10000)[[1]])
  
  post_mean <- colMeans(beta_samples)[2]
  pos_var_beta <- cov(beta_samples)  ## posterior covariance matrix for beta
  
  Z <- cbind(1, zz)
  ## the "Bayes sandwich estimator"
  Omega <- array(dim = c(2, 2, nrow(beta_samples)))
  for (i in 1:dim(Omega)[3]) {
    Omega[, , i] <- t(Z) %*% diag(as.numeric((delta - tt * exp(Z %*% beta_samples[i, ])) ^ 2)) %*% Z %*%
      solve(t(Z) %*% diag(as.numeric(tt * exp(Z %*% beta_samples[i, ]))) %*% Z)
  }
  
  Omega <- apply(Omega, c(1, 2), mean)
  
  ## the "bayesian" version of sandwich covariance matrix
  post_se <- sqrt(pos_var_beta[2, 2])
  bayes_hw_se <- sqrt((pos_var_beta %*% Omega)[2, 2])
  
  c(mle = mle, model_se = model_se, hw_se = hw_se, bayes_est = post_mean, 
    post_se = post_se, bayes_hw_se = bayes_hw_se)
}, mc.cores = detectCores())

save(sim_res, file = sprintf("./output/exp_ref_n_%s_alpha_%s_beta_%s.RData", n, alpha, beta))



