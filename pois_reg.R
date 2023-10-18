#!/bin/Rscript
## sandwich covariance matrix of $\mu_1, \mu_2$
library("sandwich")
library(R2jags)
library(parallel)

set.seed(4)
NSIM <- 1000
param_grid <- expand.grid(n = c(20, 50, 100),
                          a2 = seq(-0.5, 0.5, by = 0.25))

job_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))

n <- param_grid[job_id, "n"]
a2  <- param_grid[job_id, "a2"]
## read current parameters

sim_res <- mclapply(1:NSIM, function(ii) {
  ## generate covariates
  x <- runif(n, -3, 3)
  
  ## generate Gamma-Poisson distributed counts
  mu <- exp(x + a2 * x ^ 2)
  y <- rpois(n, mu)

  
  ## poisson regression model
  fit_pois <- glm(y ~ x, family = poisson)
  
  mle <- coef(fit_pois)[2]
  model_se <- sqrt(vcov(fit_pois)[2, 2])
  hw_se <- sqrt(vcovHC(fit_pois, "HC0")[2, 2])
  
  
  ## re-format the data for bayesian analysis with rjags
  ## censoring indicator
  
  jags_model <- textConnection("model{ 
for (i in 1:n) {
  y[i] ~ dpois(lambda[i])
  lambda[i] <- exp(beta[1] + x[i] * beta[2])
  
}

for (j in 1:2) {
  beta[j] ~ dnorm(0, 0.001)
}
	
} 
")
  
  pois_model <- jags.model(file = jags_model, 
                           data = list(x = x, y = y, n = n),
                           n.chains = 1, n.adapt = 1000)
  
  beta_samples <- as.matrix(coda.samples(pois_model, c('beta'), n.iter = 10000)[[1]])
  
  post_mean <- colMeans(beta_samples)[2]
  pos_var_beta <- cov(beta_samples)  ## posterior covariance matrix for beta
  
  X <- cbind(1, x)
  ## the "Bayes sandwich estimator"
  Omega <- array(dim = c(2, 2, nrow(beta_samples)))
  for (i in 1:dim(Omega)[3]) {
    Omega[, , i] <- t(X) %*% diag(as.numeric((y- exp(X %*% beta_samples[i, ])) ^ 2)) %*% X %*%
      solve(t(X) %*% diag(as.numeric(exp(X %*% beta_samples[i, ]))) %*% X)
  }
  
  Omega <- apply(Omega, c(1, 2), mean)
  
  ## the "bayesian" version of sandwich covariance matrix
  post_se <- sqrt(pos_var_beta[2, 2])
  bayes_hw_se <- sqrt((pos_var_beta %*% Omega)[2, 2])
  
  c(mle = mle, model_se = model_se, hw_se = hw_se, bayes_est = post_mean, 
    post_se = post_se, bayes_hw_se = bayes_hw_se)
}, mc.cores = detectCores())

save(sim_res, file = sprintf("./output/pois_reg_n_%s_a2_%s.RData", n, a2))



