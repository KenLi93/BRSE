library("sandwich")
library(R2jags)
library(dplyr)
library(xtable)
library(ggplot2)
library(ggpubr)
set.seed(4)

dat <- readRDS("nhanes_subset.rds")

n <- nrow(dat)
x <- dat[, c("MALE", "RIDAGEYR")]
y <- dat[, "BPXSY"]

## exploratory plot
plot(BPXSY ~ RIDAGEYR, data = dat, col = MALE + 1)


dat$gender <- ifelse(dat$MALE, "Male", "Female")

scatter_plot <- ggplot(dat, aes(x = RIDAGEYR, y = BPXSY, color = gender)) + 
  geom_point() +
  stat_smooth() +
  xlab("Age (years)") + ylab("Systolic blood pressure (mm Hg)") +
  theme_pubr() +
  scale_color_manual(values = c(
    "#1749FF", "#D92321", "#0AB7C9",
    "#FF6F1B", "#810094", "#378252",
    "#FF5EBF", "#3700A5", "#8F8F8F",
    "#787873"
  )) +
  theme(
    strip.text = element_text(hjust = 0.5, size = 18),
    plot.title = element_text(hjust = 0.5, size = 18),
    panel.border = element_rect(fill = NA),
    panel.grid.minor.y = element_line(),
    panel.grid.major.y = element_line(),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )

ggsave("nhanes_scatter.png", plot = scatter_plot, height = 5.3, width = 5)

## ordinary least square
lmmod <- lm(BPXSY ~ MALE + RIDAGEYR, data = dat)

## model based standard error
naive_se <- summary(lmmod)$coefficients[ ,'Std. Error']

## Huber-White robust standard error
hw_se <- sqrt(diag(vcovHC(lmmod, "HC0")))

## specify the Bayesian model
jags_model <- textConnection("model{ 
for (i in 1:n) {
  y[i] ~ dnorm(mu[i], invsigma2)
  mu[i] = x[i, ] %*% beta[] 

}

for (j in 1:3) {
  beta[j] ~ dnorm(0, 0.001)
}
	invsigma2 ~ dgamma(0.01, 0.01)
} 
")

# design matrix
X <- as.matrix(cbind(1, x))


lm_model <- jags.model(file = jags_model, 
                       data = list(y = y, x = X, n = n),
                       n.chains = 3, n.adapt = 2000)

beta_samples <- as.matrix(coda.samples(lm_model, c('beta', 'invsigma2'), n.iter = 20000)[[1]])

# compute the bayes estimates (posterior mean) and posterior variance
# the first 8,000 interations are "burn-in"
post_mean <- colMeans(beta_samples[8001:20000, ])[1:3]
pos_var_beta <- cov(beta_samples[8001:20000, ])[1:3, 1:3]


# compute the Bayesian robust standard error
Omega <- array(dim = c(3, 3, 12000))
for (i in 1:dim(Omega)[3]) {
  Omega[, , i] <- t(X) %*% diag(as.numeric((y - X %*% beta_samples[8000 + i, 1:3]) ^ 2)) %*% X %*%
    solve(t(X) %*% X) * beta_samples[8000 + i, 4]
}

Omega <- apply(Omega, c(1, 2), mean)

## the "bayesian" version of sandwich covariance matrix
post_se <- sqrt(diag(pos_var_beta))
bayes_hw_se <- sqrt(diag(pos_var_beta %*% Omega))


# tabulate the results
result_tab <- data.frame(lm = coef(lmmod),
                         naive_se = naive_se,
                         hw_se = hw_se,
                         post_mean = post_mean,
                         post_se = post_se,
                         bayes_hw_se = bayes_hw_se)
result_xtab <- xtable(result_tab, digits = 3)

print(result_xtab, file = "nhanes_result.txt")
