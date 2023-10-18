library(dplyr)
library(RColorBrewer)
cols <- brewer.pal(8, "Set1")
NSIM <- 10000

set.seed(1127)

mu <- 0; eta <- 1000 ## hyperparameters

param_grid <- expand.grid(N = seq(5, 50, 5),
                          sigma2 = c(0.1, 10),
                          credint_coverage = NA,
                          credint_avelen = NA,
                          rci_coverage = NA,
                          rci_avelen = NA,
                          brci_coverage = NA,
                          brci_avelen = NA,
                          hw_se_se = NA,
                          bayes_hw_se_se = NA)

for (ii in 1:nrow(param_grid)) {
  N <- param_grid[ii, "N"]
  sigma2 <- param_grid[ii, "sigma2"]
  
  sim_res <- replicate(NSIM, {
    X <- rnorm(N, 0, sqrt(sigma2))
    credint <- qnorm(c(0.025, 0.975), (mu + eta ^ 2 * sum(X)) / (N * eta ^ 2 + 1), sqrt(eta ^ 2 / (N * eta ^ 2 + 1)))
    hw_se <- sqrt(sum((X - mean(X)) ^ 2) / N ^ 2)
    rci <- mean(X) + qnorm(c(0.025, 0.975)) * hw_se
    bayes_hw_se <- sqrt(eta ^ 2 * sum((X - mean(X)) ^ 2) / (N * (N * eta ^ 2 + 1)) + 
                          (eta ^ 2 / (N * eta ^ 2 + 1)) ^ 2 + (eta ^ 2 * (mu - mean(X)) ^ 2) / (N * eta ^ 2 + 1) ^ 3)
    brci <- (mu + eta ^ 2 * sum(X)) / (N * eta ^ 2 + 1) + qnorm(c(0.025, 0.975)) * bayes_hw_se
      
    c(credint_cover = as.numeric(credint[1] < 0 & credint[2] > 0),
      credint_len = credint[2] - credint[1],
      rci_cover = as.numeric(rci[1] < 0 & rci[2] > 0),
      rci_len = rci[2] - rci[1],
      brci_cover = as.numeric(brci[1] < 0 & brci[2] > 0),
      brci_len = brci[2] - brci[1],
      hw_se = hw_se, bayes_hw_se = bayes_hw_se)
  })
  
  param_grid[ii, c(3:8)] <- rowMeans(sim_res[1:6, ])
  param_grid[ii, 9] <- sd(sim_res[7,])
  param_grid[ii, 10] <- sd(sim_res[8,])
}

param_grid <- arrange(param_grid, sigma2, N)

png("normal_mean_cis.png", width = 720, height = 720)
par(mfrow = c(2, 2), mar = c(5, 4.5, 2, 1))

for (ss in c(0.1, 10)) {
  plot.new()
  plot.window(xlim = c(5, 50), ylim = c(0.8, 1.01))
  axis(1, at = seq(10, 50, 10), pos = 0.8, cex.axis = 2)
  axis(2, at = seq(0.8, 1, by = 0.05), pos = 5, cex.axis = 1.8)
  title(bquote(sigma^2 ~ " = " ~ .(ss)), cex.main= 2,
        xlab = "n", ylab = "Coverage Probability",
        cex.lab = 1.8)
  
  rect(5, 0.8, 50, 1.01)
  
  for (yy in seq(0.825, 1, by = 0.025)) {
    lines(x = c(5, 50), y = c(yy, yy), col = "gray", lty = 2)
  }
  lines(credint_coverage ~ N, data = subset(param_grid, sigma2 == ss), col = cols[1], lwd = 2,
       type = "l")
  
  lines(rci_coverage ~ N, data = subset(param_grid, sigma2 == ss), col = cols[2], lwd = 2)
  lines(brci_coverage ~ N, data = subset(param_grid, sigma2 == ss), col = cols[3], lwd = 2)
  
  lines(x = c(5, 50), y = c(0.95, 0.95), col = "red", lty = 2, lwd = 2)
}

par(mar = c(5, 4.5, 0.5, 1))
for (ss in c(0.1, 10)) {
  plot.new()
  plot.window(xlim = c(5, 50), ylim = c(0, 7.5))
  axis(1, at = seq(10, 50, 10), pos = 0, cex.axis = 2)
  axis(2, at = seq(0, 7.5, by = 2.5), pos = 5, cex.axis = 1.8)
  title(xlab = "n", ylab = "Average Interval Length",
        cex.lab = 1.8)
  
  rect(5, 0, 50, 7.5)
  
  for (yy in seq(1.25, 6.25, by = 1.25)) {
    lines(x = c(5, 50), y = c(yy, yy), col = "gray", lty = 2)
  }
  lines(credint_avelen ~ N, data = subset(param_grid, sigma2 == ss), col = cols[1], lwd = 2,
        type = "l")
  
  lines(rci_avelen ~ N, data = subset(param_grid, sigma2 == ss), col = cols[2], lwd = 2)
  lines(brci_avelen ~ N, data = subset(param_grid, sigma2 == ss), col = cols[3], lwd = 2)
}
dev.off()




## plot precision for sigma2 = 0.1
ss <- 0.1
png(paste0("normal_mean_se_precision_sigma2_", ss * 10, ".png"), 
    width = 480, height = 480) 
par(mar = c(5, 4.5, 2, 1))
plot.new()
plot.window(xlim = c(5, 50), ylim = c(0, 0.05))
axis(1, at = seq(10, 50, 10), pos = 0.8, cex.axis = 2)
axis(2, at = seq(0, 0.05, by = 0.01), pos = 5, cex.axis = 1.8)
title(bquote(sigma^2 ~ " = " ~ .(ss)), cex.main= 2,
      xlab = "n", ylab = "Standard Error of estimated S.E.",
      cex.lab = 1.8)

rect(5, 0, 50, 0.051)

for (yy in seq(0, 0.05, 0.005)) {
  lines(x = c(5, 50), y = c(yy, yy), col = "gray", lty = 2)
}
lines(hw_se_se ~ N, data = subset(param_grid, sigma2 == ss), col = cols[2], lwd = 2,
      type = "l")

lines(bayes_hw_se_se ~ N, data = subset(param_grid, sigma2 == ss), col = cols[3], lwd = 2)
dev.off()



## plot precision for sigma2 = 10
  ss <- 10
  png(paste0("normal_mean_se_precision_sigma2_", ss * 10, ".png"), 
      width = 480, height = 480) 
  par(mar = c(5, 4.5, 2, 1))
  plot.new()
  plot.window(xlim = c(5, 50), ylim = c(0, 0.5))
  axis(1, at = seq(10, 50, 10), pos = 0.8, cex.axis = 2)
  axis(2, at = seq(0, 0.5, by = 0.1), pos = 5, cex.axis = 1.8)
  title(bquote(sigma^2 ~ " = " ~ .(ss)), cex.main= 2,
        xlab = "n", ylab = "Standard Error of estimated S.E.",
        cex.lab = 1.8)
  
  rect(5, 0, 50, 0.51)
  
  for (yy in seq(0, 0.5, 0.05)) {
    lines(x = c(5, 50), y = c(yy, yy), col = "gray", lty = 2)
  }
  lines(hw_se_se ~ N, data = subset(param_grid, sigma2 == ss), col = cols[2], lwd = 2,
        type = "l")
  
  lines(bayes_hw_se_se ~ N, data = subset(param_grid, sigma2 == ss), col = cols[3], lwd = 2)
  dev.off()
}




pdf("ci_legend.pdf", width = 2.5, height = 2.5)
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c("Cred. Int.", "RCI", "BRCI"), fill = cols[1:3], bty = "n", cex = 2)
dev.off()
