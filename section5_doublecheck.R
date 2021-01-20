set.seed(1234)
source("functions.R")
library(ggplot2)
library(reshape2)

### Section 5.1 Results
## Distance to the truth based on epsilon allocation for each quantile
# Stepwise KNG using current data method
total_eps = 1
reps = 100
n = 5000
runs = 1000
lambda = 0.1
a0 = 4
b0 = 10
tau = c(seq(0.05, 0.95, 0.05), 0.99)
allocate_eps = seq(0.5, 0.9, 0.1)

true_beta = matrix(b0, nrow = 2, ncol = length(tau))
invert_quantile = function (x) log(1-x)/(-lambda)
true_beta[1,] = sapply(tau, invert_quantile) + a0
colnames(true_beta) = paste("q", tau*100, sep = "")
true_beta = apply(true_beta, 2, rep, reps)

dist_intercept = matrix(NA, nrow = length(allocate_eps), ncol = length(tau))
dist_slope = matrix(NA, nrow = length(allocate_eps), ncol = length(tau))

for (j in 1:length(allocate_eps)){
  median_eps = allocate_eps[j]
  scale_vec = rep(0, length(tau))
  beta_ans = NA
  for (i in 1:reps){
    print(paste("Allocate eps:", median_eps, "reps: ", i))
    x1 = rexp(n, lambda)
    x2 = a0 + b0*x1 + rexp(n, lambda)
    data = cbind(x1, x2)
    
    ans = stepwiseKNG(data = data, total_eps = total_eps, tau = tau, method = "currentdata", 
                      lower_accept = 0.2, upper_accept = 0.25,
                      median_eps = median_eps, nonneg = TRUE)
    beta_ans = rbind(beta_ans, ans[[1]])
    scale_vec = ans[[2]]
  }
  beta_ans = beta_ans[-1, ]
  beta_ans = abs(beta_ans - true_beta)
  select_intercept = rep(c(TRUE, FALSE), reps)
  select_slope = rep(c(FALSE, TRUE), reps)
  dist_intercept[j, ] = apply(beta_ans[select_intercept, ], 2, mean)
  dist_slope[j, ] = apply(beta_ans[select_slope, ], 2, mean)
}


png("results/doublecheck_stepwise_intercept_distance_data.png", width = 600, height = 400)
colnames(dist_intercept) = tau
rownames(dist_intercept) = allocate_eps
dist_intercept_data = dist_intercept
dist_intercept = melt(dist_intercept)
dist_intercept$Var1 = as.factor(dist_intercept$Var1)
colnames(dist_intercept) = c("Eps", "Quantile", "Value")
ggplot(dist_intercept, aes(x = Quantile, log(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau)  +
  ggtitle("Distance to the True Intercept by Epsilon Allocation to the Median") + 
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation")
dev.off()

png("results/doublecheck_stepwise_slope_distance_data.png", width = 600, height = 400)
colnames(dist_slope) = tau
rownames(dist_slope) = allocate_eps
dist_slope_data = dist_slope
dist_slope = melt(dist_slope)
dist_slope$Var1 = as.factor(dist_slope$Var1)
colnames(dist_slope) = c("Eps", "Quantile", "Value")
ggplot(dist_slope, aes(x = Quantile, log(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ggtitle("Distance to the True Slope by Epsilon Allocation to the Median") + 
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation")
dev.off()


# Stepwise KNG using Koenker method
total_eps = 1
reps = 100
n = 5000
runs = 1000
lambda = 0.1
a0 = 4
b0 = 10
tau = c(seq(0.05, 0.95, 0.05), 0.99)
allocate_eps = seq(0.5, 0.9, 0.1)

true_beta = matrix(b0, nrow = 2, ncol = length(tau))
invert_quantile = function (x) log(1-x)/(-lambda)
true_beta[1,] = sapply(tau, invert_quantile) + a0
colnames(true_beta) = paste("q", tau*100, sep = "")
true_beta = apply(true_beta, 2, rep, reps)

dist_intercept = matrix(NA, nrow = length(allocate_eps), ncol = length(tau))
dist_slope = matrix(NA, nrow = length(allocate_eps), ncol = length(tau))

for (j in 1:length(allocate_eps)){
  median_eps = allocate_eps[j]
  scale_vec = rep(0, length(tau))
  beta_ans = NA
  for (i in 1:reps){
    print(paste("Allocate eps:", median_eps, "reps: ", i))
    x1 = rexp(n, lambda)
    x2 = a0 + b0*x1 + rexp(n, lambda)
    data = cbind(x1, x2)
    
    ans = stepwiseKNG(data = data, total_eps = total_eps, tau = tau, method = "koenker", 
                      lower_accept = 0.2, upper_accept = 0.25, median_eps = median_eps, 
                      nonneg = TRUE)
    beta_ans = rbind(beta_ans, ans[[1]])
    scale_vec = ans[[2]]
  }
  beta_ans = beta_ans[-1, ]
  beta_ans = abs(beta_ans - true_beta)
  select_intercept = rep(c(TRUE, FALSE), reps)
  select_slope = rep(c(FALSE, TRUE), reps)
  dist_intercept[j, ] = apply(beta_ans[select_intercept, ], 2, mean)
  dist_slope[j, ] = apply(beta_ans[select_slope, ], 2, mean)
}


png("results/doublecheck_stepwise_intercept_distance_koenker.png", width = 600, height = 400)
colnames(dist_intercept) = tau
rownames(dist_intercept) = allocate_eps
dist_intercept_koenker = dist_intercept
dist_intercept = melt(dist_intercept)
dist_intercept$Var1 = as.factor(dist_intercept$Var1)
colnames(dist_intercept) = c("Eps", "Quantile", "Value")
ggplot(dist_intercept, aes(x = Quantile, log(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau)  +
  ggtitle("Distance to the True Intercept by Epsilon Allocation to the Median") + 
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation")
dev.off()

png("results/doublecheck_stepwise_slope_distance_koenker.png", width = 600, height = 400)
colnames(dist_slope) = tau
rownames(dist_slope) = allocate_eps
dist_slope_koenker = dist_slope
dist_slope = melt(dist_slope)
dist_slope$Var1 = as.factor(dist_slope$Var1)
colnames(dist_slope) = c("Eps", "Quantile", "Value")
ggplot(dist_slope, aes(x = Quantile, log(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ggtitle("Distance to the True Slope by Epsilon Allocation to the Median") + 
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation")
dev.off()
save(list = c("dist_intercept_koenker", "dist_slope_koenker", "dist_intercept_data", "dist_slope_data"),
     file = "results/doublecheck_stepwise_distance.Rdata")