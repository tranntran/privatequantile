library(doParallel)
library(snow)
library(doRNG)
source("functions.R")
num_cores=detectCores()-1


workers=makeCluster(num_cores,type="SOCK",outfile="log.txt")
registerDoParallel(workers)

#registerDoRNG(1)

# Stepwise KNG using current data method
total_eps = 1
reps = 2 #remember to change reps
n = 5000
runs = 1000
lambda = 0.1
a0 = 4
b0 = 10
tau = c(seq(0.05, 0.95, 0.05), 0.99)
allocate_eps = seq(0.1, 0.9, 0.1)

true_beta = matrix(b0, nrow = 2, ncol = length(tau))
invert_quantile = function (x) log(1-x)/(-lambda)
true_beta[1,] = sapply(tau, invert_quantile) + a0
colnames(true_beta) = paste("q", tau*100, sep = "")
true_beta = apply(true_beta, 2, rep, reps)

distance_median_eps_stepd = function(median_eps){
  dist_intercept = matrix(NA, nrow = length(allocate_eps), ncol = length(tau))
  dist_slope = matrix(NA, nrow = length(allocate_eps), ncol = length(tau))
  scale_vec = rep(0, length(tau))
  beta_ans = matrix(NA, nrow = 2*reps, ncol = length(tau))
  # for each median eps allocation, scale is calculated once and can be reused
  for (i in 1:reps){
    print(paste("Allocate eps:", median_eps, "reps: ", i))
    x1 = rexp(n, lambda)
    x2 = a0 + b0*x1 + rexp(n, lambda)
    data = cbind(x1, x2)
    
    ans = stepwiseKNG(data = data, total_eps = total_eps, tau = tau, 
                      method = "currentdata", lower_accept = 0.2, 
                      upper_accept = 0.25, median_eps = median_eps, 
                      nonneg = TRUE, scale = scale_vec)
    beta_ans[c(i*2-1, i*2), ] = ans[[1]]
    scale_vec = ans[[2]]
  }
  
  beta_ans_is = abs(beta_ans - true_beta)
  select_intercept = rep(c(TRUE, FALSE), reps)
  select_slope = rep(c(FALSE, TRUE), reps)
  dist_intercept = apply(beta_ans_is[select_intercept, ], 2, mean)
  dist_slope = apply(beta_ans_is[select_slope, ], 2, mean)
  sd_intercept = apply(beta_ans_is[select_intercept, ], 2, sd)
  sd_slope = apply(beta_ans_is[select_slope, ], 2, sd)
  
  beta_ans_l2 = (beta_ans - true_beta)^2
  beta_ans_l2 = rowsum(beta_ans_l2, rep(1:reps, each = 2))
  beta_ans_l2 = t(apply(beta_ans_l2, 1, sqrt))
  dist_l2 = apply(beta_ans_l2, 2, mean)
  sd_l2 = apply(beta_ans_l2, 2, sd)
  
  return(list(dist_intercept, dist_slope, dist_l2, sd_intercept, sd_slope, sd_l2))
  
}

ctype = rbind

# Out = foreach(median_eps = allocate_eps,.combine=ctype, .errorhandling='stop') %dopar% 
#   distance_median_eps_stepd(median_eps)

Out1 = foreach(median_eps = allocate_eps,.combine=ctype, .errorhandling='stop',
              .options.RNG=1) %dorng% rexp(5, 1)
load("out.Rdata")

foreach(median_eps = allocate_eps,.combine=ctype, .errorhandling='stop',
        .options.RNG=1) %dopar% {
          set.seed(123)
          rexp(5, 1)
        }

stepwise_data_int = matrix(unlist(Out[,1]), nrow = 9, byrow = TRUE)
stepwise_data_slope = matrix(unlist(Out[,2]), nrow = 9, byrow = TRUE)
stepwise_data_l2 = matrix(unlist(Out[,3]), nrow = 9, byrow = TRUE)

sd_int = matrix(unlist(Out[,4]), nrow = 9, byrow = TRUE)
sd_slope = matrix(unlist(Out[,5]), nrow = 9, byrow = TRUE)
sd_l2 = matrix(unlist(Out[,6]), nrow = 9, byrow = TRUE)

load("../output/testseed1.Rdata")
sd_l2 == sd_l21
beta_ans

save(list = c("stepwise_data_int1", "stepwise_data_slope1", "stepwise_data_l21",
              "sd_int1", "sd_slope1", "sd_l21"), file = "../output/testseed1.Rdata")

stopCluster(workers)
