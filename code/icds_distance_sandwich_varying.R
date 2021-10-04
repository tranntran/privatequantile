# This code compare the performance of the sandwich varying function based on
# the level of privacy budget used on the anchor quantiles vs other quantiles.
# The code will compute the mean and std of the distance to the true parameters,
# for each amount of allocation 10%, 20%,..., 90% of the total budget.

library(doParallel)
library(snow)
library(doRNG)
source("functions_final.R") #change directory
num_cores=detectCores()-1 #use all available core but 1

workers=makeCluster(num_cores,type="SOCK",outfile="log.txt")
registerDoParallel(workers)


total_eps = 0.5
reps = 100 #change reps
n = 5000
runs = 10000
lambda = 0.1
a0 = 4
b0 = 3
tau = c(seq(0.05, 0.95, 0.05), 0.99)
main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
allocate_eps = c(1/length(tau), seq(0.1, 0.9, 0.1))

true_beta = matrix(b0, nrow = 2, ncol = length(tau))
invert_quantile = function (x) log(1-x)/(-lambda)
true_beta[1,] = sapply(tau, invert_quantile) + a0
colnames(true_beta) = paste("q", tau*100, sep = "")
true_beta = apply(true_beta, 2, rep, reps)

# Privacy budget allocated to the median is 0.7, based on the result of file
# icds_distance_stepwise_data.R.
distance_maintau_eps = function(main_tau_eps){
  scale_vec = rep(0, length(tau))
  beta_ans = matrix(NA, nrow = 2*reps, ncol = length(tau))
  for (i in 1:reps){
    print(paste("Allocate eps:", main_tau_eps, "reps: ", i))
    x1 = rexp(n, lambda)
    x2 = a0 + b0*x1 + rexp(n, lambda)
    data = cbind(x1, x2)
    mod = "x2 ~ x1"
    
    ans = sandwichKNG(data = data, total_eps = total_eps, median_eps = 0.8, #fill in 
                main_tau_eps = main_tau_eps, tau = tau, main_tau = main_tau, 
                scale = 0.005, sw_scale = 0.001, nbatch = runs, method = "varying", 
                ub = 1000, lb = 0, formula = mod)
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

t = 123
Out = foreach(main_tau_eps = allocate_eps,.combine=ctype, .errorhandling='stop',
              .options.RNG = t) %dorng% distance_maintau_eps(main_tau_eps)


sandwich_varying_int = matrix(unlist(Out[,1]), nrow = 10, byrow = TRUE)
sandwich_varying_slope = matrix(unlist(Out[,2]), nrow = 10, byrow = TRUE)
sandwich_varying_l2 = matrix(unlist(Out[,3]), nrow = 10, byrow = TRUE)

sd_int = matrix(unlist(Out[,4]), nrow = 10, byrow = TRUE)
sd_slope = matrix(unlist(Out[,5]), nrow = 10, byrow = TRUE)
sd_l2 = matrix(unlist(Out[,6]), nrow = 10, byrow = TRUE)

filename = paste("../output/sandwich_varying_sd_", t, "_e", total_eps,".Rdata", sep = "")

save(list = c("sandwich_varying_int", "sandwich_varying_slope", "sandwich_varying_l2",
              "sd_int", "sd_slope", "sd_l2"), file = filename)

stopCluster(workers)
