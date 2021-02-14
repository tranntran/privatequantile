source("code/functions.R")

total_eps = 1
reps = 1000 #change reps
n = 5000
runs = 1000
lambda = 0.1
a0 = 4
b0 = 10
#tau = c(seq(0.05, 0.95, 0.05), 0.99)
main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95)
tau = seq(0.05, 0.95, 0.05)

true_beta = matrix(b0, nrow = 2, ncol = length(tau))
invert_quantile = function (x) log(1-x)/(-lambda)
true_beta[1,] = sapply(tau, invert_quantile) + a0
colnames(true_beta) = paste("q", tau*100, sep = "")

x1 = rexp(n, lambda)
x2 = a0 + b0*x1 + rexp(n, lambda)
data = cbind(x1, x2)

scale_vec = rep(0, length(tau))

ans9
stepwiseKNG(data, total_eps, median_eps = 0.3, tau = main_tau, nonneg = TRUE,
            method = "currentdata")



utility_logit = function(data, inds){
  data = cbind(data, inds)
  mod = glm(inds ~ ., data = data, family = binomial(link = "logit"))
  preds = predict(mod, type = "response")
  score = sum((preds-0.5)^2)/nrow(data)
  return(score)
}

reps = 100
utility = matrix(NA, nrow = reps, ncol = 2)
for (i in 1:reps){
  ans3 = sandwichKNG(data = data, total_eps = 0.1, tau = tau, 
                     main_tau = main_tau, median_eps = 0.7, scale = scale_vec, 
                     method = "currentdata", lower_accept = 0.2, upper_accept = 0.25,
                     main_tau_eps = 0.3, nonneg = TRUE)[[1]]
  ans9 = sandwichKNG(data = data, total_eps = 0.1, tau = tau, 
                     main_tau = main_tau, median_eps = 0.7, scale = scale_vec, 
                     method = "currentdata", lower_accept = 0.2, upper_accept = 0.25,
                     main_tau_eps = 0.9, nonneg = TRUE)[[1]]
  sample = sample(1:length(tau), 5000, TRUE)
  all3 = cbind(rep(1, n), x1) %*%ans3
  all9 = cbind(rep(1, n), x1) %*%ans9
  syn3 = all3[cbind(1:n, sample)]
  syn9 = all9[cbind(1:n, sample)]
  
  inds = c(rep(1, n), rep(0, n))
  syndata3 = as.data.frame(rbind(data, cbind(x1, syn3)))
  syndata9 = as.data.frame(rbind(data, cbind(x1, syn3)))
  
  utility[i, 1] = utility_logit(syndata3, inds)
  utility[i, 2] = utility_logit(syndata9, inds)
  
  plot(density(x2))
  lines(density(syn3), col = "red")
  lines(density(syn9), col = "blue")
}
apply(utility, 2, mean)
