# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
# example(stan_model, package = "rstan", run.dontrun = TRUE)
library("rstan")
library(coda)
library(quantreg)


rm(list = ls())
#set.seed(525252)


n = 5000
lambda = 0.1
a0 = 4
b0 = 2
a1 = 3 
b1 = 2
b2 = 1
x1 = rexp(n, lambda)
#x1 = runif(n, -100, 100)
x2 = a0 + b0*x1 + rexp(n, lambda)
X = cbind(rep(1, n), x1)
sumX = apply(X, 2, sum)
fit_np = rq(x2~ x1, tau = 0.5)
#R = max(abs(x2))
#divided by max(Y)
R = 1

data = list(N = n, dim = 2, eps = 1, tau = 0.5, y = as.matrix(x2), x = as.matrix(X), 
            sumX = t(as.matrix(sumX)))
fit = stan(file = 'kng.stan', data = data, warmup = 1000, iter = 10000, chains = 4, 
           control = list(max_treedepth = 10), cores = 5)
#, init = list(beta = as.matrix(c(10,2)))
#http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
#http://mc-stan.org/misc/warnings.html#bulk-ess 
print(fit)
plot(fit)
#fit1 = fit
#pairs(fit, pars = c("mu", "tau", "lp__"))
save(fit, file = "stan_heavytailed_depth10.Rdata")
#load("stan_heavytailed_depth12.Rdata")
la <- extract(fit, permuted = TRUE) # return a list of arrays 
pairs(fit, pars = c("beta", "lp__"), las = 1)
beta <- matrix(la$beta, ncol = 2)
#beta <- beta*R
acf(beta[,1])
acf(beta[,2])
effectiveSize(beta[,1])
effectiveSize(beta[,2])



plot(density(x2), ylim = c(0, 0.04), main = "Density of Predicted Values Y")
lines(density(1+2.66*x1), col = "red")
lines(density(predict(fit_np)), col = "blue")












# fit <- stan(file = stan_model1, data = stan_data, warmup = 500, iter = 1000, chains = 4, cores = 2, thin = 1)
### return an array of three dimensions: iterations, chains, parameters 
# a <- extract(fit, permuted = FALSE) 
#fit_test = stan(file = 'normal.stan', data = list(mu = 0, sigma2 = 1))
# schools_dat <- list(J = 8, 
#                     y = c(28,  8, -3,  7, -1,  1, 18, 12),
#                     sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
### use S3 functions on stanfit objects
# a2 <- as.array(fit)
# m <- as.matrix(fit)
# d <- as.data.frame(fit)