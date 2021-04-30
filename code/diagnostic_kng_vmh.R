library(quantreg)
# Grid Sampling
# use grid sampling to observe the behavior of KNG when range of data (X) changes
# distribution of intercept and slope from grid sampling seems to match distribution
# from mcmc

# generate the data
rm(list = ls())
set.seed(525252)
n = 5000
runs = 1000
lambda = 0.1
a0 = 4
b0 = 2
a1 = 3 
b1 = 2
b2 = 1
# when lambda >= 1, KNG works ok; but when lambda reduces to 0.01 KNG has more 
# variability in its estimate
x1 = rexp(n, lambda)
# when range of runif is up to 10, it works ok; when range increases to 100 or more
# it seems to be harder for kng to get to the true estimates
#x1 = runif(n, 0, 100)
x2 = a0 + b0*x1 + rexp(n, lambda)


# prepare data for kng
Y = as.matrix(x2)
R = max(abs(Y))
#R = 1
Y = Y/R
# Some observations based on results of grid sampling
# When x1 ~ unif, normalize it will result in much more accuracy and less variability
# When x1 ~ exp, take log then normalize it will produce more stable estimates

#X = as.matrix(cbind(rep(1, n), log(x1)))
X = as.matrix(cbind(rep(1, n), log(x1)/max(abs(log(x1)))))
#X = as.matrix(cbind(rep(1, n), x1/max(abs(x1))))
sumX = apply(X = X, 2, FUN = sum)

nbatch = 100000 # number of runs in mcmc
tau = 0.5
ep = 0.1
rq(x2 ~ log(x1), tau)
#rq(x2 ~ x1, tau)


# sample from this function
logA = function(beta){
  #print(beta)
  left = cbind(Y, X)%*% c(1, -beta)
  lessEq = (left <= 0)
  ans = -(ep/2) * norm(-tau*sumX + t(X)%*%lessEq, type = "M") / ((1-tau)*2*max(abs(X))) -
    1/2*beta%*%beta
  return(ans)
}


# grid sampling
grid_length = 1000
beta = cbind(seq(-500, 500, length.out = grid_length), seq(-500, 500, length.out = grid_length))
beta = beta/R
beta_grid = expand.grid(x = beta[,1], y = beta[,2])


# evaluate the log of the target density at each of 1000^2 pairs of slopes and intercepts on the grid.
log_target = matrix(sapply(1 : nrow(beta_grid), 
                           function(k){logA(c(beta_grid[k, 1], beta_grid[k, 2]))}), 
                    nrow = grid_length, ncol = grid_length)

# back to the likelihood values (from log-likelihood) in a way that the largest likelihood value is 1.
target_z = exp(log_target - max(log_target))

# deriving marginal of intercept
beta0_marg = rowSums(target_z)

# sample from marginal "with replacement"
beta0_index = sample(1:grid_length, size = 1e3, prob = beta0_marg, replace = TRUE)
beta0_sample = beta[,1][beta0_index]

# sample from conditional given each sampled value of intercepts
beta1_sample = sapply(beta0_index, function(i){
  sample(beta[,2], size = 1, prob = target_z[i, ])
})


# jitter the samples to blur the grid effect
beta0_sample = jitter(beta0_sample, amount = 0.00001)
beta1_sample = jitter(beta1_sample, amount = 0.00001)


# scatter plot of the samples on the contour plot (true distribution)
contour(beta[, 1]*R, beta[, 2]*R/max(abs(log(x1))), target_z, nlevels = 5, xlab = "", ylab = "", 
        xlim = c(-10, 20), ylim = c(-10, 20))
points(beta0_sample*R, beta1_sample*R/max(abs(log(x1))), pch = 46, cex = 2)
mtext(side = 1, text = "Intercept", line = 2, cex = 1.3)
mtext(side = 2, text = "Slope", line = 1.9, cex = 1.2)

# histogram of x with the true marginal density superimposed
hist(beta0_sample*R, 20, prob = TRUE, main = "Histogram of Intercepts")

# histogram of y with the true marginal density superimposed
hist(beta1_sample*R/max(abs(log(x1))), 20, prob = TRUE, main = "Histogram of Slopes")

quantile(beta0_sample*R, c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(beta1_sample*R/max(abs(log(x1))), c(0.1, 0.25, 0.5, 0.75, 0.90))


####################################################################################
# MCMC
rm(list = ls())
set.seed(1)
n = 5000
runs = 1000
lambda = 0.1
a0 = 4
b0 = 2
a1 = 3 
b1 = 2
b2 = 1
x1 = rexp(n, lambda)
x2 = a0 + b0*x1 + rexp(n, lambda)
data = cbind(x1, x2)


# prepare data for kng
Y = as.matrix(x2)
R = max(abs(Y))
Y = Y/R
X = as.matrix(cbind(rep(1, n), log(x1)/max(abs(log(x1)))))
sumX = apply(X = X, 2, FUN = sum)

nbatch = 100000 # number of runs in mcmc
tau = 0.5
ep = 0.1
m = 2 # m stands for dimension
init = rep(9, m)/R # set starting value
rq(x2 ~ log(x1), tau)


logA = function(beta){
  left = cbind(Y, X)%*% c(1, -beta)
  lessEq = (left <= 0)
  ans = -(ep/2) * norm(-tau*sumX + t(X)%*%lessEq, type = "M") / ((1-tau)*2*max(abs(X))) -
    1/2*beta%*%beta
  return(ans)
}


metrop = function(logA, init, Y, X, nbatch = 10000, scale, nonneg = FALSE){
  dim = length(init)
  
  out = list(accept = 0, batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim))
  U = matrix(runif(nbatch*dim, min = 0, max = 1), nrow = nbatch, ncol = dim)
  Prop = matrix(rnorm(nbatch*dim, m = 0, s = scale), nrow = nbatch, ncol = dim)
  # scale_intercept = 1
  # scale_slope = 0.01
  # Prop = cbind(rnorm(nbatch, 0, scale_intercept), rnorm(nbatch, 0, scale_slope))
  
  for(r in 1:nbatch){
    if(r == 1){
      out$batch[1,] = init
      oldLogA = logA(init)
    } else {
      out$batch[r,] = out$batch[r-1,]
    }
    ###
    for(i in 1:dim){
      oldVector = out$batch[r,]
      
      newVector = oldVector
      newVector[i] = newVector[i] + Prop[r,i]
      newLogA = logA(newVector)
      check = TRUE
      #check = all(quantile(Y, 0.01) <= X%*%newVector) & all(X%*%newVector <= quantile(Y, 0.99))
      #check = all(min(Y) <= X%*%newVector) & all(X%*%newVector <= max(Y))
      # nonnegative constraints can be replace using range of Y constraints ^
      # if (nonneg) {
      #   check = check & (newVector[i] > 0) 
      # }
      
      #TestValue = exp((newLogA - oldLogA))
      if(log(U[r,i]) < newLogA - oldLogA & check){
        out$batch[r,] = newVector
        out$accept = out$accept + 1/(nbatch*dim)
        oldLogA = newLogA
      }
    }
  }
  
  return(out)
}

library(coda)
par(mfrow = c(1, 2))
out = metrop(logA = logA, init = init, X = X, Y = Y, nbatch = 100000, scale = 2/R)
out$accept*100
plot(1:100000, out$batch[,1]*R, type = "l", xlab = "runs", ylab = "Intercept")
plot(1:100000, out$batch[,2]*R/max(abs(log(x1))), type = "l", xlab = "runs", ylab = "Slope")

#mcmc diagnostic looks a bit concerning
acf(out$batch[,1])
acf(out$batch[,2])
effectiveSize(out$batch[,1])
effectiveSize(out$batch[,2])
