n = 5000
lambda = 0.1
X = as.matrix(cbind(rep(1, n), rexp(n, lambda)))
y_init = as.matrix(4 + 10*X[,2] + rexp(n, lambda))
Y = y_init/max(abs(y_init))
theta = c(10, 10)/max(Y)
eps = 1
tau = 0.5
logf = function(theta){
  incl_sum = ((Y - X%*%theta) <= 0)
  sumX = apply(X, 2, sum)
  ans = -eps/(4*(1-tau)*max(X)) * norm(-tau*sumX + t(incl_sum)%*%X, type = "m") 
  #- 1/2*theta%*%theta
  return(ans)
}

vmh = function(logf, init = rep(0, ncol(X)), nbatch = 10000, scale, nonneg = FALSE){
  dim = length(init)
  batch = matrix(0, nrow = nbatch, ncol = dim)
  batch[1, ] = init
  for (i in 2:nbatch){
    batch[i, ] = batch[i-1, ]
    for (j in 1:dim){
      propose_val = batch[i, ]
      propose_val[j] = rnorm(1, batch[i-1, j], scale[j])
      u = runif(1, 0, 1)
      if (log(u) < logf(propose_val) - logf(batch[i, ])){
        batch[i, ] = propose_val
      }
    }
  }
  return(batch)
}

reps = 1000
output = matrix(NA, ncol = 2, nrow = reps)
for (i in 1:reps){
  print(i)
  out = vmh(logf, scale = c(1/max(y_init), 1/max(y_init)))
  out = out*max(y_init)
  plot(1:10000, out[,1], type = "l", xlab = "Value", ylab = "Runs", main = "Intercept")
  plot(1:10000, out[,2], type = "l", xlab = "Value", ylab = "Runs", main = "Slope")
  output[i, ] = out[10000,]
  
}
hist(output[,1], breaks = 25, xlab = "Value", ylim = c(0, 600), 
     main = paste("Histogram of Intercept over", reps, "Runs \n with VMH"))
hist(output[,2], breaks = 25, xlab = "Value", ylim = c(0, 650),
     main = paste("Histogram of Slope over", reps, "Runs \n with VMH"))

apply(output, 2, hist, breaks = 50)
#par(mfrow = c(1, 2))


#library(coda)
heidel.diag(out)
geweke.diag(out)
