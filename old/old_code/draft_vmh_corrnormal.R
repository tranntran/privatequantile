library(mvtnorm)
library(quantreg)
library(coda)
rm(list = ls())
n = 5000
lambda = 0.1 #0.01
X = as.matrix(cbind(rep(1, n), rexp(n, lambda)))
y_init = as.matrix(4 + 5*X[,2] + rexp(n, lambda))
#R = max(abs(y_init))
R = 1
Y = y_init/R
eps = 1
tau = 0.5
logf = function(theta){
  incl_sum = ((Y - X%*%theta) <= 0)
  sumX = apply(X, 2, sum)
  ans = -eps/(4*(1-tau)*max(abs(X))) * norm(-tau*sumX + t(incl_sum)%*%X, type = "m")
  #- 0.01*theta%*%theta
  return(ans)
}

amh = function(logf, init = c(1, 1)/R, nbatch = 100000, scale, nonneg = FALSE){
  dim = length(init)
  batch = matrix(0, nrow = nbatch, ncol = dim)
  accept = matrix(0, nrow = nbatch, ncol = dim)
  batch[1, ] = init
  s = 5
  sigma = diag(scale)/s # first assume it is indep
  ct = 0
  mean_acc = NA
  s_ans = NA
  for (i in 2:nbatch){
    batch[i, ] = batch[i-1, ]
    ct_temp = 0
    for (j in 1:dim){
      propose_val = batch[i, ]
      propose_val[j] = rmvnorm(1, t(batch[i,]), sigma)[j]
      u = runif(1, 0, 1)
      check = TRUE
      #logfval[i] = logf(c(propose_val)) - logf(batch[i-1, ])
      #check = all(min(Y) <= X%*%propose_val) & all(X%*%propose_val <= max(Y))
      if ((log(u) < logf(c(propose_val)) - logf(batch[i, ])) & check){
        batch[i, ] = propose_val
        accept[i, j] = 1
        ct_temp = ct_temp + 1
        # if (ct_temp == dim) {
        #   ct = ct + 1 
        # }
        #print(ct)
      }
    }
    if (ct_temp > 0){
      # sigma = cor(na.omit(batch)) + diag(c(runif(1, 0, 0.25), runif(1, 0, 0.25)))#/100#/2000 #2000 #100
      # sigma[1, 1]= sigma[1, 1]*10
      # sigma[2, 2]= sigma[2, 2]/10
      # #sigma = sigma/max(abs(sigma))/200
      # sigma = sigma*0.00005#0.00005
      temp_acc = mean(accept[1:i,])
      if (temp_acc > 0.6) {
        s = s/2
      } else if (temp_acc < 0.3) {
        s = s*2
      }
      sigma = cov(na.omit(batch))/max(abs(y_init))/s#12
      print(sigma)
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i,]))
    s_ans = cbind(s_ans, max(abs(y_init))*s)
  }
  mean_acc = mean_acc[-1]
  s_ans = s_ans[-1]
  return(list(batch, mean(accept), mean_acc, sigma, s_ans))
}

test = rq(y_init ~ X[,2], tau = 0.5)
test
nbatch = 100000
out = amh(logf, init = coef(test), scale = c(0.1/R, 0.1/R), nbatch = nbatch)
out[[2]]
plot(out[[3]], type = "l")
plot(out[[5]], type = "l")
out[[4]]
tail(out[[1]])
out[[1]] = out[[1]]*R
tail(out[[1]])
plot(1:nbatch, out[[1]][,1], type = "l", ylab = "Value", xlab = "Runs", main = "Intercept")

plot(1:nbatch, out[[1]][,2], type = "l", ylab = "Value", xlab = "Runs", main = "Slope")
print(effectiveSize(out[[1]][,1]))
print(effectiveSize(out[[1]][,2]))
print(effectiveSize(out[[1]][(nbatch*0.2):nbatch,1]))
print(effectiveSize(out[[1]][(nbatch*0.2):nbatch,2]))
tab = rbind(c(effectiveSize(out[[1]][,1]), effectiveSize(out[[1]][,2])),
            c(effectiveSize(out[[1]][(nbatch*0.2):nbatch,1]), effectiveSize(out[[1]][(nbatch*0.2):nbatch,2])))

colnames(tab) = c("Intercept", "Slope")
rownames(tab) = c("Without Burnins", "With Burnins")
tab
acf(out[[1]][,1])
acf(out[[1]][,2])
acf(out[[1]][(nbatch*0.2):nbatch,1])


acf(out[[1]][(nbatch*0.2):nbatch,2])