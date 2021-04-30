#https://m-clark.github.io/docs/ld_mcmc/index_onepage.html#am
library(mvtnorm)
library(quantreg)
library(coda)
library(matrixcalc)
library(corpcor)

x1 = rexp(n, lambda)
x2 = 4 + 3*x1 + rexp(n, lambda)
data = cbind(x1, x2)

rm(list = ls())
n = 5000
lambda = 0.1 #0.01
X = as.matrix(cbind(rep(1, n), rexp(n, lambda)))
#M = max(X[,2])
#X[,2] = X[,2]/M
y_init = as.matrix(4 + 3*X[,2] + rexp(n, lambda))
#R = max(abs(y_init))
R = 1
Y = y_init/R
#theta = c(10, 10)/max(Y)
eps = 1
tau = 0.5
logf = function(theta){
  incl_sum = ((Y - X%*%theta) <= 0)
  sumX = apply(X, 2, sum)
  ans = -eps/(4*(1-tau)*max(abs(X))) * norm(-tau*sumX + t(incl_sum)%*%X, type = "m") #- 0.00001*theta%*%theta
  return(ans)
}

amh = function(logf, init = c(1, 1)/R, nbatch = 100000, scale, nonneg = FALSE){
  logfval = rep(NA, nrow = 100000)
  dim = length(init)
  batch = matrix(0, nrow = nbatch, ncol = dim)
  accept = matrix(0, nrow = nbatch)
  batch[1, ] = init
  #s = #1/20000 #5 # first assume it is indep
  ct = 0
  mean_acc = NA
  s_ans = NA
    scale_sigma = 1/20000
  sigma = diag(scale)*scale_sigma
  for (i in 2:nbatch){
    propose_val = rmvnorm(1, t(batch[i-1,]), sigma)
    # propose_val = rep(0, dim)
    # for (j in 1:dim){
    #   propose_val[j] = rnorm(1, batch[i-1, j], scale[j])
    # }
    u = runif(1, 0, 1)
    check = TRUE
    logfval[i] = logf(c(propose_val)) - logf(batch[i-1, ])
    check = all(0<= X%*%t(propose_val)) #& all(X%*%t(propose_val) <= max(Y))
    if ((log(u) < logf(c(propose_val)) - logf(batch[i-1, ])) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
      #print(ct)
    } else {
      batch[i, ] = batch[i-1, ]
    }
    if (ct == 10){
      sigma = cor(na.omit(batch)) + diag(dim)*0.00002 #diag(c(runif(1, 0, 0.25), runif(1, 0, 0.25)))#/100#/2000 #2000 #100
      sigma[1, 1]= sigma[1, 1]#/2
      sigma[1, 2] = sigma[1, 2]#*sqrt(2) #7
      sigma[2, 1] = sigma[2, 1]#*sqrt(2)
      sigma[2, 2]= sigma[2, 2]#*2
      #sigma = sigma*0.1
      temp_acc = mean(accept[1:i])
      if (temp_acc < 0.2){ #prev 0.1 and 0.2
        scale_sigma = scale_sigma/2#*0.1
      } else if (temp_acc > 0.5) {
        print(temp_acc)
        scale_sigma = scale_sigma*2
      }
      print(is.positive.semi.definite(sigma))
      sigma = sigma*scale_sigma #+ diag(ncol(sigma))*scale_sigma*0.01
      # if (!is.positive.semi.definite(sigma)){
      #   sigma = make.positive.definite(sigma, tol = 1e-10)
      # }
      # #sigma = sigma/max(abs(sigma))/200
      # sigma = sigma*0.00005#0.00005
      # temp_acc = mean(accept[1:i])
      # if (temp_acc > 0.4) {
      #   s = s/2
      # } else if (temp_acc < 0.2) {
      #   s = s*2
      # }
      # sigma = cov(na.omit(batch))/max(abs(y_init))/s#12
      print(sigma)
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    s_ans = cbind(s_ans, scale_sigma)
  }
  mean_acc = mean_acc[-1]
  s_ans = s_ans[-1]
  return(list(batch, mean(accept), mean_acc, sigma, s_ans))
}

test = rq(y_init ~ X[,2], tau = tau)
test
nbatch = 10000 #100000
# init = coef(test)
out = amh(logf, init = coef(test) , scale = c(0.1/R, 0.1/R), nbatch = nbatch)
out[[2]]
plot(out[[3]], type = "l")
plot(out[[5]], type = "l")
out[[4]]
tail(out[[1]])
out[[1]] = out[[1]]*R
tail(out[[1]])
par(mfrow = c(2, 2))
plot(1:nbatch, out[[1]][,1], type = "l", ylab = "Value", xlab = "Runs", main = "Intercept")
plot(1:nbatch, out[[1]][,2], type = "l", ylab = "Value", xlab = "Runs", main = "Slope")
acf(out[[1]][,1], main = "Intercept")
acf(out[[1]][,2], main = "Slope")

print(effectiveSize(out[[1]][,1]))
print(effectiveSize(out[[1]][,2]))
print(effectiveSize(out[[1]][(nbatch*0.2):nbatch,1]))
print(effectiveSize(out[[1]][(nbatch*0.2):nbatch,2]))

print(effectiveSize(out[[1]][0:(nbatch*0.2),1]))
print(effectiveSize(out[[1]][0:(nbatch*0.2),2]))


tab = rbind(c(effectiveSize(out[[1]][,1]), effectiveSize(out[[1]][,2])),
            c(effectiveSize(out[[1]][(nbatch*0.2):nbatch,1]), effectiveSize(out[[1]][(nbatch*0.2):nbatch,2])))

colnames(tab) = c("Intercept", "Slope")
rownames(tab) = c("Without Burnins", "With Burnins")


tab

acf(out[[1]][(nbatch*0.2):nbatch,1])
acf(out[[1]][(nbatch*0.2):nbatch,2])

acf(out[[1]][0:(nbatch*0.5),1])
acf(out[[1]][0:(nbatch*0.5),2])


# gelman.diag(rbind(as.mcmc.list(as.mcmc(out1)), as.mcmc.list(as.mcmc(out[[1]][,1]))))
# gelman.diag(rbind(as.mcmc.list(as.mcmc(out2)), as.mcmc.list(as.mcmc(out[[1]][,2]))))
# 
# out1 = out[[1]][,1]
# out2 = out[[1]][,2]

out = KNG(init = coef(test), ep = eps, tau = tau, sumX = sumX, X = X, Y = Y)


#sigma = summary(test, se = "bootstrap", covariance = TRUE)$cov/R
library(coda)
par(mfrow = c(1, 1))
eps = 1
reps = 1
output = matrix(NA, ncol = 2, nrow = reps)
for (i in 1:reps){
  print(i)
  out = amh(logf, scale = c(1/R, 1/R))
  print(out[[2]])
  #logfval_preprocessed = out[[3]]
  out[[1]] = out[[1]]*R
  plot(1:100000, out[[1]][,1], type = "l", ylab = "Value", xlab = "Runs", main = "Intercept")
  plot(1:100000, out[[1]][,2], type = "l", ylab = "Value", xlab = "Runs", main = "Slope")
  output[i, ] = out[[1]][100000,]
  print(effectiveSize(out[[1]][,1]))
  print(effectiveSize(out[[1]][,2]))
  acf(out[[1]][,1])
  acf(out[[1]][,2])
  # print(heidel.diag(out))
  # print(geweke.diag(out))
  
}
output
hist(na.omit(output[,1]), breaks = 25, xlab = "Value", ylim = c(0, 50),
     main = paste("Histogram of Intercept over", reps, "Runs \n with AMH"))
hist(na.omit(output[,2]), breaks = 25, xlab = "Value", ylim = c(0, 50),
     main = paste("Histogram of Slope over", reps, "Runs \n with AMH"))

save(output, file = "amh_eps1_output_preprocessed.Rdata")
#par(mfrow = c(1, 2))

runs = 1000
lambda = 0.01
a0 = 4
b0 = 2
a1 = 3 
b1 = 2
b2 = 1
x1 = rexp(n, lambda)
x2 = a0 + b0*x1 + rexp(n, lambda)
x3 = a1 + b1*x1 + b2*x2 + rexp(n, lambda)
data = cbind(x1, x2)
cor(data)
