library(mvtnorm)
library(quantreg)
library(coda)
library(matrixcalc)
library(corpcor)
library(ggplot2)
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
logf = function(theta, eps){
  incl_sum = ((Y - X%*%theta) <= 0)
  sumX = apply(X, 2, sum)
  ans = -eps/(4*(1-tau)*max(abs(X))) * norm(-tau*sumX + t(incl_sum)%*%X, type = "m") #- 0.00001*theta%*%theta
  return(ans)
}

amh = function(logf, init = c(1, 1)/R, nbatch = 100000, scale, eps, nonneg = FALSE){
  #logfval = rep(NA, nrow = 100000)
  dim = length(init)
  batch = matrix(0, nrow = nbatch, ncol = dim)
  accept = matrix(0, nrow = nbatch)
  batch[1, ] = init
  s = 5
  sigma = diag(scale)/s # first assume it is indep
  ct = 0
  mean_acc = NA
  s_ans = NA
  scale_sigma = 1
  for (i in 2:nbatch){
    propose_val = rmvnorm(1, t(batch[i-1,]), sigma)
    # propose_val = rep(0, dim)
    # for (j in 1:dim){
    #   propose_val[j] = rnorm(1, batch[i-1, j], scale[j])
    # }
    u = runif(1, 0, 1)
    check = TRUE
    #logfval[i] = logf(c(propose_val)) - logf(batch[i-1, ])
    check = all(0<= X%*%t(propose_val)) #& all(X%*%t(propose_val) <= max(Y))
    if ((log(u) < logf(c(propose_val), eps) - logf(batch[i-1, ], eps)) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
      #print(ct)
    } else {
      batch[i, ] = batch[i-1, ]
    }
    if (ct == 10){
      sigma = cor(na.omit(batch)) + diag(dim)*0.01 #diag(c(runif(1, 0, 0.25), runif(1, 0, 0.25)))#/100#/2000 #2000 #100
      sigma[1, 1]= sigma[1, 1]
      sigma[1, 2] = sigma[1, 2]#/sqrt(3) #7
      sigma[2, 1] = sigma[2, 1]#/sqrt(3)
      sigma[2, 2]= sigma[2, 2]#/3
      #sigma = sigma*0.1
      temp_acc = mean(accept[1:i])
      if (temp_acc < 0.2){
        scale_sigma = scale_sigma/2#*0.1
      } else if (temp_acc > 0.25) {
        scale_sigma = scale_sigma*2
      }
      #print(is.positive.semi.definite(sigma))
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
      #print(sigma)
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    s_ans = cbind(s_ans, max(abs(y_init))*s)
  }
  mean_acc = mean_acc[-1]
  s_ans = s_ans[-1]
  return(list(batch, mean(accept), mean_acc, sigma, s_ans))
}

# test = rq(y_init ~ X[,2], tau = tau)
# test
# nbatch = 10000 #100000
# # init = coef(test)
# out = amh(logf, init = c(5, 1) , scale = c(0.1/R, 0.1/R), nbatch = nbatch, eps = eps)
# out[[2]]
# plot(out[[3]], type = "l")
# plot(out[[5]], type = "l")
# out[[4]]
# tail(out[[1]])
# out[[1]] = out[[1]]*R
# tail(out[[1]])
# plot(1:nbatch, out[[1]][,1], type = "l", ylab = "Value", xlab = "Runs", main = "Intercept")
# plot(1:nbatch, out[[1]][,2], type = "l", ylab = "Value", xlab = "Runs", main = "Slope")
# print(effectiveSize(out[[1]][,1]))
# print(effectiveSize(out[[1]][,2]))
# print(effectiveSize(out[[1]][(nbatch*0.3):nbatch,1]))
# print(effectiveSize(out[[1]][(nbatch*0.3):nbatch,2]))
# 
# print(effectiveSize(out[[1]][0:(nbatch*0.5),1]))
# print(effectiveSize(out[[1]][0:(nbatch*0.5),2]))
# 
# tab = rbind(c(effectiveSize(out[[1]][,1]), effectiveSize(out[[1]][,2])),
#             c(effectiveSize(out[[1]][(nbatch*0.2):nbatch,1]), effectiveSize(out[[1]][(nbatch*0.2):nbatch,2])))
# 
# colnames(tab) = c("Intercept", "Slope")
# rownames(tab) = c("Without Burnins", "With Burnins")
# tab
# acf(out[[1]][,1])
# acf(out[[1]][,2])
# acf(out[[1]][(nbatch*0.3):nbatch,1])
# acf(out[[1]][(nbatch*0.3):nbatch,2])
# acf(out[[1]][0:(nbatch*0.5),1])
# acf(out[[1]][0:(nbatch*0.5),2])


nbatch = 10000
reps = 50
beta1 = matrix(NA, ncol = 2, nrow = reps)
beta05 = matrix(NA, ncol = 2, nrow = reps)
beta01 = matrix(NA, ncol = 2, nrow = reps)
dist1 = matrix(NA, ncol = 2, nrow = reps)
dist05 = matrix(NA, ncol = 2, nrow = reps)
dist01 = matrix(NA, ncol = 2, nrow = reps)
invert_quantile = function(tau) log(1-tau)/(-lambda)
truth = c(4 + invert_quantile(tau), 3)
for (i in 1:reps){
  print(i)
  out1 = amh(logf, init = c(5, 1) , scale = c(0.1/R, 0.1/R), nbatch = nbatch, eps = 1)
  val1 = tail(out1[[1]], 1) 
  beta1[i, ] = val1
  dist1[i, ] = val1-truth
  
  out05 = amh(logf, init = c(5, 1) , scale = c(0.1/R, 0.1/R), nbatch = nbatch, eps = 0.5)
  val05 = tail(out05[[1]], 1) 
  beta05[i, ] = val05
  dist05[i, ] = val05-truth
  
  out01 = amh(logf, init = c(5, 1) , scale = c(0.1/R, 0.1/R), nbatch = nbatch, eps = 0.1)
  val01 = tail(out01[[1]], 1) 
  beta01[i, ] = val01
  dist01[i, ] = val01-truth
}

apply(dist01, 2, median)
apply(dist05, 2, median)
apply(dist1, 2, median)
l2_dist01 = sqrt(apply(dist01^2, 1, sum))
l2_dist05 = sqrt(apply(dist05^2, 1, sum))
l2_dist1 = sqrt(apply(dist1^2, 1, sum))
median(l2_dist01)
median(l2_dist05)
median(l2_dist1)
tab_plot = rbind(mean(l2_dist01), mean(l2_dist05), mean(l2_dist1))
tab_plot = cbind(tab_plot, rbind(sd(l2_dist01)/reps, sd(l2_dist05)/reps, sd(l2_dist1))/reps)
tab_plot = cbind(tab_plot, c(0.1, 0.5, 1))
colnames(tab_plot) = c("Mean", "SE", "Eps")
tab_plot = as.data.frame(tab_plot)
tab_plot$Eps = as.factor(tab_plot$Eps)
#tab_plot = as.numeric(tab_plot)
ggplot(tab_plot, aes(x = Eps, y = Mean, color = Eps)) + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.1) +
  geom_point()

