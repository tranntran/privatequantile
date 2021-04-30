metrop = function(logA, init, nbatch = 1000, scale, nonneg = FALSE){
  dim = length(init)
  
  out = list(accept = 0, batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim))
  U = matrix(runif(nbatch*dim, min = 0, max = 1), nrow = nbatch, ncol = dim)
  Prop = matrix(rnorm(nbatch*dim, m = 0, s = scale), nrow = nbatch, ncol = dim)
  
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
      if (nonneg) {
        check = check & (newVector[i] > 0) 
      }
      
      TestValue = exp((newLogA - oldLogA))
      if(log(U[r,i]) < newLogA - oldLogA & check){
        out$batch[r,] = newVector
        out$accept = out$accept + 1/(nbatch*dim)
        oldLogA = newLogA
      }
    }
  }
  
  return(out)
}

#par(mfrow = c(1, 2))
plot(1:100000, out$batch[,1]*R, type = "l", xlab = "runs", ylab = "Intercept")
plot(1:100000, out$batch[,2]*R, type = "l", xlab = "runs", ylab = "Slope")
acf(out$batch[,1])
acf(out$batch[,2])

getScale = function(logA, m, nbatch, start_scale = 0.1, lower_accept = 0.2, upper_accept = 0.25, nonneg){
  scale = start_scale
  prevBelow = 0
  prevAbove = Inf
  count = 0
  
  init = rep(0, m)
  out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, nonneg = nonneg)
  start_accept = out$accept
  while((out$accept < lower_accept | out$accept > upper_accept) & (count <= 10)){
    if (out$accept < lower_accept) {
      prevAbove = scale 
    } else if (out$accept > upper_accept) {
      prevBelow = scale 
    }
    
    if (prevAbove < Inf) {
      scale = (1/2)*(prevBelow + prevAbove) 
    } else {
      scale = 2*scale 
    }
    
    out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, nonneg = nonneg)
    count = count + 1
  }
  
  if (count == 11) {
    return(list(0, start_accept))
  }
  
  return(list(scale, start_accept))
}



KNG = function(ep, tau, sumX, X, Y, nbatch = 1000, scale = 0, start_scale = 0.1,
               lower_accept = 0.1, upper_accept = 0.2, nonneg = FALSE){
  m = ncol(X)
  init = rep(0, m)
  logA = function(beta){
    left = cbind(Y, X)%*% c(1, -beta)
    lessEq = (left <= 0)
    ans = -(ep/2) * norm(-tau*sumX + t(X)%*%lessEq, type = "M") / ((1-tau)*2*110) - 1/2*beta%*%beta
    #print(-tau*sumX + t(X)%*%lessEq)
    return(ans)
  }
  
  count = 0
  while(scale == 0){
    ans = getScale(logA = logA, m = m, nbatch = nbatch, start_scale = start_scale, 
                   lower_accept = lower_accept, upper_accept = upper_accept, nonneg = nonneg)
    scale = ans[[1]]
    start_accept = ans[[2]]
    count = count + 1
    if (count == 2){
      if (start_accept < lower_accept){
        message("Unable to get scale to converge. Automatically decrease starting scale.")
        start_scale = start_scale*0.5
      } else if (start_accept > upper_accept){
        message("Unable to get scale to converge. Automatically increase starting scale.")
        start_scale = start_scale*5
      }
      count = 0
    }
  }
  
  
  out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, nonneg = nonneg)
  beta_kng = t(tail(out$batch, n = 1))
  
  return(list(beta_kng, scale, out$accept))
}


originalKNG = function(data, total_eps, tau, nbatch = 1000, scale = rep(0, length(tau)), 
                       start_scale = 0.001, lower_accept = 0.2, upper_accept = 0.25, 
                       nonneg = FALSE){
  ep = total_eps/length(tau)
  scale_kng = 0
  i = ncol(data)
  Y = data[,i]
  R = max(abs(Y))
  Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  sumX = apply(X = X, 2, FUN = sum)
  #m = ncol(X) - 1
  
  accept_rate = rep(NA, length(tau))
  scale_output = rep(NA, length(tau))
  ans = NA
  for (i in 1:length(tau)){
    print(tau[i])
    curr_scale = scale[i]
    temp = KNG(ep = ep, tau = tau[i], sumX = sumX, X = X, Y = Y, 
               nbatch = nbatch, scale = curr_scale, start_scale = start_scale, 
               upper_accept = upper_accept, lower_accept = lower_accept, 
               nonneg = nonneg)
    ans = cbind(ans, temp[[1]]*R)
    scale_output[i] = temp[[2]]
    accept_rate[i] = temp[[3]]
  }
  ans = ans[, -1]
  return(list(ans, scale_output, accept_rate))
}

ep = 0.25
n = 5000
runs = 1000
lambda = 0.1
a0 = 4
a1 = 3 
b0 = 10
b1 = 2
b2 = 1

x1 = rexp(n, lambda)
x2 = a0 + b0*x1 + rexp(n, lambda)
x3 = a1 + b1*x1 + b2*x2 + rexp(n, lambda)
data = cbind(x1, x2)
reps = 2500
true_beta = matrix(b0, nrow = 2, ncol = length(0.5))
invert_quantile = function (x) log(1-x)/(-lambda)
true_beta[1,] = sapply(0.5, invert_quantile) + a0
true_beta = apply(as.matrix(true_beta), 1, rep, reps)

originalKNG_0.1 = matrix(NA, ncol = 2, nrow = reps)
originalKNG_1 = matrix(NA, ncol = 2, nrow = reps)
originalKNG_0.01 = matrix(NA, ncol = 2, nrow = reps)
for (i in 1:reps){
  print(i)
  originalKNG_0.1[i, c(1:2)] = originalKNG(data = data, total_eps = 0.1, tau = 0.5, 
                                           start_scale = 0.01, scale = 1/max(x2))[[1]]
  
  originalKNG_1[i, c(1:2)] = originalKNG(data = data, total_eps = 1, tau = 0.5, 
                                         start_scale = 0.01, scale = 1/max(x2))[[1]]
  originalKNG_0.01[i, c(1:2)] = originalKNG(data = data, total_eps = 0.01, tau = 0.5, 
                                         start_scale = 0.01, scale = 1/max(x2))[[1]]
  
}
dist_0.1 = sqrt(rowSums((originalKNG_0.1 - true_beta)^2))
dist_1 = sqrt(rowSums((originalKNG_1 - true_beta)^2))
dist_0.01 = sqrt(rowSums((originalKNG_0.01 - true_beta)^2))
quantile(dist_0.1, na.rm = TRUE)
quantile(dist_1, na.rm = TRUE)
quantile(dist_0.01, na.rm = TRUE)
mean(dist_1, na.rm = TRUE)
sd(dist_0.01, na.rm = TRUE)
