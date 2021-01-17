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
      if(U[r,i] < TestValue & check){
        out$batch[r,] = newVector
        out$accept = out$accept + 1/(nbatch*dim)
        oldLogA = newLogA
      }
    }
  }
  
  return(out)
}

getScale = function(logA, m, nbatch, start_scale = 0.1, lower_accept = 0.1, upper_accept = 0.2, nonneg){
  scale = start_scale
  prevBelow = 0
  prevAbove = Inf
  count = 0
  
  init = rep(0, m)
  out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, nonneg = nonneg)

  while((out$accept < lower_accept | out$accept > upper_accept) & (count <= 10)){
    if (out$accept < lower_accept) {
      prevAbove = scale 
    } else if (out$accept > .2) {
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
    message('Scale did not converge. If this continues, adjust start_scale or acceptance rate range.')
    return(0)
  }
  
  message("Scale", scale)
  return(scale)
}



KNG = function(ep, tau, sumX, X, Y, nbatch = 1000, scale = 0, start_scale = 0.1,
               lower_accept = 0.1, upper_accept = 0.2, nonneg = FALSE){
  m = ncol(X)
  init = rep(0, m)
  logA = function(beta){
    left = cbind(Y, X)%*% c(1, -beta)
    lessEq = (left <= 0)
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*1)+ (-1/2)*(beta%*%beta)
    return(ans)
  }
  while(scale == 0){
    scale = getScale(logA = logA, m = m, nbatch = nbatch, start_scale = start_scale, 
                     lower_accept = lower_accept, upper_accept = upper_accept, nonneg = nonneg)
  }
  
  out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, nonneg = nonneg)
  beta_kng = t(tail(out$batch, n = 1))
  
  return(list(beta_kng, scale))
}



#consider deleting m and change init to (initial/check)beta
#let check_data be NA or NULL by default
#change logD to logA
constrMetrop = function(logA, init, nbatch = 1000, scale, check_data, nonneg, 
                             method = c("koenker", "currentdata", "newdata"), type = c("upper", "lower")){
  method = match.arg(method)
  type = match.arg(type)
  
  dim = length(init)
  out = list(accept = 0, batch = matrix(rep(0, nbatch*dim), nrow = nbatch, ncol = dim))
  if (method == "koenker"){
    dim = 1
  }
  U = matrix(runif(nbatch*dim, min = 0, max = 1), nrow = nbatch, ncol = dim)
  Prop = matrix(rnorm(nbatch*dim, m = 0, s = scale), nrow = nbatch, ncol = dim)
  
  
  for (r in 1:nbatch) {
    if(r == 1){
      out$batch[1,] = init
      oldLogD = logA(init)
    } else {
      out$batch[r,] = out$batch[r-1,]
    }
    
    for (i in 1:dim){
      oldVector = out$batch[r,]
      newVector = oldVector
      newVector[i] = newVector[i] - Prop[r,i]
      
      
      newLogD = logA(newVector)
      if (method == "koenker") {
        check = ifelse(type == "lower", newVector[i] < init[i], newVector[i] > init[i])
      } else if (method == "currentdata" | method == "newdata") {
        check = ifelse(type == "lower",
                       all(check_data %*% as.matrix(newVector) < check_data %*% as.matrix(init)),
                       all(check_data %*% as.matrix(newVector) > check_data %*% as.matrix(init)))
      }
      if (nonneg){
        check = check & newVector[i] >= 0
      }
      
      TestValue = exp(newLogD - oldLogD)
      if (U[r,i] < TestValue & check) {
        out$batch[r,] = newVector
        out$accept = out$accept + 1/(nbatch*dim)
        oldLogD = newLogD
      }
    }
  }
  
  return(out)
}

#should scale be recalculated for each time?
#what is the use of m in this function
#allows nbatch to change
constrGetScale = function(logA, init, nbatch, start_scale, lower_accept = 0.2, 
                               upper_accept = 0.25, check_data, nonneg, 
                               method = c("koenker", "currentdata", "newdata"), type = c("upper", "lower")){
  type = match.arg(type)
  method = match.arg(method)
  
  scale = start_scale
  prevBelow = 0
  prevAbove = Inf
  count = 0
  
  out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale,
                          check_data = check_data, nonneg = nonneg, method = method, type = type)
  start_accept = out$accept
  print(out$accept)
  while((out$accept < lower_accept | out$accept > upper_accept) & (count <= 10)){
    if (out$accept < lower_accept) {
      prevAbove = scale
    } else if(out$accept > upper_accept) {
      prevBelow = scale
    }

    if (prevAbove < Inf) {
      scale = (1/2)*(prevBelow+prevAbove)
    } else {
      scale = 2*scale
    }
    
    out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale, 
                            check_data = check_data, nonneg = nonneg, method = method, type = type)
    count = count + 1
  }
  
  if (count == 11) {
    message('Scale did not converge. If this continues, adjust start_scale or acceptance rate range.', out$accept)
    return(list(0, start_accept))
  }
  
  message("Scale", scale)
  print(out$accept)
  return(list(scale, start_accept))
}

#what does blen do? should it be deleted?
#manipulate check_data here
constrKNG = function(ep, tau, sumX, X, Y, init, nbatch = 1000, scale = 0, start_scale = 0.5/max(Y),
                          lower_accept = 0.2, upper_accept = 0.25, check_data = NULL, nonneg = FALSE,
                          method = c("koenker", "currentdata", "newdata"), type = c("upper", "lower")){
  if (is.null(check_data)){
    if (method == "currentdata"){
      check_data = X
    } else if (method == "newdata"){
      message("No data input. Method changed to Koenker.")
    }
  } else {
    if (method == "newdata") {
      check_data = as.matrix(check_data)
      check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
    }
  }
  
  logA = function(beta) {
    left = cbind(Y, X) %*% c(1,-beta)
    lessEq = (left <= 0)
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*1) + (-1/2)*(beta%*%beta)
    return(ans)
  }
  
  count = 0
  while (scale == 0) {
    ans = constrGetScale(logA = logA, init = init, nbatch = nbatch, start_scale = start_scale, 
                                lower_accept = lower_accept, upper_accept = upper_accept,
                                check_data = check_data, nonneg = nonneg, method = method, type = type)
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
  
  out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale, check_data = check_data,
                          nonneg = nonneg, method = method, type = type)
  beta_kng = t(tail(out$batch, n = 1))
  return(list(beta_kng , scale, out$accept))
}

#instead of having lower_scale and upper_scale, scale should be input as a vector, with
#the same order as tau. If the same scale is used for all the quantiles, then it should be a
#vector of the same value
stepwiseKNG = function(data, total_eps, median_eps = 0.7, tau, scale = rep(0, length(tau)),
                      nbatch = 1000, start_scale = NULL, median_start_scale = 0.1,
                      lower_accept = 0.2, upper_accept = 0.25, check_data = NULL, d = 1e-04,
                      method = c("koenker", "currentdata", "newdata"), nonneg = FALSE){
  method = match.arg(method)
  data = as.matrix(data)
  ep = total_eps*(1-median_eps)/(length(tau)-1)
  scale_kng = 0
  i = ncol(data)
  Y = data[,i]
  R = max(abs(Y))
  Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  sumX = apply(X = X, 2, FUN = sum)
  m = ncol(X) - 1
  names(scale) = tau
  median_start_scale = 0.1/R
  if (m == 0 & method == "newdata"){
    method = "currentdata"  
    #to avoid not having any data to input in the case of private quantile
    #needs to be better optimized for the package
    #one option is asking user to input checkdata as a design matrix
  }
  
  out = KNG(ep = total_eps*median_eps, tau = 0.5, sumX = sumX, X = X, Y = Y, 
            nbatch = nbatch, scale = scale["0.5"], start_scale = median_start_scale, 
            upper_accept = 0.2, lower_accept = 0.1, nonneg = nonneg)
  median_beta_kng = out[[1]]*R
  scale_output = out[[2]]
  ans = median_beta_kng
  
  tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
  tau_upper = sort(tau[tau > 0.5])
  if (is.null(start_scale)) {
    start_scale = 0.1/R
  }
  
  if (length(tau_lower) > 0){
    check_beta = median_beta_kng
    for (i in 1:length(tau_lower)){
      new_tau = tau_lower[i]
      print(new_tau)
      out = constrKNG(ep = ep, tau = new_tau, sumX = sumX, X = X, Y = Y, init = c(check_beta)/R,
                      scale = scale[names(scale) == new_tau], start_scale = start_scale, 
                      check_data = check_data, method = method, type = "lower", nonneg = nonneg)
      proposed_beta = out[[1]]*R
      message("Acceptance rate of ", new_tau, ": ", out[[3]])
      
      check_beta = proposed_beta
      ans = cbind(proposed_beta, ans)
      scale_output = cbind(out[[2]], scale_output)

    }
  }
  
  if (length(tau_upper) > 0){
    check_beta = median_beta_kng
    for (i in 1:length(tau_upper)){
      new_tau = tau_upper[i]
      print(new_tau)
     
      out = constrKNG(ep = ep, tau = new_tau, sumX = sumX, X = X, Y = Y, init = c(check_beta)/R,
                      scale = scale[names(scale) == new_tau], start_scale = start_scale, 
                      check_data = check_data, method = method, type = "upper", nonneg = nonneg)
      proposed_beta = out[[1]]*R
      message("Acceptance rate of ", new_tau, ": ", out[[3]])
        
      check_beta = proposed_beta
      ans = cbind(ans, proposed_beta)
      scale_output = cbind(scale_output, out[[2]])
    }
  }
  
  tau = sort(tau)
  colnames(ans) = tau
  colnames(scale_output) = tau
  rownames(scale_output) = "Scale"
  
  return(list(ans, scale_output))
}


constrMetropSandwich = function(logA, init, nbatch = 1000, scale, lowerbeta, upperbeta, check_data, 
                                method = c("koenker", "currentdata", "newdata"), nonneg){
  method = match.arg(method)
  
  dim = length(init)
  out = list(accept = 0,batch = matrix(rep(0,nbatch*dim),nrow=nbatch,ncol=dim))
  if (method == "koenker") {
    dim = 1
  }
  U = matrix(runif(nbatch*dim, min = 0, max = 1), nrow = nbatch, ncol = dim)
  Prop = matrix(rnorm(nbatch*dim, m = 0, s = scale), nrow = nbatch, ncol = dim)
  
  for(r in 1:nbatch){
    if (r == 1){
      out$batch[1,] = init
      oldlogA = logA(init)
    } else {
      out$batch[r,] = out$batch[r-1,]
    }
    
    for (i in 1:dim) {
      oldVector = out$batch[r,]
      newVector = oldVector
      newVector[i] = newVector[i] - Prop[r,i]
      
      
      newlogA = logA(newVector)
      if (method == "koenker"){
        check = (newVector[i] < upperbeta[i]) & (newVector[i] > lowerbeta[i])  
      } else {
        check = all(check_data %*% as.matrix(lowerbeta) < check_data %*% as.matrix(newVector)) &
          all(check_data %*% as.matrix(newVector) < check_data %*% as.matrix(upperbeta))
      }
      
      if (nonneg){
        check = check & newVector[i] >= 0
      }
      
      TestValue = exp((newlogA - oldlogA))
      if(U[r,i] < TestValue & check){
        out$batch[r,] = newVector
        out$accept = out$accept + 1/(nbatch*dim)
        oldlogA = newlogA
      }
    }
  }
  return(out)
}

constrGetScaleSandwich = function(logA, init, nbatch, start_scale, lowerbeta, upperbeta, lower_accept = 0.2, 
                          upper_accept = 0.25, check_data, nonneg, method = c("koenker", "currentdata", "newdata")){
  method = match.arg(method)
  
  scale = start_scale
  prevBelow = 0
  prevAbove = Inf
  count = 0
  
  out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale, lowerbeta = lowerbeta,
                             upperbeta = upperbeta, check_data = check_data, nonneg = nonneg, method = method)
  print(out$accept)
  start_accept = out$accept
  while((out$accept < lower_accept | out$accept > upper_accept) & (count <= 10)){
    if (out$accept < lower_accept) {
      prevAbove = scale
    } else if(out$accept > upper_accept) {
      prevBelow = scale
    }
    
    if (prevAbove < Inf) {
      scale = (1/2)*(prevBelow+prevAbove)
    } else {
      scale = 2*scale
    }
    
    out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale, lowerbeta = lowerbeta,
                               upperbeta = upperbeta, check_data = check_data, nonneg = nonneg, method = method)
    print(out$accept)
    count = count + 1
  }
  
  if (count == 11) {
    message('Scale did not converge. If this continues, adjust start_scale or acceptance rate range.')
    print(out$accept)
    return(list(0, start_accept))
  }
  
  return(list(scale, start_accept))
}

constrKNGSandwich = function(ep, tau, sumX, X, Y, init, nbatch = 1000, scale = 0, start_scale = 0.5/max(Y),
                             lowerbeta, upperbeta, lower_accept = 0.2, upper_accept = 0.25, check_data = NULL, 
                             nonneg = FALSE, method = c("koenker", "currentdata", "newdata")){
  if (is.null(check_data)){
    if (method == "currentdata"){
      check_data = X
    } else if (method == "newdata"){
      message("No data input. Method changed to Koenker.")
    }
  } else {
    if (method == "newdata") {
      check_data = as.matrix(check_data)
      check_data = as.matrix(cbind(rep(1, nrow(check_data)), check_data))
    }
  }
  
  logA = function(beta) {
    left = cbind(Y, X) %*% c(1,-beta)
    lessEq = (left <= 0)
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*1) + (-1/2)*(beta%*%beta)
    return(ans)
  }
  
  count = 0
  while (scale == 0) {
    ans = constrGetScaleSandwich(logA = logA, init = init, nbatch = nbatch, start_scale = start_scale, 
                                   lowerbeta = lowerbeta, upperbeta = upperbeta, lower_accept = lower_accept, 
                                   upper_accept = upper_accept, check_data = check_data, nonneg = nonneg, method = method)
    print(ans)
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
  
  out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale, lowerbeta = lowerbeta,
                             upperbeta = upperbeta, check_data = check_data, nonneg = nonneg, method = method)
  beta_kng = t(tail(out$batch, n = 1))
  return(list(beta_kng, scale, out$accept))
}

#differnt lower and upper scale helps find scale easier (due to the long tail)
#scale needs to be a vector of the same order as tau
#quantile vector does not need to be in an order
sandwichKNG = function(data, total_eps, median_eps = 0.7, main_tau_eps = 0.5,
                       tau, main_tau, scale = NULL, nbatch = 1000, 
                       median_start_scale = 0.1, main_tau_start_scale = NULL, sandwich_start_scale = NULL,
                       lower_accept = 0.2, upper_accept = 0.25, check_data = NULL, d = 1e-04,
                       method = c("koenker", "currentdata", "newdata"), nonneg = FALSE){
  main_tau_fac = as.factor(main_tau)
  sandwich_tau = tau[!tau %in% main_tau_fac]
  
  if(is.null(scale)){
    scale = rep(0, length(tau))
  }
  
  data = as.matrix(data)
  i = ncol(data)
  Y = as.matrix(data[,i])
  R = max(abs(Y))
  Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  sumX = apply(X = X, 2, FUN = sum)
  m = ncol(X) - 1
  if (m == 0 & method == "newdata"){
    method = "currentdata"  
    #to avoid not having any data to input in the case of private quantile
    #needs to be better optimized for the package
    #one option is asking user to input checkdata as a design matrix
  }
  names(scale) = tau
  main_tau_scale = scale[names(scale) %in% main_tau]
  sandwich_scale = scale[!names(scale) %in% main_tau]
  main_tau_order = tau[tau %in% main_tau]

  
  out = stepwiseKNG(data = data, total_eps = total_eps*main_tau_eps, median_eps = median_eps, 
                  tau = main_tau_order, scale = main_tau_scale, nbatch = nbatch, 
                  start_scale = main_tau_start_scale, median_start_scale = median_start_scale, 
                  lower_accept = lower_accept, upper_accept = upper_accept, check_data = check_data, 
                  d = d, method = method, nonneg = nonneg)
  b = out[[1]]
  scale_output = out[[2]]
  main_tau = sort(main_tau)
  names(scale_output) = main_tau
  colnames(b) = main_tau

  eps_sandwich = (1 - main_tau_eps) * total_eps / length(sandwich_tau)
  if (is.null(sandwich_start_scale)) {
    sandwich_start_scale = 0.15/R
  }
  
  update_tau = main_tau
  sandwich_tau_lower = sort(sandwich_tau[sandwich_tau < 0.5], decreasing = TRUE)
  sandwich_tau_upper = sort(sandwich_tau[sandwich_tau > 0.5])
  if (length(sandwich_tau_lower) > 0){
    sandwich_scale_lower = rep(NA, length(sandwich_tau_lower))
    names(sandwich_scale_lower) = sandwich_tau_lower
    for (i in 1:length(sandwich_tau_lower)){
      curr_tau = sandwich_tau_lower[i]
      print(curr_tau)
      update_tau = sort(c(update_tau, curr_tau))
      idx = which(update_tau == curr_tau)
      lowertau = update_tau[idx-1]
      uppertau = update_tau[idx+1]
      
      lowerbeta = b[, which(colnames(b) == lowertau)]
      upperbeta = b[, which(colnames(b) == uppertau)]
      #main_tau_eps
      out = constrKNGSandwich(ep = eps_sandwich, tau = curr_tau, sumX = sumX, X = X, Y = Y, 
                                     init = c(upperbeta)/R, nbatch = nbatch, 
                                     scale = sandwich_scale[names(sandwich_scale) == curr_tau], 
                                     start_scale = sandwich_start_scale, lowerbeta = c(lowerbeta)/R, 
                                     upperbeta = c(upperbeta)/R, lower_accept = lower_accept,
                                     upper_accept = upper_accept, check_data = check_data, 
                                     nonneg = nonneg, method = method)
      b_sandwich = out[[1]]*R
      sandwich_scale_lower[i] = out[[2]]
      b = cbind(b, b_sandwich)
      colnames(b)[ncol(b)] = curr_tau
      if(m == 0){
        b = t(as.matrix(b[, match(update_tau, colnames(b))]))
      } else {
        b = as.matrix(b[, match(update_tau, colnames(b))])
      }
    }
  }
  
  if (length(sandwich_tau_upper) > 0){
    sandwich_scale_upper = rep(NA, length(sandwich_tau_upper))
    names(sandwich_scale_upper) = sandwich_tau_upper
    for (i in 1:length(sandwich_tau_upper)){
      curr_tau = sandwich_tau_upper[i]
      print(curr_tau)
      update_tau = sort(c(update_tau, curr_tau))
      idx = which(update_tau == curr_tau)
      lowertau = update_tau[idx-1]
      uppertau = update_tau[idx+1]
      
      lowerbeta = b[, which(colnames(b) == lowertau)]
      upperbeta = b[, which(colnames(b) == uppertau)]
      #good sandwich startscale 1/R
      out = constrKNGSandwich(ep = eps_sandwich, tau = curr_tau, sumX = sumX, X = X, Y = Y, 
                                     init = c(lowerbeta)/R, nbatch = nbatch, 
                                     scale = sandwich_scale[names(sandwich_scale) == curr_tau], 
                                     start_scale = 1/R, lowerbeta = c(lowerbeta)/R, 
                                     upperbeta = c(upperbeta)/R, lower_accept = lower_accept,
                                     upper_accept = upper_accept, check_data = check_data, 
                                     nonneg = nonneg, method = method)
      b_sandwich = out[[1]]*R
      sandwich_scale_upper[i] = out[[2]]
      b = cbind(b, b_sandwich)
      colnames(b)[ncol(b)] = curr_tau
      if(m == 0){
        b = t(as.matrix(b[, match(update_tau, colnames(b))]))
      } else {
        b = as.matrix(b[, match(update_tau, colnames(b))])
      }
    }
  }
  
  scale_output = c(scale_output, sandwich_scale_upper, sandwich_scale_lower)
  scale_output = scale_output[order(names(scale_output))]
  
  return(list(b, scale_output))
  
}


