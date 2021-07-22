rm(list = ls())
library(mvtnorm)
metrop = function(logA, init, nbatch = 10000, scale = 1e-4, nonneg = FALSE, X = X,
                  Y = Y, lower_accept = 0.2, upper_accept = 0.5,
                  update_after = 10, adjust_scale_by = 2){
  dim = length(init)
  
  batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim)
  u = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
  accept = matrix(0, nrow = nbatch)
  
  ct = 0
  mean_acc = NA
  scale_ans = NA
  sigma = diag(rep(1, dim))*scale
  
  batch[1, ] = init
  for (i in 2:nbatch){
    propose_val = mvtnorm::rmvnorm(1, t(batch[i-1,]), sigma)
    check = TRUE
    if (nonneg){
      check = check & all(X %*% t(propose_val) > min(Y))
      check = check & all(X %*% t(propose_val) < max(Y))
    }
    
    if ((log(u[i]) < logA(c(propose_val)) - logA(batch[i-1, ])) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
    } else {
      batch[i, ] = batch[i-1, ]
    }
    
    if (ct == update_after){
      sigma = cor(na.omit(batch)) + diag(dim)*0.00002 
      temp_acc = mean(accept[1:i])
      if (temp_acc < lower_accept){
        scale = scale/adjust_scale_by
      } else if (temp_acc > upper_accept) {
        #print(temp_acc)
        scale = scale*adjust_scale_by
      }
      sigma = sigma*scale
      #print(sigma)
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    scale_ans = cbind(scale_ans, scale)
  }
  mean_acc = mean_acc[-1]
  scale_ans = scale_ans[-1]
  
  out = list(batch = batch, accept_rate = mean(accept), scale = scale_ans, sigma = sigma)
  return(out)
}


KNG = function(init, ep, tau, sumX, X, Y, nbatch = 10000, scale = 1e-4, nonneg = FALSE,
               lower_accept = 0.2, upper_accept = 0.5, update_after = 10, 
               adjust_scale_by = 2){
  logA = function(beta){
    left = cbind(Y, X)%*% c(1, -beta)
    lessEq = (left <= 0)
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*1*max(X))
    #+ (-1/2)*(beta%*%beta)
    return(ans)
  }
  
  out = metrop(logA = logA, init = init, nbatch = nbatch, scale = scale, nonneg = nonneg,
               X = X, Y = Y, lower_accept = lower_accept, upper_accept = upper_accept,
               update_after = update_after, adjust_scale_by = adjust_scale_by)
  
  return(out)
}



#consider deleting m and change init to (initial/check)beta
#let check_data be NA or NULL by default
#change logD to logA
constrMetrop = function(logA, init, nbatch = 10000, scale = 1e-4, check_beta, check_data, 
                        Y = Y, nonneg, method = c("koenker", "currentdata", "newdata"), 
                        type = c("upper", "lower"), lower_accept = 0.2, upper_accept = 0.5,
                        update_after = 10, adjust_scale_by = 2){
  method = match.arg(method)
  type = match.arg(type)
  
  dim = length(init)
  
  batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim)
  u = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
  accept = matrix(0, nrow = nbatch)
  
  ct = 0
  mean_acc = NA
  scale_ans = NA
  sigma = diag(rep(1, dim))*scale
  
  batch[1, ] = init
  for (i in 2:nbatch){
    check = TRUE
    
    if (method == "koenker") {
      propose_val = batch[i-1, ]
      propose_val[1] = propose_val[1] + rnorm(1, 0, sigma)
      check = check & ifelse(type == "lower", propose_val[1] < check_beta[1], propose_val[1] > check_beta[1])
      # if (nonneg){
      #   check = check & (propose_val[1] >= 0)
      # }
      # add proposal here
      # add non neg here
      
    } else if (method == "currentdata" | method == "newdata") {
      propose_val = t(mvtnorm::rmvnorm(1, t(batch[i-1,]), sigma))
      check = check & ifelse(type == "lower",
                             all(check_data %*% propose_val < check_data %*% t(check_beta)),
                             all(check_data %*% propose_val > check_data %*% t(check_beta)))
      # if (nonneg){
      #   check = check & all(check_data %*% as.matrix(propose_val))
      # }
    }
    # print(dim(check_data))
    # print(as.matrix(propose_val))
    # print(t(t(propose_val)))
    check = check & max(as.matrix(check_data) %*% as.matrix(propose_val))<= max(Y)
    check = check & min(as.matrix(check_data) %*% as.matrix(propose_val)) >= min(Y)
    
    if ((log(u[i]) < logA(c(propose_val)) - logA(batch[i-1, ])) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
    } else {
      batch[i, ] = batch[i-1, ]
    }
    
    if (ct == update_after){
      # add case here for koenker use variance instead
      if (method == "koenker") {
        sigma = var(batch[,1])
      } else {
        sigma = cor(na.omit(batch)) + diag(dim)*0.00002 
        
      }
      temp_acc = mean(accept[1:i])
      if (temp_acc < lower_accept){
        scale = scale/adjust_scale_by
      } else if (temp_acc > upper_accept) {
        #print(temp_acc)
        scale = scale*adjust_scale_by
      }
      sigma = sigma*scale
      print(sigma)
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    scale_ans = cbind(scale_ans, scale)
  }
  mean_acc = mean_acc[-1]
  scale_ans = scale_ans[-1]
  
  out = list(batch = batch, accept_rate = mean(accept), scale = scale_ans, sigma = sigma)
  
  return(out)
}



#what does blen do? should it be deleted?
#manipulate check_data here
constrKNG = function(init, ep, tau, sumX, X, Y, nbatch = 10000, scale = 1e-4,
                     check_beta, check_data = NULL, nonneg = FALSE, 
                     method = c("koenker", "currentdata", "newdata"), 
                     type = c("upper", "lower"), lower_accept = 0.2, upper_accept = 0.5, 
                     update_after = 10, adjust_scale_by = 2){
  if (is.null(check_data)){
    if (method == "newdata"){
      message("No data input. Method changed to Koenker.")
    } else {
      check_data = X
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
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*1*max(X)) 
    #+ (-1/2)*(beta%*%beta)
    return(ans)
  }
  
  out = constrMetrop(logA = logA, init = init, nbatch = nbatch, scale = scale,
                     check_beta = check_beta, check_data = check_data, Y = Y,
                     nonneg = nonneg, method = method, type = type, 
                     lower_accept = lower_accept, upper_accept = upper_accept,
                     update_after = update_after, adjust_scale_by = adjust_scale_by)
  return(out)
}

#instead of having lower_scale and upper_scale, scale should be input as a vector, with
#the same order as tau. If the same scale is used for all the quantiles, then it should be a
#vector of the same value
stepwiseKNG = function(data, total_eps, median_eps = NULL, tau, scale = 1e-4,
                       nbatch = 10000, method = c("koenker", "currentdata", "newdata"), 
                       nonneg = FALSE, check_data = NULL, lower_accept = 0.2, 
                       upper_accept = 0.5, update_after = 10, adjust_scale_by = 2,
                       formula = NULL){
  method = match.arg(method)
  print(method)
  #data = as.matrix(data)
  if(is.null(median_eps)){
    median_eps = ifelse(method == "koenker", 0.7, 0.4)
  }
  ep = total_eps*(1-median_eps)/(length(tau)-1)
  i = ncol(data)
  Y = data[,i]
  #R = max(abs(Y))
  #Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  sumX = apply(X = X, 2, FUN = sum)
  m = ncol(X) - 1
  #names(scale) = tau
  
  if (is.null(formula)) {
    vars = colnames(data)
    formula = paste(vars[i], " ~ .")
  }
  
  
  if (m == 0 & method == "newdata"){
    method = "currentdata"  
    #to avoid not having any data to input in the case of private quantile
    #needs to be better optimized for the package
    #one option is asking user to input checkdata as a design matrix
  }
  
  #nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
  out = KNG(init = rep(1, m+1), ep = total_eps*median_eps, tau = 0.5, sumX = sumX, X = X, 
            Y = Y, nbatch = nbatch, scale = scale, nonneg = nonneg, upper_accept = upper_accept, 
            lower_accept = lower_accept, update_after = update_after, adjust_scale_by = adjust_scale_by)
  median_beta_kng = tail(out[[1]], 1)#*R
  accept_rate = out[[2]]
  scale_output = tail(out[[3]], 1)
  ans = t(median_beta_kng)
  
  tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
  tau_upper = sort(tau[tau > 0.5])
  
  if (length(tau_lower) > 0){
    check_beta = median_beta_kng
    for (i in 1:length(tau_lower)){
      new_tau = tau_lower[i]
      print(new_tau)
      nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = new_tau)
      out = constrKNG(init = check_beta, ep = ep, tau = new_tau, sumX = sumX, X = X, Y = Y, 
                      nbatch = nbatch, scale = scale, check_beta = check_beta, 
                      check_data = check_data, nonneg = nonneg, method = method, 
                      type = "lower", lower_accept = lower_accept, upper_accept = upper_accept, 
                      update_after = update_after, adjust_scale_by = adjust_scale_by)
      proposed_beta = tail(out[[1]], 1)#*R
      message("Acceptance rate of ", new_tau, ": ", out[[2]])
      
      check_beta = proposed_beta
      ans = cbind(t(proposed_beta), ans)
      scale_output = cbind(tail(out[[3]], 1), scale_output)
      accept_rate = cbind(out[[2]], accept_rate)
      
    }
  }
  
  if (length(tau_upper) > 0){
    check_beta = median_beta_kng
    for (i in 1:length(tau_upper)){
      new_tau = tau_upper[i]
      print(new_tau)
      nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = new_tau)
      out = constrKNG(init = check_beta, ep = ep, tau = new_tau, sumX = sumX, X = X, Y = Y, 
                      nbatch = nbatch, scale = scale, check_beta = check_beta, 
                      check_data = check_data, nonneg = nonneg, method = method, 
                      type = "upper", lower_accept = lower_accept, upper_accept = upper_accept, 
                      update_after = update_after, adjust_scale_by = adjust_scale_by)
      proposed_beta = tail(out[[1]], 1)#*R
      message("Acceptance rate of ", new_tau, ": ", out[[2]])
      
      check_beta = proposed_beta
      ans = cbind(ans, t(proposed_beta))
      scale_output = cbind(scale_output, tail(out[[3]], 1))
      accept_rate = cbind(accept_rate, out[[2]])
    }
  }
  
  tau = sort(tau)
  colnames(ans) = tau
  colnames(scale_output) = tau
  rownames(scale_output) = "Scale"
  colnames(accept_rate) = tau
  
  return(list(ans, scale_output, accept_rate))
}


constrMetropSandwich = function(logA, init, nbatch = 10000, scale = 1e-4, lowerbeta, upperbeta, 
                                check_data, Y = Y, nonneg, method = c("koenker", "currentdata", "newdata"),
                                lower_accept = 0.2, upper_accept = 0.5, update_after = 10, 
                                adjust_scale_by = 2){
  
  method = match.arg(method)
  
  dim = length(init)
  
  batch = matrix(rep(0, nbatch*dim),nrow = nbatch, ncol = dim)
  u = matrix(runif(nbatch, min = 0, max = 1), nrow = nbatch)
  accept = matrix(0, nrow = nbatch)
  
  ct = 0
  mean_acc = NA
  scale_ans = NA
  sigma = diag(rep(1, dim))*scale
  
  batch[1, ] = init
  for (i in 2:nbatch){
    check = TRUE
    
    if (method == "koenker") {
      propose_val = batch[i-1, ]
      propose_val[1] = propose_val[1] + rnorm(1, 0, sigma)
      check = check & (propose_val[1] < upperbeta[1]) & (propose_val[1] > lowerbeta[1])
      
    } else if (method == "currentdata" | method == "newdata") {
      propose_val = t(mvtnorm::rmvnorm(1, t(batch[i-1,]), sigma))
      check = check & 
        all(check_data %*% propose_val < check_data %*% upperbeta) &
        all(check_data %*% propose_val > check_data %*% lowerbeta)
    }
    
    check = check & (max(as.matrix(check_data) %*% as.matrix(propose_val)) <= max(Y))
    check = check & (min(as.matrix(check_data) %*% as.matrix(propose_val)) >= min(Y))
    
    if ((log(u[i]) < logA(c(propose_val)) - logA(batch[i-1, ])) & check){
      batch[i, ] = propose_val
      accept[i,] = 1
      ct = ct + 1
    } else {
      batch[i, ] = batch[i-1, ]
    }
    
    if (ct == update_after){
      # add case here for koenker use variance instead
      if (method == "koenker") {
        sigma = var(batch[,1])
      } else {
        sigma = cor(na.omit(batch)) + diag(dim)*0.00002 
        
      }
      temp_acc = mean(accept[1:i])
      if (temp_acc < lower_accept){
        scale = scale/adjust_scale_by
      } else if (temp_acc > upper_accept) {
        #print(temp_acc)
        scale = scale*adjust_scale_by
      }
      sigma = sigma*scale
      print(sigma)
      ct = 0
    }
    mean_acc = cbind(mean_acc, mean(accept[1:i]))
    scale_ans = cbind(scale_ans, scale)
  }
  mean_acc = mean_acc[-1]
  scale_ans = scale_ans[-1]
  
  out = list(batch = batch, accept_rate = mean(accept), scale = scale_ans, sigma = sigma)
  return(out)
}


constrKNGSandwich = function(init, ep, tau, sumX, X, Y, nbatch = 1000, scale = 1e-4,
                             lowerbeta, upperbeta, check_data = NULL, nonneg = FALSE, 
                             method = c("koenker", "currentdata", "newdata"),
                             lower_accept = 0.2, upper_accept = 0.5, 
                             update_after = 10, adjust_scale_by = 2){
  if (is.null(check_data)){
    if (method == "newdata"){
      message("No data input. Method changed to Koenker.")
    } else {
      check_data = X
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
    ans = -(ep/2) * max(abs(-tau*sumX + t(X)%*%lessEq)) / ((1-tau)*2*1*max(X)) 
    #+ (-1/2)*(beta%*%beta)
    return(ans)
  }
  
  out = constrMetropSandwich(logA = logA, init = init, nbatch = nbatch, scale = scale,
                             lowerbeta = lowerbeta, upperbeta = upperbeta, 
                             check_data = check_data, Y = Y, nonneg = nonneg, method = method,
                             lower_accept = lower_accept, upper_accept = upper_accept,
                             update_after = update_after, adjust_scale_by = adjust_scale_by)
  return(out)
}

#differnt lower and upper scale helps find scale easier (due to the long tail)
#scale needs to be a vector of the same order as tau
#quantile vector does not need to be in an order
sandwichKNG = function(data, total_eps, median_eps = NULL, main_tau_eps = NULL,
                       tau, main_tau, scale = 1e-4, nbatch = 10000, 
                       method = c("koenker", "currentdata", "newdata"), nonneg = FALSE,
                       check_data = NULL, lower_accept = 0.2, upper_accept = 0.5, 
                       update_after = 10, adjust_scale_by = 2, formula = NULL){
  
  main_tau_fac = as.factor(main_tau)
  sandwich_tau = tau[!tau %in% main_tau_fac]
  main_tau_order = tau[tau %in% main_tau_fac]
  
  if(is.null(scale)){
    scale = rep(0, length(tau))
  }
  
  if(is.null(main_tau_eps)){
    main_tau_eps = ifelse(method == "koenker", 0.8, 0.9)
  }
  
  data = as.matrix(data)
  i = ncol(data)
  Y = as.matrix(data[,i])
  #R = max(abs(Y))
  #Y = Y/R
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
  # names(scale) = tau
  # main_tau_scale = scale[names(scale) %in% main_tau_fac]
  # sandwich_scale = scale[!names(scale) %in% main_tau_fac]
  
  if (is.null(formula)) {
    vars = colnames(data)
    formula = paste(vars[i], " ~ .")
  }
  
  out = stepwiseKNG(data = data, total_eps = total_eps*main_tau_eps, median_eps = median_eps, 
                    tau = main_tau_order, scale = scale, nbatch = nbatch, 
                    method = method, nonneg = nonneg, check_data = check_data, 
                    lower_accept = lower_accept, upper_accept = upper_accept, 
                    update_after = update_after, adjust_scale_by = adjust_scale_by,
                    formula = formula)
  b = out[[1]]
  scale_output = out[[2]]
  accept_rate = out[[3]]
  main_tau = sort(main_tau)
  names(scale_output) = main_tau
  names(accept_rate) = main_tau
  colnames(b) = main_tau
  
  eps_sandwich = (1 - main_tau_eps) * total_eps / length(sandwich_tau)
  # if (is.null(sandwich_start_scale)) {
  #   sandwich_start_scale = 0.15/R
  # }
  
  update_tau = main_tau
  sandwich_tau_lower = sort(sandwich_tau[sandwich_tau < 0.5], decreasing = TRUE)
  sandwich_tau_upper = sort(sandwich_tau[sandwich_tau > 0.5])
  if (length(sandwich_tau_lower) > 0){
    sandwich_scale_lower = rep(NA, length(sandwich_tau_lower))
    names(sandwich_scale_lower) = sandwich_tau_lower
    accept_rate_lower = rep(NA, length(sandwich_tau_lower))
    names(accept_rate_lower) = sandwich_tau_lower
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
      out = constrKNGSandwich(init = upperbeta, ep = eps_sandwich, tau = curr_tau, 
                              sumX = sumX, X = X, Y = Y, nbatch = nbatch, scale = scale, 
                              lowerbeta = lowerbeta, upperbeta = upperbeta, 
                              check_data = check_data, nonneg = FALSE, method = method,
                              lower_accept = lower_accept, upper_accept = upper_accept, 
                              update_after = update_after, adjust_scale_by = adjust_scale_by)
      
      b_sandwich = tail(out[[1]], 1)
      sandwich_scale_lower[i] = tail(out[[3]], 1)
      accept_rate_lower[i] = out[[2]]
      b = cbind(b, t(b_sandwich))
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
    accept_rate_upper = rep(NA, length(sandwich_tau_upper))
    names(accept_rate_upper) = sandwich_tau_upper
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
      out = constrKNGSandwich(init = lowerbeta, ep = eps_sandwich, tau = curr_tau, 
                              sumX = sumX, X = X, Y = Y, nbatch = nbatch, scale = scale, 
                              lowerbeta = lowerbeta, upperbeta = upperbeta, 
                              check_data = check_data, nonneg = FALSE, method = method,
                              lower_accept = lower_accept, upper_accept = upper_accept, 
                              update_after = update_after, adjust_scale_by = adjust_scale_by)
      b_sandwich = tail(out[[1]], 1)
      sandwich_scale_upper[i] = tail(out[[3]], 1)
      accept_rate_upper[i] = out[[2]]
      b = cbind(b, t(b_sandwich))
      colnames(b)[ncol(b)] = curr_tau
      if(m == 0){
        b = t(as.matrix(b[, match(update_tau, colnames(b))]))
      } else {
        b = as.matrix(b[, match(update_tau, colnames(b))])
      }
    }
  }
  
  #add case when there are no lower/upper quantiles
  scale_output = c(scale_output, sandwich_scale_upper, sandwich_scale_lower)
  scale_output = scale_output[order(names(scale_output))]
  accept_rate = c(accept_rate, accept_rate_lower, accept_rate_upper)
  accept_rate = accept_rate[order(names(accept_rate))]
  
  return(list(b, scale_output, accept_rate))
  
}


#rm(list = ls())
# n = 5000
# lambda = 0.1
# x1 = rexp(n, lambda)
# x2 = 4 + 3*x1 + rexp(n, lambda)
# data = cbind(x1, x2)
# #something is weird with the constraints
# # scale = 0.03 for stepwise koenker
# out = stepwiseKNG(data = data, total_eps = 0.1, tau = c(0.05, 0.25, 0.5, 0.75, 0.95), method = "currentdata",
#                   scale = 0.05, lower_accept = 0.2, upper_accept = 0.6, nonneg = TRUE, median_eps = 0.4)
# # out1 = originalKNG(data = data, total_eps = 0.1, tau = c(0.05, 0.25, 0.5, 0.75, 0.95),
# #                    scale = 0.005, upper_accept = 1, nonneg = TRUE, formula = "x2 ~ .")
# out
# quantreg::rq(x2 ~ x1, tau = c(0.05, 0.25, 0.5, 0.75, 0.95))
# out1
# out2 = sandwichKNG(data = data, total_eps = 0.1, median_eps = 0.5, main_tau_eps = 0.8,
#                    tau = seq(0.05, 0.95, 0.05), main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95), scale = 0.03,
#                    method = "currentdata", nonneg = TRUE, lower_accept = 0.2, upper_accept = 0.5)
# out2



utility_logit = function(data, inds){
  data = cbind(data, inds)
  mod = glm(inds ~ ., data = data, family = binomial(link = "logit"))
  preds = predict(mod, type = "response")
  score = sum((preds-0.5)^2)/nrow(data)
  return(score)
}

utility_logit_inter = function(data, inds){
  data = cbind(data, inds)
  mod = glm(inds ~ .^2, data = data, family = binomial(link = "logit"))
  preds = predict(mod, type = "response")
  score = sum((preds-0.5)^2)/nrow(data)
  return(score)
}

syndata = function(beta_result, x, select_quantile){
  #  print(head(select_quantile))
  allsyn = x%*%beta_result
  coord = cbind(c(1:nrow(x)), select_quantile)
  ans = allsyn[coord]
  return(ans)
}

originalKNG = function(data, total_eps, tau, nbatch = 10000, scale = 1e-4, 
                       lower_accept = 0.2, upper_accept = 0.5, nonneg = FALSE, 
                       formula = NULL, update_after = 10, adjust_scale_by = 2){
  ep = total_eps/length(tau)
  scale_kng = 0
  i = ncol(data)
  Y = data[,i]
  R = max(abs(Y))
  Y = Y/R
  X = as.matrix(cbind(rep(1, nrow(data)), data))
  X = as.matrix(X[, -ncol(X)])
  sumX = apply(X = X, 2, FUN = sum)
  m = ncol(X) - 1
  
  accept_rate = rep(NA, length(tau))
  scale_output = rep(NA, length(tau))
  ans = NA
  for (i in 1:length(tau)){
    print(tau[i])
    curr_scale = scale[i]
    #nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = tau[i])
    temp = KNG(init = rep(1, m+1), ep = ep, tau = tau[i], sumX = sumX, X = X, 
               Y = Y, nbatch = nbatch, scale = scale, nonneg = nonneg,
               lower_accept = lower_accept, upper_accept = upper_accept, 
               update_after = update_after, adjust_scale_by = adjust_scale_by)
    ans = cbind(ans, t(tail(temp[[1]], 1)))
    scale_output[i] = tail(temp[[3]], 1)
    accept_rate[i] = temp[[2]]
  }
  ans = ans[, -1]
  return(list(ans, scale_output, accept_rate))
}

# beta_res = list()
# beta_res[[1]] = matrix(1:6, nrow = 2, byrow = TRUE)
# beta_res[[2]] = matrix(seq(10, 60, 10), nrow = 2, byrow = TRUE)
# lapply(beta_res, syndata, x = Y, select_quantile = c(1:3))
# lapply(beta_res, sum)
# test = mapply(syndata, beta_result = beta_res, x = y, MoreArgs = list(sample(1:3, 3)),SIMPLIFY = FALSE)
# test = mapply(cbind, y, test, SIMPLIFY = FALSE)


library(quantreg)
library(reshape2)
library(ggplot2)
t = 1
set.seed(t)
reps = 1
ut_logit = matrix(NA, nrow = 6, ncol = reps) #6
ut_logit_inter = matrix(NA, nrow = 6, ncol = reps)
main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
tau = c(seq(0.05, 0.95, 0.05), 0.99)

ep = 0.25
n = 5000
runs = 10000
lambda = 0.1
a0 = 4
a1 = 3 
b0 = 10
b1 = 2
b2 = 1

#add pdf
#why ggplot is not plotting
#what is this error No id variables; using all as measure variables
par(mfrow = c(1,1))
scale_x1 = rep(list(c(rep(0.01, 10), rep(0.03, 10))), 5)
accept_x1 = rep(list(rep(0, length(tau))), 5)
scale_x2 = rep(list(c(rep(0.0005, 9), 0.001, rep(0.0005, 10))), 5)
accept_x2 = rep(list(rep(0, length(tau))), 5)
scale_x3 = rep(list(c(rep(0.0001, 10), rep(0.0001, 10))), 5)
accept_x3 = rep(list(rep(0, length(tau))), 5)
for (j in 1:reps){
  print(j)
  x1 = rexp(n, lambda)
  x2 = a0 + b0*x1 + rexp(n, lambda)
  #x3 = a1 + b1*x1 + b2*x2 + rexp(n, lambda)
  #vars = c("x1", "x2", "x3")
  vars = c("x1", "x2")
  
  for (k in 1:length(vars)){
    syn_var = vars[k]
    print(syn_var)
    if (syn_var == "x1"){
      fml = "x1 ~ 1"
      data = as.data.frame(x1)
      all_beta = list()
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = runs,
                         scale = 0.05, lower_accept = 0, upper_accept = 1,
                         nonneg = FALSE, formula = fml, update_after = 10, 
                         adjust_scale_by = 2)
      all_beta[[1]] = temp[[1]]
      accept_x1[[1]] = temp[[3]]
      scale_x1[[1]] = temp[[2]]
      
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 1/length(tau), 
                         tau = tau, scale = 0.03, nbatch = runs, method = "koenker", 
                         nonneg = TRUE, lower_accept = 0.2, upper_accept = 0.5, 
                         update_after = 10, adjust_scale_by = 2, formula = fml)
      all_beta[[2]] = temp[[1]]
      scale_x1[[2]] = temp[[2]]
      accept_x1[[2]] = temp[[3]]
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 1/length(tau), 
                         tau = tau, scale = 0.05, nbatch = runs, method = "currentdata", 
                         nonneg = TRUE, lower_accept = 0.2, upper_accept = 0.5, 
                         update_after = 10, adjust_scale_by = 2, formula = fml)
      
      all_beta[[3]] = temp[[1]]
      scale_x1[[3]] = temp[[2]]
      accept_x1[[3]] = temp[[3]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 1/length(main_tau),
                         main_tau_eps = length(main_tau)/length(tau), tau = tau, 
                         main_tau = main_tau, scale = 0.03, nbatch = runs, method = "koenker", 
                         nonneg = TRUE, lower_accept = 0.2, upper_accept = 0.5, 
                         update_after = 10, adjust_scale_by = 2, formula = fml)
      
      all_beta[[4]] = temp[[1]]
      scale_x1[[4]] = temp[[2]]
      accept_x1[[4]] = temp[[3]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 1/length(main_tau),
                         main_tau_eps = length(main_tau)/length(tau), tau = tau, 
                         main_tau = main_tau, scale = 0.05, nbatch = runs, method = "currentdata", 
                         nonneg = TRUE, lower_accept = 0.2, upper_accept = 0.5, 
                         update_after = 10, adjust_scale_by = 2, formula = fml)
      all_beta[[5]] = temp[[1]]
      scale_x1[[5]] = temp[[2]]
      accept_x1[[5]] = temp[[3]]
      all_beta[[6]] = coef(rq(x1 ~ 1, tau)) # change back to 6
      
      
      
      X = rep(list(matrix(1, nrow = n)), 6) # change back to 6
      #X = rep(list(matrix(1, nrow = n)), 1)
      synx1 = mapply(syndata, beta_result = all_beta, x = X, 
                     MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
      synall = synx1
      synx1[[7]] = x1 # change back to 7
      # names(synx1) = c("Original KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope", 
      #                  "Non-Private", "Raw Data")
      names(synx1) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                       "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
                       "Raw Data")
      # synx1[[2]] = x1
      # names(synx1) = c("Original", "Truth")
      plotdata = melt(synx1)
      plotdata$L1 = as.factor(plotdata$L1)
      print(ggplot(plotdata, aes(x=value, fill = L1)) + facet_wrap(~L1, ncol = 4) +
              stat_density(geom = "area", bw = 5, alpha = 0.5, size = 1) +
              scale_color_brewer(palette="Dark2") + theme_minimal() +
              theme(legend.position=c(0.9,0.2)) +
              ggtitle(paste("Density of variable", syn_var, "- Rep", j)))
      
    } else {
      if (syn_var == "x2"){
        data = as.data.frame(cbind(x1, x2))
        mod = "x2 ~ x1"
      } else {
        data = as.data.frame(cbind(x1, x2, x3))
        mod = "x3 ~ x1 + x2"
      }
      all_beta = list()
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = 10000,
                         scale = 1, 
                         lower_accept = 0, upper_accept = 1,
                         nonneg = FALSE, formula = mod, update_after = 10, 
                         adjust_scale_by = 2)
      all_beta[[1]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[1]] = temp[[2]]
        accept_x2[[1]] = temp[[3]]
      } else {
        scale_x3[[1]] = temp[[2]]
        accept_x3[[1]] = temp[[3]]
      }
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.5, 
                         tau = tau, scale = 3, nbatch = 10000, method = "koenker", 
                         nonneg = TRUE, lower_accept = 0.2, upper_accept = 0.5, 
                         update_after = 10, adjust_scale_by = 2, formula = mod)
      all_beta[[2]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[2]] = temp[[2]]
        accept_x2[[2]] = temp[[3]]
      } else {
        scale_x3[[2]] = temp[[2]]
        accept_x3[[2]] = temp[[3]]
      }
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.4, 
                         tau = tau, scale = 3, nbatch = 10000, method = "currentdata", 
                         nonneg = TRUE, lower_accept = 0.2, upper_accept = 0.5, 
                         update_after = 10, adjust_scale_by = 2, formula = mod)
      all_beta[[3]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[3]] = temp[[2]]
        accept_x2[[3]] = temp[[3]]
      } else {
        scale_x3[[3]] = temp[[2]]
        accept_x3[[3]] = temp[[3]]
      }
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.7, main_tau_eps = 0.8,
                         tau = tau, main_tau = main_tau, scale = 3, nbatch = runs, 
                         method = "koenker", nonneg = TRUE, lower_accept = 0.2, 
                         upper_accept = 0.5, update_after = 10, adjust_scale_by = 2, 
                         formula = mod)
      
      all_beta[[4]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[4]] = temp[[2]]
        accept_x2[[4]] = temp[[3]]
      } else {
        scale_x3[[4]] = temp[[2]]
        accept_x3[[4]] = temp[[3]]
      }
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.4, main_tau_eps = 0.9,
                         tau = tau, main_tau = main_tau, scale = 3, nbatch = runs, 
                         method = "currentdata", nonneg = TRUE, lower_accept = 0.2, 
                         upper_accept = 0.5, update_after = 10, adjust_scale_by = 2, 
                         formula = mod)
      all_beta[[5]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[5]] = temp[[2]]
        accept_x2[[5]] = temp[[3]]
      } else {
        scale_x3[[5]] = temp[[2]]
        accept_x3[[5]] = temp[[3]]
      }
      
      all_beta[[6]] = coef(rq(mod, tau)) # change back to 6
      X = mapply(cbind, rep(list(matrix(1, nrow = n)), 6), synall, SIMPLIFY = FALSE) # 6
      #X = mapply(cbind, synall, rep(list(matrix(1, nrow = n)), 1), SIMPLIFY = FALSE)
      syn = mapply(syndata, beta_result = all_beta, x = X, 
                   MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
      synall = mapply(cbind, synall, syn, SIMPLIFY = FALSE)
      syn[[7]] = data[, ncol(data)] #change to 7
      # names(syn) = c("Original KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope", 
      #                  "Non-Private", "Raw Data")
      names(syn) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                     "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
                     "Raw Data")
      # syn[[2]] = data[, ncol(data)]
      # names(syn) = c("Original", "Truth")
      plotdata = melt(syn)
      plotdata$L1 = as.factor(plotdata$L1)
      print(ggplot(plotdata, aes(x=value, fill = L1)) + facet_wrap(~L1, ncol = 4) +
              stat_density(geom = "area", bw = ifelse(syn_var == "x2", 40, 80), alpha = 0.5, size = 1) + 
              scale_color_brewer(palette="Dark2") + theme_minimal() +
              coord_cartesian(xlim=c(-50, 1000)) +
              theme(legend.position=c(0.9,0.2)) +
              ggtitle(paste("Density of variable", syn_var, "- Rep", j)) +
              labs(fill="Methods")) 
      
    }
    
  }
  
  synall_name = lapply(synall, `colnames<-`, vars)
  inds = c(rep(0,n), rep(1,n))
  ut.data = lapply(synall_name, rbind, data)
  
  ut_logit[, j] = sapply(ut.data, utility_logit, inds)
  ut_logit_inter[, j] = t(sapply(ut.data, utility_logit_inter, inds))
}

tab3 = cbind(ut_logit, ut_logit_inter)
rownames(tab3) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope", 
                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
colnames(tab3) = c("W/out interaction", "W/ interaction")
tab3

tab1 = apply(ut_logit, 1, quantile)
colnames(tab1) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope", 
                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
tab1
rownames(tab1) = c("Min", "Q1", "Median", "Q3", "Max")

library(gridExtra)
png(paste("utility_", t, "_edited.png", sep =""), height = 30*nrow(tab1), 
    width = 150*ncol(tab1))
grid.table(round(tab1, 5))
dev.off()




tab2 = apply(ut_logit_inter, 1, quantile)
colnames(tab2) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope", 
                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
tab2
rownames(tab2) = c("Min", "Q1", "Median", "Q3", "Max")
png(paste("utility_inter_", t, "_edited.png", sep =""), height = 30*nrow(tab2), 
    width = 150*ncol(tab2))
grid.table(round(tab2, 5))
dev.off()

filename = paste("output/utility_", t, "_edited.Rdata", sep = "")
save(list = c("ut_logit", "ut_logit_inter", "tab1", "tab2"), file = filename)
