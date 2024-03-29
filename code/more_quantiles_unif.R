# This code does not need any input to run.
# The outout of this code is used in compare_pMSE_utility.R

# This code creates synthetic versions of simulated data using
# KNG methods, pMSE mechanism, and non-private quantile regression.
# The simulated data consists of 3 variables (x1, x2, and x3)
# and 5000 observations. 

# The code for pMSE mechanism is adapted from "pMSE Mechanism: 
# Differentially Private Synthetic Data with Maximal Distributional 
# Similarity" (Snoke and Slavković, 2018).

rm(list = ls())

## pMSE mechanism
# Note: Based on the current set up synParam[1, 4, 8] are rates 
# and need to be positive, hence the abs().
densitySamp = function(synParam, ep, origData, nSyn, nVar, nObs, lambda){
  synParam[c(1, 4, 8)] = abs(synParam[c(1, 4, 8)])
  util = getPMSE_multivar(origData, synParam, nSyn, nVar, nObs)
  return(- (sum(synParam ^ 2 / (2 * lambda))) - (nObs * ep * util[1] / 2)) ## multivariate normal prior, diag. sigma, fixed var.
}

getPMSE_multivar = function(origData, synParam, nSyn, nVar, nObs){
  synX = vector("list", nSyn)
  pmse = rep(NA, nSyn)
  for(a in 1:nSyn){
    synX[[a]] = syn_data(synParam, nVar, nObs)
    combData = data.frame(cbind(rbind(origData, synX[[a]]), rep(0:1, each = nrow(origData))))
    combData[, nVar + 1] = as.factor(combData[, nVar + 1])
    colnames(combData) = c(paste("var", 1:nVar, sep = ""), "Y")
    
    # nonsubsampling
    testUtil = rpart(Y ~ ., data = combData, method = "class", 
                     control = rpart.control(maxdepth = 5, cp = 0.01, minsplit = 2, minbucket = 1, maxsurrogate = 0))
    testProp = predict(testUtil)[, 2]
    pmse[a] = sum((testProp - 0.5) ^ 2) / nrow(combData)
  }
  output = c(mean(pmse), var(pmse))
  return(output)
}

# Note: this function may not be generalizable beyond 3 variables, depending on
# the simulated data setup. Double check param set up before making changes.
syn_data = function(param, nVar, nObs){
  x = matrix(NA, nrow = nObs, ncol = nVar)
  paramCount = 0
  for(a in 1:nVar){
    if(a == 1){
      x[, 1] = rexp(nObs, param[1])
    } else{
      x[, a] = syn_exp(param[(1 + paramCount):(paramCount + a + 1)], x[, 1:(a - 1), drop = F], nObs)
    }
    if (a == 1){
      paramCount = paramCount + a
    } else {
      paramCount = paramCount + a + 1
    }
  }
  return(x)
}

syn_exp = function(param, pred_mat, nObs){
  output = cbind(1, pred_mat) %*% param[-length(param)] + rexp(nObs, param[length(param)])
  return(output)
}


## implementation of KNG (Reimherr and Awan, 2019)
originalKNG = function(data, total_eps, tau, nbatch = 10000, scale = 1e-4, 
                       lower_accept = 0, upper_accept = 1, 
                       formula = NULL, update_after = 10, adjust_scale_by = 2){
  ep = total_eps/length(tau)
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
    curr_scale = scale[i]
    temp = KNG(init = rep(0, ncol(X)), ep = ep, tau = tau[i], sumX = sumX, X = X, 
               Y = Y, nbatch = nbatch, scale = scale/R,
               lower_accept = lower_accept, upper_accept = upper_accept, 
               update_after = update_after, adjust_scale_by = adjust_scale_by)
    ans = cbind(ans, t(tail(temp[[1]], 1)*R))
    scale_output[i] = tail(temp[[3]], 1)
    accept_rate[i] = temp[[2]]
  }
  ans = ans[, -1]
  return(list(ans, scale_output, accept_rate))
}


# this function is used to get synthetic values from
# a vector of randomly selected quantiles
syndata = function(beta_result, x, select_quantile){
  allsyn = x%*%beta_result
  coord = cbind(c(1:nrow(x)), select_quantile)
  ans = allsyn[coord]
  return(ans)
}

# this function is for generating synthetic versions for simulated data
compare_methods = function(data1, holdout_dat = holdout_dat,
                           tau = c(seq(0.01, 0.47, 0.02), 0.5, seq(0.53, 0.99, 0.02)),
                           main_tau = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), 
                           n = 5000, toteps = 1, runs = 1000){
  
  vars = c("x1", "x2", "x3")
  
  ## find synthetic data using pMSE mechnism
  # check is the mcmc acceptance rate and is used to ensure it moves around
  check = 0
  while(check == 0){
    pmse = mcmc::metrop(densitySamp, initial = c(1e-4, 4, 3, 1e-4, 3, 2, 1, 1e-4), 
                        nbatch = 100, scale = 0.4, ep = toteps, origData = data1, nSyn = 25,
                        nVar = length(vars), nObs = n, lambda = 100000)
    check = pmse$accept
  }
  # generate data from the last iteration of the mcmc
  pmse_res = tail(pmse$batch, 1)
  pmse_res[, c(1, 4, 8)] = abs(pmse_res[, c(1, 4, 8)])
  syn_pmse = syn_data(pmse_res, length(vars), n)
  
  ## generate synthetic data using KNG methods and non-private quantile regression
  # we will generate synthetic x1 first, then x2, and finally x3
  for (k in 1:length(vars)){
    syn_var = vars[k]
    print(syn_var)
    if (syn_var == "x1"){
      # we allocate more eps to x1 as it can be more bumpy (due to having on 49 
      # possible values) and it is used to generate x2 and x3
      ep = 0.5
      # set formula
      fml = "x1 ~ 1"
      # subset data
      data = as.data.frame(data1[, 1])
      colnames(data) = "x1"
      all_beta = list()
      
      # generate private quantiles using KNG mechanism
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = runs,
                         scale = 5e-4, formula = fml)
      all_beta[[1]] = temp[[1]]
      
      # generate private quantiles using stepwise fixed slope KNG
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.25, 
                         tau = tau, scale = 0.015, change_scale = 0.07, 
                         change_quantile = 0.75, nbatch = runs, 
                         method = "fixed", lb = 0, ub = 1000, formula = fml)
      all_beta[[2]] = temp[[1]]
      
      # generate private quantiles using stepwise varying slope KNG
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.25, 
                         tau = tau, scale = 0.015, change_scale = 0.07, 
                         change_quantile = 0.75, nbatch = runs, 
                         method = "varying", lb = 0, ub = 1000, formula = fml)
      all_beta[[3]] = temp[[1]]
      
      # generate private quantiles using sandwich fixed slope KNG
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.25,
                         main_tau_eps = 0.6, tau = tau, main_tau = main_tau, 
                         scale = 0.2, change_scale = 0.25, change_quantile = 0.83, 
                         sw_scale = 0.02, sw_change_scale = 0.1,
                         sw_change_quantile = 0.8, nbatch = runs, method = "fixed", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[4]] = temp[[1]]
      
      # generate private quantiles using sandwich varying slope KNG
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.25,
                         main_tau_eps = 0.6, tau = tau, main_tau = main_tau, 
                         scale = 0.22, change_scale = 0.25, change_quantile = 0.8, 
                         sw_scale = 0.03, sw_change_scale = 0.1,
                         sw_change_quantile = 0.80, nbatch = runs, method = "varying", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[5]] = temp[[1]]
      
      # generate non-private quantiles using quantile regression
      all_beta[[6]] = coef(rq(x1 ~ 1, tau,  data = data))
      
      # generate synthetic data using randomly selected quantiles
      # we sample quantiles using different probabilities because
      # the quantiles are not evenly spaced out
      X = rep(list(matrix(1, nrow = n)), 6)
      synx1 = mapply(syndata, beta_result = all_beta, x = X,
                     MoreArgs = list(sample(1:length(tau), n, replace = TRUE, 
                                            prob = c(0.01, rep(0.02, 23), 0.03, 0.03, rep(0.02, 22), 0.03))), 
                     SIMPLIFY = FALSE)
      
      # synall is the final output
      synall = synx1
      
      # plotting synthetic data out can help with tuning the scales
      # synx1[[7]] = c(data[,1])
      # synx1[[8]] = syn_pmse[,1]
      # names(synx1) = c("KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
      #                  "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
      #                  "Raw Data", "pMSE Mechanism")
      # for (i in 1:5){
      #   plotdata = melt(synx1[c(i, 6, 7)])
      #   print(ggplot(plotdata,aes(x=value, fill= L1)) + geom_density(alpha=0.51, bw = 2))
      # 
      # }
      
    } else {
      # for x2 and x3, the eps used is 0.25
      ep = 0.25
      
      # subset data and specify model depending on whether the synthesizing variable 
      # is x2 or x3
      if (syn_var == "x2"){
        data = data = as.data.frame(data1[, c(1, 2)])
        colnames(data) = c("x1", "x2")
        mod = "x2 ~ x1"
      } else {
        data = as.data.frame(data1)
        mod = "x3 ~ x1 + x2"
      }
      
      all_beta = list()
      
      # generate private coef (beta) using the KNG mechanism (Reimherr and Awan, 2019)
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = 10000,
                         scale = ifelse(syn_var == "x2", 1e-5, 1e-5), formula = mod)
      all_beta[[1]] = temp[[1]]
      
      # generate private coef (beta) using stepwise fixed slope KNG
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.8, tau = tau, 
                         scale = ifelse(syn_var == "x2", 0.008, 0.003), 
                         change_scale = ifelse(syn_var == "x2", 0.02, 0.02), 
                         change_quantile = ifelse(syn_var == "x2", 0.8, 0.7), 
                         nbatch = 10000, method = "fixed", lb = 0, 
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod,
                         check_data = synall[[2]])
      all_beta[[2]] = temp[[1]]
      
      # generate private coef (beta) using stepwise varying slope KNG
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.8, tau = tau,
                         scale = ifelse(syn_var == "x2", 2e-9, 1e-10),
                         change_scale = ifelse(syn_var == "x2", 6e-9, 2e-10),  #5e-9
                         change_quantile = ifelse(syn_var == "x2", 0.55, 0.6),
                         nbatch = 10000, method = "varying", lb = 0, 
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[3]])
      all_beta[[3]] = temp[[1]]
      
      # generate private coef (beta) using sandwich fixed slope KNG
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.8, main_tau_eps = 0.8,
                         tau = tau, main_tau = main_tau, 
                         scale = ifelse(syn_var == "x2", 0.22, 0.2), 
                         change_scale = ifelse(syn_var == "x2", 0.2, 0.2), #0.15
                         sw_scale = ifelse(syn_var == "x2", 0.03, 0.03), 
                         sw_change_scale = ifelse(syn_var == "x2", 0.1, 0.1),
                         change_quantile = 0.7, sw_change_quantile = 0.7,
                         nbatch = runs, method = "fixed", lb = 0, check_data = synall[[4]],
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod)
      
      all_beta[[4]] = temp[[1]]
      
      # generate private coef (beta) using sandwich varying slope KNG
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.8, main_tau_eps = 0.8,
                         tau = tau, main_tau = main_tau, 
                         scale = ifelse(syn_var == "x2", 1e-7, 2e-8),
                         change_scale = ifelse(syn_var == "x2", 2e-7, 5e-9), #1.5e-5
                         sw_scale = ifelse(syn_var == "x2", 3e-8, 1e-8), 
                         sw_change_scale = ifelse(syn_var == "x2", 2e-7, 5e-9), #1.5e-6
                         change_quantile = 0.7, sw_change_quantile = 0.65,
                         nbatch = runs, method = "varying", lb = 0, 
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[5]])
      all_beta[[5]] = temp[[1]]
      
      # generate private coef (beta) using quantile regression
      all_beta[[6]] = coef(rq(mod, tau, data = data))
      
      # generate synthetic data using randomly selected quantiles
      # we sample the quantiles with different probability due to them not
      # being evenly spaced
      X = mapply(cbind, rep(list(matrix(1, nrow = n)), 6), synall, SIMPLIFY = FALSE)
      syn = mapply(syndata, beta_result = all_beta, x = X,
                   MoreArgs = list(sample(1:length(tau), n, replace = TRUE, 
                                          prob = c(0.01, rep(0.02, 23), 0.03, 0.03, rep(0.02, 22), 0.03))), 
                   SIMPLIFY = FALSE)
      # bind the data to the final output
      synall = mapply(cbind, synall, syn, SIMPLIFY = FALSE)
      
      # ploting the synthetic data can help tune the scales.
      # syn[[7]] = data[, k]
      # syn[[8]] = syn_pmse[, k]
      # names(syn) = c("KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
      #                  "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
      #                  "Raw Data", "pMSE Mechanism")
      # 
      # for (i in 1:5){
      #   plotdata = melt(syn[c(i, 6, 7)])
      #   print(ggplot(plotdata,aes(x=value, fill= L1)) + geom_density(alpha=0.51))
      # }
      
    }
    
  }
  
  # synall, the final output, is a list of 9. The order of synall is KNG, stepwise
  # fixed slope KNG, stepwise varying slope KNG, sandwich fixed slope KNG, sandwich
  # varying slope kNG, non-private, pMSE mechanism, training data (data used to generate
  # synthetic data), and finally testing (holdout) data
  synall[[7]] = syn_pmse
  synall[[8]] = data1
  synall[[9]] = holdout_dat
  synall_name = lapply(synall, `colnames<-`, vars)
  
  return(synall_name)
}


# parallel code to speed up the run time
library(doParallel)
num_cores=detectCores()-1 # use all available core but 1
workers=makeCluster(num_cores,type="SOCK",outfile="log.txt")
registerDoParallel(workers)

# run 100 reps, the ith rep is run with seed i
oper <- foreach(i=1:100, .combine=rbind, .multicombine=TRUE, 
                .init=list()) %dopar% {
                  source("functions_final_scale.R") #code/functions_final_scale.R
                  library(rmutil)
                  library(rpart)
                  library(mcmc)
                  library(quantreg)
                  library(reshape2)
                  library(ggplot2)
                  set.seed(i)
                  # create training and testing data with 5000 obs each
                  # x1 = rexp(0.1)
                  # x2 = 4 + 3x1 + rexp(0.1)
                  # x3 = 3 + 2x1 + x2 + rexp(0.1)
                  n = 5000
                  data1 = syn_data(param = c(0.1, 4, 3, 0.1, 3, 2, 1, 0.1), nVar = 3, nObs = n*2)
                  vars = c("x1", "x2", "x3")
                  colnames(data1) = vars
                  
                  # randomly sample training and hold out data
                  holdout = sample(1:nrow(data1), n)
                  holdout_dat = data1[holdout, ]
                  data1 = data1[-holdout, ]
                  
                  # get synthetic data for all methods
                  out = compare_methods(data1 = data1, holdout_dat = holdout_dat, toteps = 1)
                }

stopCluster(workers)

# if run on pc
# save(oper, file = "./output/data_eps1_50q_unif.Rdata")

# if run on icds cluster
save(oper, file = "../output/data_eps1_50q_unif.Rdata")