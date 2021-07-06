# This code is used to generate data for the Really Useful Synthetic Data Framework.
# It consists of 10 training * 10 synthesizer * 10 different syn dataset with
# 2 variables x1 = rexp(0.1) and x2 = 4 + 3x1 + rexp(0.1).

# Output: An Rdata file consists of 1000 synthetic data for each method.
# Analysis: Use code in useful_syn_data_quality.R to analyze the output from this code.

rm(list = ls())
#library(MCMCpack)
#library(matrixcalc)
#library(Matrix)
#library(scales)

densitySamp = function(synParam, ep, origData, nSyn, nVar, nObs, lambda){
  # synParam[1] = ifelse(synParam[1] <= 0, 1e-3, synParam[1])
  # synParam[length(synParam)] = ifelse(synParam[length(synParam)] <= 0, 1e-3, synParam[length(synParam)]) 
  synParam[1] = abs(synParam[1])
  synParam[length(synParam)] = abs(synParam[length(synParam)])
  util = getPMSE_multivar(origData, synParam, nSyn, nVar, nObs)
  
  return(- (sum(synParam ^ 2 / (2 * lambda))) - (nObs * ep * util[1] / 2)) ## multivariate normal prior, diag. sigma, fixed var.
}

getPMSE_multivar = function(origData, synParam, nSyn, nVar, nObs){
  # paramCount = 0
  # for(a in 1:nVar){
  #   if(synParam[paramCount + a + 1] <= 0){
  #     return(Inf)
  #   } else{
  #     paramCount = paramCount + a + 1
  #   }
  # }
  
  #abs(synParam[length(synParam)])
  synX = vector("list", nSyn)
  pmse = rep(NA, nSyn)
  for(a in 1:nSyn){
    synX[[a]] = syn_data(synParam, nVar, nObs)
    combData = data.frame(cbind(rbind(origData, synX[[a]]), rep(0:1, each = nrow(origData))))
    combData[, nVar + 1] = as.factor(combData[, nVar + 1])
    colnames(combData) = c(paste("var", 1:nVar, sep = ""), "Y")
    
    # subsampling option
    #nTrees = 50
    #subSamp = sample(1:nObs, nObs, replace = F)
    #subSamp2 = sample((nObs + 1):(nObs * 2), nObs, replace = F)
    #subCart = vector("list", nTrees)
    #predProb = matrix(NA, nrow = nrow(combData), ncol = nTrees)
    #for(b in 1:nTrees){
    #    subCart[[b]] = rpart(Y ~ ., data = combData[c(subSamp[(1 + (b - 1) * (nObs / nTrees)):(b *  (nObs / nTrees))], 
    #                                                  subSamp2[(1 + (b - 1) *  (nObs / nTrees)):(b *  (nObs / nTrees))]), 
    #                                                c(sample(1:nVar, 1), ncol(combData)), drop = F], 
    #                         method = "class", control = rpart.control(maxdepth = 1, cp = 0, minsplit = 2, minbucket = 1))
    #    predProb[, b] = predict(subCart[[b]], newdata = combData)[, 2]
    #}
    #pmse[a] = mean(apply(predProb, 2, function(vec){ mean((vec - 0.5) ^ 2) }))
    
    # nonsubsampling
    testUtil = rpart(Y ~ ., data = combData, method = "class", 
                     control = rpart.control(maxdepth = 5, cp = 0.01, minsplit = 2, minbucket = 1, maxsurrogate = 0))
    testProp = predict(testUtil)[, 2]
    pmse[a] = sum((testProp - 0.5) ^ 2) / nrow(combData)
  }
  output = c(mean(pmse), var(pmse))
  return(output)
}

syn_data = function(param, nVar, nObs){
  x = matrix(NA, nrow = nObs, ncol = nVar)
  paramCount = 0
  for(a in 1:nVar){
    if(a == 1){
      x[, 1] = rexp(nObs, param[1])
    } else{
      x[, a] = syn_exp(param[(1 + paramCount):(paramCount + a + 1)], x[, 1:(a - 1), drop = F], nObs)
    }
    paramCount = paramCount + a
  }
  return(x)
}

syn_exp = function(param, pred_mat, nObs){
  output = cbind(1, pred_mat) %*% param[-length(param)] + rexp(nObs, param[length(param)])
  return(output)
}

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

originalKNG = function(data, total_eps, tau, nbatch = 10000, scale = 1e-4, 
                       lower_accept = 0, upper_accept = 1, 
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
syndata = function(beta_result, x, select_quantile){
  #  print(head(select_quantile))
  allsyn = x%*%beta_result
  coord = cbind(c(1:nrow(x)), select_quantile)
  ans = allsyn[coord]
  return(ans)
}


# reps = 2
# toteps = 1 #0.25
# n = 5000
# runs = 1000

# x1 = rexp(n, l)
# x2 = a0 + b0*x1 + rexp(n, l)
#set.seed(1)
# main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
# tau = c(seq(0.05, 0.95, 0.05), 0.99)
# ut_logit = matrix(NA, nrow = 7, ncol = reps) #6
# ut_logit_inter = matrix(NA, nrow = 7, ncol = reps)

#ans_pmse = NULL
compare_methods = function(data1, holdout_dat, tau = c(seq(0.05, 0.95, 0.05), 0.99), 
                           main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99), 
                           n = 5000, toteps = 1, runs = 1000, t){
  
  vars = c("x1", "x2")
  ep = toteps/2
  
  check = 0
  while(check == 0){
    pmse = mcmc::metrop(densitySamp, initial = c(1e-4, 4, 3, 1e-4), nbatch = 100, scale = 0.4, 
                        ep = toteps, origData = data1, nSyn = 25, nVar = 2, nObs = n, lambda = 100000)
    check = pmse$accept
  }
  pmse_res = tail(pmse$batch, 1)
  ans_pmse = pmse_res
  pmse_res[, c(1, 4)] = abs(pmse_res[, c(1, 4)])
  
  
  for (k in 1:length(vars)){
    syn_var = vars[k]
    print(syn_var)
    if (syn_var == "x1"){
      fml = "x1 ~ 1"
      data = as.data.frame(data1[, 1])
      colnames(data) = "x1"
      all_beta = list()
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = runs,
                         scale = 0.01, formula = fml)
      all_beta[[1]] = temp[[1]]
      
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 1/length(tau), 
                         tau = tau, scale = 0.03, nbatch = runs, method = "fixed", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[2]] = temp[[1]]
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 1/length(tau), 
                         tau = tau, scale = 0.05, nbatch = runs, method = "varying", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[3]] = temp[[1]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 1/length(main_tau),
                         main_tau_eps = length(main_tau)/length(tau), tau = tau, 
                         main_tau = main_tau, scale = 0.1, sw_scale = 0.03,
                         nbatch = runs, method = "fixed", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[4]] = temp[[1]]
      
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 1/length(main_tau),
                         main_tau_eps = length(main_tau)/length(tau), tau = tau, 
                         main_tau = main_tau, scale = 0.1, sw_scale = 0.03,
                         nbatch = runs, method = "varying", lb = 0, ub = 1000, 
                         formula = fml)
      all_beta[[5]] = temp[[1]]
      
      all_beta[[6]] = coef(rq(x1 ~ 1, tau,  data = data))
      
      X = rep(list(matrix(1, nrow = n)), 6)
      synall = rep(list(rep(list(), 10)), 7)
      for (l in 1:10){
        synx1 = mapply(syndata, beta_result = all_beta, x = X,
                       MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
        for (k in 1:6){
          synall[[k]][[l]] = synx1[[k]]
        }
      }
      
    } else {
      
      data = as.data.frame(data1)
      mod = "x2 ~ x1"
      
      all_beta = list()
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = 10000,
                         scale = 0.001, formula = mod)
      all_beta[[1]] = temp[[1]]
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.6, #0.003
                         tau = tau, scale = 0.03, nbatch = 10000, method = "fixed", 
                         lb = 0, ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod,
                         check_data = synall[[2]][[1]])
      all_beta[[2]] = temp[[1]]
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.8, #0.003
                         tau = tau, scale = 5e-5, nbatch = 10000, method = "varying", 
                         lb = 0, ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[3]][[1]])
      all_beta[[3]] = temp[[1]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.6, main_tau_eps = 0.7,
                         tau = tau, main_tau = main_tau, scale = 0.06, sw_scale = 0.02,
                         nbatch = runs, method = "fixed", lb = 0, check_data = synall[[4]][[1]],
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod)
      
      all_beta[[4]] = temp[[1]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.8, main_tau_eps = 0.8,
                         tau = tau, main_tau = main_tau, scale = 0.00005, sw_scale = 0.00001,
                         nbatch = runs, method = "varying", lb = 0, 
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[5]][[1]])
      all_beta[[5]] = temp[[1]]
      
      all_beta[[6]] = coef(rq(mod, tau, data = data)) # change back to 6
      
      for (l in 1:6){
        X = mapply(cbind, rep(list(matrix(1, nrow = n)), 10), synall[[l]], SIMPLIFY = FALSE)
        syn = lapply(X, syndata, beta_result = all_beta[[l]], 
                     select_quantile = sample(1:length(tau), n, replace = TRUE))
        synall[[l]] = mapply(cbind, synall[[l]], syn, SIMPLIFY = FALSE)
        
      }
      
      
    }
    
  }
  
  for (l in 1:10){
    synall[[7]][[l]] = syn_data(pmse_res, 2, n)
  }
  
  for (l in 1:7){
    synall[[l]] = lapply(synall[[l]], `colnames<-`, vars)
  }
  
  synall[[8]] = data1
  colnames(synall[[8]]) = vars
  
  synall[[9]] = holdout_dat
  colnames(synall[[9]]) = vars
  
  return(synall)
}

library(doParallel)
num_cores=5#detectCores()-1 #use all available core but 1

workers=makeCluster(num_cores,type="SOCK",outfile="log.txt")
registerDoParallel(workers)
oper <- foreach(i=1:10, .combine=rbind, .multicombine=TRUE, 
                .init=list()) %dopar% {
                  source("code/functions_final.R") #change source before uploading to cluster
                  library(rmutil)
                  library(rpart)
                  library(mcmc)
                  library(quantreg)
                  library(reshape2)
                  library(ggplot2)
                  set.seed(i)
                  l = 0.1
                  a0 = 4
                  b0 = 3
                  n = 5000
                  data1 = syn_data(param = c(l, a0, b0, l), nVar = 2, nObs = n*2)
                  holdout = sample(1:nrow(data1), n)
                  holdout_dat = data1[holdout, ]
                  data1 = data1[-holdout, ]
                  vars = c("x1", "x2")
                  colnames(data1) = vars
                  out = list()
                  for (j in 1:10){ #change end
                    out[[j]] = compare_methods(data1 = data1, holdout_dat = holdout_dat,
                                               t = i, toteps = 1)
                  }
                  out
                }

stopCluster(workers)

save(oper, file = "./output/useful_syn_data_eps1.Rdata")

##################################################################################
# This section is for calculating the utility for the synthetic data generated above.

n = 5000
inds = c(rep(0,n), rep(1,n))
tab_logit = NULL
tab_logit_inter = NULL
for (i in 1:25){
  df = oper[i,]
  data = df[[8]]
  ut.data = lapply(oper[2,], rbind, data)
  ut.data = lapply(ut.data, as.data.frame)
  tab_logit = cbind(tab_logit, sapply(ut.data, utility_logit, inds))
  tab_logit_inter = cbind(tab_logit_inter, sapply(ut.data, utility_logit_inter, inds))
}

tab = apply(tab_logit, 1, quantile)
tab = tab[, -8]
colnames(tab) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                  "Sandwich Fixed Slope", "Sandwich Varying Slope",
                  "NonPrivate", "pMSE Mechanism")
tab = t(round(tab, 5))
colnames(tab) = c("Min", "Q1", "Median", "Q3", "Max")

library(gridExtra)
png(paste("./plot/utility_pMSE_seed1to25.png", sep =""), height = 35*nrow(tab),
    width = 100*ncol(tab))
grid.table(tab)
dev.off()


tab = apply(tab_logit_inter, 1, quantile)
tab = tab[, -8]
colnames(tab) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                  "Sandwich Fixed Slope", "Sandwich Varying Slope",
                  "NonPrivate", "pMSE Mechanism")
tab = t(round(tab, 5))
colnames(tab) = c("Min", "Q1", "Median", "Q3", "Max")
png(paste("./plot/utility_inter_pMSE_seed1to25.png", sep =""), height = 35*nrow(tab),
    width = 100*ncol(tab))
grid.table(tab)
dev.off()









# n = 5000
# inds = c(rep(0,n), rep(1,n))
# df = oper[2,]
# data = df[[8]]
# ut.data = lapply(oper[2,], rbind, data)
# ut.data = lapply(ut.data, as.data.frame)
# sapply(ut.data, utility_logit, inds)
# 
# round(sapply(ut.data, utility_logit_inter, inds), 4)

###########################################################################################
# dpPmseFunc = function(x_orig, nVar, nObs, totEps, batchSize){
#   nSyn = 25; lambda = 100000
#   
#   # eps 1
#   #parScale = 0.35 ## unlimited
#   parScale = 0.4 ## five cut
#   #parScale = 1 ## one cut
#   #parScale = 0.65 ## two cut
#   # eps 0.5
#   #parScale = 0.5 ## unlimited
#   parScale = 0.55 ## five cut
#   #parScale = 1.75 ## one cut
#   #parScale = 1 ## two cut
#   ## eps 0.25
#   #parScale = 1.5 ## unlimited
#   parScale = 1.6 ## five cut
#   #parScale = 2.25 ## one cut
#   #parScale = 1.75 ## two cut
#   
#   #setScale = c(0.6, 0.75)
#   #for(d in 2:nVar){
#   #    setScale = c(setScale, c(1.5, rep(1, (d - 1)), 1.5))
#   #}
#   temp = metrop(densitySamp, init = c(2, 10, -2.5, 0.5, 3), nbatch = batchSize, scale = parScale, 
#                 ep = totEps / 10, origData = x_orig, nSyn = nSyn, nVar = nVar, nObs = nObs, lambda = lambda)
#   #tempFollow = metrop(densitySamp, init = temp$batch[250, ], nbatch = 250, scale = parScale, 
#   #                    ep = totEps / 10, origData = x_orig, nSyn = nSyn, nVar = nVar, nObs = nObs, lambda = lambda)
#   #test = MCMC.parallel(densitySamp, n = 10, init = c(2, 10, -2.5, 0.5, 3), n.chain = 4, n.cpu = 2, packages = c("synthpop", "rpart"), scale = parScale, adapt = F, 
#   #              ep = totEps / 10, origData = x_orig, nSyn = nSyn, nVar = nVar, nObs = nObs, lambda = lambda)
#   #temp = metrop(densitySamp, init = c(2, 10, -2.5, 0.5, 3, 0.75, 1, 2, 1), nbatch = 100, scale = parScale, 
#   #              ep = totEps / 10, origData = x_orig, nSyn = nSyn, nVar = nVar, nObs = nObs, lambda = lambda)
#   
#   dpParam = temp$batch[sample(1:batchSize, 10, replace = T), ]
#   #colMeans(dpParam)
#   #colMeans(temp$batch)
#   
#   dpSyn = list("data", "analyses", "mcoef", "mvar")
#   dpSyn$data = dpSyn$analyses = vector("list", 10)
#   dpSyn$mcoef = dpSyn$mvar = matrix(NA, nrow = 10, ncol = nVar)
#   for(i in 1:10){
#     dpSyn$data[[i]] = syn_data(dpParam[i, ], nVar, nObs)
#     dpSyn$analyses[[i]] = lm(dpSyn$data[[i]][, 2] ~ dpSyn$data[[i]][, 1])
#     #dpSyn$analyses[[i]] = lm(dpSyn$data[[i]][, 3] ~ dpSyn$data[[i]][, 1] + dpSyn$data[[i]][, 2])
#     dpSyn$mcoef[i, ] = dpSyn$analyses[[i]]$coef
#     dpSyn$mvar[i, ] = summary(dpSyn$analyses[[i]])$coefficients[, 2] ^ 2
#   }
#   dpPmseOutput = cbind(colMeans(dpSyn$mcoef), sqrt(colMeans(dpSyn$mvar)))
#   #dpSyn$mcoef
#   #dpSyn$mvar
#   
#   return(list("output" = dpPmseOutput, "accept" = temp$accept, "time" = temp$time))
# }



# ep = 1 #0.25
# n = 5000
# runs = 10000
# l = 0.1
# a0 = 4
# a1 = 3 
# b0 = 3
# b1 = 2
# b2 = 1
# x1 = rexp(n, l)
# x2 = a0 + b0*x1 + rexp(n, l)
# set.seed(1)
# data = syn_data(param = c(l, a0, b0, l), nVar = 2, nObs = 5000)
# 
# batchSize = 1000
# parScale = 0.4 #1.6
# nVar = 2
# nSyn = 25
# nObs = 5000
# lambda = 100000
# temp = mcmc::metrop(densitySamp, initial = c(runif(1), a0, b0, runif(1)), nbatch = batchSize, scale = parScale, 
#               ep = ep, origData = data, nSyn = nSyn, nVar = nVar, nObs = nObs, lambda = lambda)
# temp$accept
# 
# tail(temp$batch)
# syn = syn_data(abs(temp$batch), 2, 5000)
# colnames(syn) = c("x1", "x2")
# lm(x2~x1, data = as.data.frame(syn))
# 
# ut = rbind(data, syn)
# colnames(ut) = c("x1", "x2")
# utility_logit = function(data, inds){
#   data = cbind(data, inds)
#   mod = glm(inds ~ ., data = data, family = binomial(link = "logit"))
#   preds = predict(mod, type = "response")
#   score = sum((preds-0.5)^2)/nrow(data)
#   return(score)
# }
# utility_logit(as.data.frame(ut), inds = rep(1:0, each = n))
