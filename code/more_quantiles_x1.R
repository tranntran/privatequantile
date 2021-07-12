# This file is similar to more_quantiles.R, but with more quantiles for generating synthetic x1.
# This file is also used to experiment more eps on x1.


rm(list = ls())


# Note: Based on the current set up synParam[1, 4, 8] are rates and need to be
# positive, hence the abs()
densitySamp = function(synParam, ep, origData, nSyn, nVar, nObs, lambda){
  # synParam[1] = ifelse(synParam[1] <= 0, 1e-3, synParam[1])
  # synParam[length(synParam)] = ifelse(synParam[length(synParam)] <= 0, 1e-3, synParam[length(synParam)]) 
  synParam[c(1, 4, 8)] = abs(synParam[c(1, 4, 8)])
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

# Note: this function may not be generalizable beyond 3 variables, depending on
# the simulated data setup. Double check param set up before use.
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


compare_methods = function(data1, holdout_dat = holdout_dat,
                           tau_x1 = c(seq(0.05, 0.47, 0.01), 0.5, seq(0.53, 0.99, 0.01), 0.995),
                           tau = c(seq(0.05, 0.95, 0.05), 0.99),
                           #tau = c(seq(0.01, 0.74, 0.01), seq(0.75, 0.995, 0.05)),
                           #tau = c(seq(0.05, 0.995, 0.005)), 
                           main_tau_x1 = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.97, 0.99, 0.995), 
                           main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99),
                           n = 5000, toteps = 1, runs = 1000, t){
  
  vars = c("x1", "x2", "x3")
  #vars = c("x1", "x2")
  #ep = toteps/length(vars)
  
  check = 0
  while(check == 0){
    pmse = mcmc::metrop(densitySamp, initial = c(1e-4, 4, 3, 1e-4, 3, 2, 1, 1e-4), 
                        nbatch = 100, scale = 0.4, ep = toteps, origData = data1, nSyn = 25,
                        nVar = length(vars), nObs = n, lambda = 100000)
    check = pmse$accept
  }
  pmse_res = tail(pmse$batch, 1)
  ans_pmse = pmse_res
  pmse_res[, c(1, 4, 8)] = abs(pmse_res[, c(1, 4, 8)])
  syn_pmse = syn_data(pmse_res, length(vars), n)
  
  
  for (k in 1:length(vars)){
    syn_var = vars[k]
    print(syn_var)
    if (syn_var == "x1"){
      fml = "x1 ~ 1"
      data = as.data.frame(data1[, 1])
      colnames(data) = "x1"
      all_beta = list()
      temp = originalKNG(data = data, total_eps = 0.5, tau = tau_x1, nbatch = runs,
                         scale = 1e-3, formula = fml)
      all_beta[[1]] = temp[[1]]
      
      temp = stepwiseKNG(data = data, total_eps = 0.5, median_eps = 0.25, 
                         tau = tau_x1, scale = 0.006, change_scale = 0.03, change_quantile = 0.7,
                         nbatch = runs, method = "fixed", lb = 0, ub = 1000, formula = fml)
      all_beta[[2]] = temp[[1]]
      
      
      temp = stepwiseKNG(data = data, total_eps = 0.5, median_eps = 0.25, 
                         tau = tau_x1, scale = 0.007, change_scale = 0.04, 
                         change_quantile = 0.70,
                         nbatch = runs, method = "varying", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[3]] = temp[[1]]
      
      temp = sandwichKNG(data = data, total_eps = 0.5, median_eps = 0.25,
                         main_tau_eps = 0.6, tau = tau_x1, 
                         main_tau = main_tau_x1, scale = 0.1, change_scale = 0.1, 
                         sw_change_quantile = 0.70,
                         nbatch = runs, method = "fixed", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[4]] = temp[[1]]
      
      
      temp = sandwichKNG(data = data, total_eps = 0.5, median_eps = 0.25,
                         main_tau_eps = 0.6, tau = tau_x1, 
                         sw_change_quantile = 0.70,
                         nbatch = runs, method = "varying", lb = 0, ub = 1000, 
                         formula = fml)
      all_beta[[5]] = temp[[1]]
      
      
      all_beta[[6]] = coef(rq(x1 ~ 1, tau_x1,  data = data))
      
      X = rep(list(matrix(1, nrow = n)), 6)
      synx1 = mapply(syndata, beta_result = all_beta, x = X,
                     MoreArgs = list(sample(1:length(tau_x1), n, replace = TRUE)), SIMPLIFY = FALSE)
      synall = synx1
      synx1[[7]] = c(data[,1]) # change back to 7
      synx1[[8]] = syn_pmse[,1]
      names(synx1) = c("KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                       "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
                       "Raw Data", "pMSE Mechanism")
      for (i in 2:5){
        plotdata = melt(synx1[c(i, 6, 7)])
        print(ggplot(plotdata,aes(x=value, fill= L1)) + geom_density(alpha=0.51))
        
      }
      
      
      # plotdata$L1 = factor(plotdata$L1, levels = c("Raw Data", "Non-Private", "pMSE Mechanism", "KNG",
      #                                              "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
      #                                              "Sandwich-Fixed Slope", "Sandwich-Varying Slope"))
      # filename = paste(paste0('seed', t), paste0('e', ep*2), syn_var, sep = '_')
      # png(paste("./plot/simulations/pMSE/",paste0('eps', ep*2) , "/simulation_", filename, ".png", sep = ""),
      #    width = 700, height = 400)
      
      
      # print(ggplot(plotdata, aes(x=value, fill = L1)) + facet_wrap(~L1, ncol = 4) +
      #         stat_density(geom = "area", bw = 5, alpha = 0.5, size = 1) +
      #         scale_color_brewer(palette="Dark2") + theme_minimal() +
      #         theme(legend.position='none') +
      #         labs(fill='Methods')  + coord_cartesian(xlim=c(-50, 100)) +
      #         ggtitle(paste("Density of variable", syn_var, "- Eps", round(ep, 2))))
      # dev.off()
      
    } else {
      
      if (syn_var == "x2"){
        data = data = as.data.frame(data1[, c(1, 2)])
        colnames(data) = c("x1", "x2")
        mod = "x2 ~ x1"
      } else {
        data = as.data.frame(data1)
        mod = "x3 ~ x1 + x2"
      }
      
      all_beta = list()
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = 10000,
                         scale = 0.0001, formula = mod)
      all_beta[[1]] = temp[[1]]
      
      #0.006 #double check this
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.8, #0.003
                         tau = tau, scale = 0.02, change_scale = 0.08, 
                         change_quantile = 0.7, nbatch = 10000, method = "fixed", 
                         lb = 0, ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod,
                         check_data = synall[[2]])
      all_beta[[2]] = temp[[1]]
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.8, #0.003
                         tau = tau,scale = 6e-9, change_scale = 2e-8,  #5e-9
                         change_quantile = 0.6, nbatch = 10000, method = "varying", 
                         lb = 0, ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[3]])
      all_beta[[3]] = temp[[1]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.6, main_tau_eps = 0.7,
                         tau = tau, main_tau = main_tau, scale = 0.2, change_scale = 0.25, 
                         change_quantile = 0.7, sw_scale = 0.04, sw_change_scale = 0.1,
                         sw_change_quantile = 0.7,
                         nbatch = runs, method = "fixed", lb = 0, check_data = synall[[4]],
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod)
      
      all_beta[[4]] = temp[[1]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.8, main_tau_eps = 0.8,
                         tau = tau, main_tau = main_tau, scale = 2e-5, change_scale = 7e-6, #1.5e-5
                         change_quantile = 0.75, sw_scale = 1.5e-6, sw_change_scale = 3e-6, #2e-6
                         sw_change_quantile = 0.65,
                         nbatch = runs, method = "varying", lb = 0, 
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[5]])
      all_beta[[5]] = temp[[1]]
      
      all_beta[[6]] = coef(rq(mod, tau, data = data)) # change back to 6
      
      X = mapply(cbind, rep(list(matrix(1, nrow = n)), 6), synall, SIMPLIFY = FALSE)
      syn = mapply(syndata, beta_result = all_beta, x = X,
                   MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
      synall = mapply(cbind, synall, syn, SIMPLIFY = FALSE)
      
      syn[[7]] = data[,2]
      syn[[8]] = syn_pmse[,2]
      names(syn) = c("KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                     "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
                     "Raw Data", "pMSE Mechanism")
      
      for (i in 2:5){
        plotdata = melt(syn[c(i, 6, 7)])
        print(ggplot(plotdata,aes(x=value, fill= L1)) + geom_density(alpha=0.51))
      }
      
      
      # syn[[7]] = data[, ncol(data)] #change to 7
      # syn[[8]] = syn_pmse[, 2]
      # names(syn) = c("KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
      #                "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
      #                "Raw Data", "pMSE Mechanism")
      # 
      # plotdata = melt(syn)
      # plotdata$L1 = factor(plotdata$L1, levels = c("Raw Data", "Non-Private", "pMSE Mechanism", "KNG",
      #                                              "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
      #                                              "Sandwich-Fixed Slope", "Sandwich-Varying Slope"))
      # 
      # filename = paste(paste0('seed', t), paste0('e', ep*2), syn_var, sep = '_')
      # png(paste("./plot/simulations/pMSE/",paste0('eps', ep*2) , "/simulation_", filename, ".png", sep = ""), 
      #     width = 700, height = 400)
      # print(ggplot(plotdata, aes(x=value, fill = L1)) + facet_wrap(~L1, ncol = 4) +
      #         stat_density(geom = "area", bw = ifelse(syn_var == "x2", 40, 80), alpha = 0.5, size = 1) + 
      #         scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(c(-250, 500))+
      #         #coord_cartesian(xlim=c(-500, 1000)) +
      #         theme(legend.position='none') +
      #         ggtitle(paste("Density of variable", syn_var, "- Eps", round(ep, 2))) +
      #         labs(fill="Methods")) 
      # dev.off()
      
    }
    
  }
  
  synall[[7]] = syn_pmse
  synall[[8]] = data1
  synall[[9]] = holdout_dat
  synall_name = lapply(synall, `colnames<-`, vars)
  # inds = c(rep(0,n), rep(1,n))
  # ut.data = lapply(synall_name, rbind, data)
  
  return(synall_name)
  # ut_logit[, j] = sapply(ut.data, utility_logit, inds)
  # ut_logit_inter[, j] = t(sapply(ut.data, utility_logit_inter, inds))
  # return(list(ut_logit, ut_logit_inter))
}

library(doParallel)
num_cores=detectCores()-1 #use all available core but 1

workers=makeCluster(num_cores,type="SOCK",outfile="log.txt")
registerDoParallel(workers)
oper <- foreach(i=1:25, .combine=rbind, .multicombine=TRUE, 
                .init=list()) %dopar% {
                  source("code/functions_final_scale.R")
                  library(rmutil)
                  library(rpart)
                  library(mcmc)
                  library(quantreg)
                  library(reshape2)
                  library(ggplot2)
                  set.seed(i)
                  # x1 = rexp(0.1)
                  # x2 = 4 + 3x1 + rexp(0.1)
                  # x3 = 3 + 2x1 + x2 + rexp(0.1)
                  n = 5000
                  data1 = syn_data(param = c(0.1, 4, 3, 0.1, 3, 2, 1, 0.1), nVar = 3, nObs = n*2)
                  vars = c("x1", "x2", "x3")
                  colnames(data1) = vars
                  holdout = sample(1:nrow(data1), n)
                  holdout_dat = data1[holdout, ]
                  data1 = data1[-holdout, ]
                  #compare_methods(data1 = data1, holdout_dat = holdout_dat, t = i, toteps = 1)
                  out = compare_methods(data1 = data1, holdout_dat = holdout_dat, t = i, toteps = 1)
                }

stopCluster(workers)

save(oper, file = "./output/pMSE/data_eps1_50q.Rdata")
