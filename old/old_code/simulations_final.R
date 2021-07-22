rm(list = ls())
source("functions_final.R")

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
t = 100
set.seed(t)
reps = 10 #change reps
ut_logit = matrix(NA, nrow = 6, ncol = reps) #6
ut_logit_inter = matrix(NA, nrow = 6, ncol = reps)
main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
tau = c(seq(0.05, 0.95, 0.05), 0.99)

ep = 1/3
n = 5000
runs = 10000
lambda = 0.1
a0 = 4
a1 = 3 
b0 = 3
b1 = 2
b2 = 1

par(mfrow = c(1,1))
scale_x1 = rep(list(c(rep(0.01, 10), rep(0.03, 10))), 5)
accept_x1 = rep(list(rep(0, length(tau))), 5)
scale_x2 = rep(list(c(rep(0.0005, 9), 0.001, rep(0.0005, 10))), 5)
accept_x2 = rep(list(rep(0, length(tau))), 5)
scale_x3 = rep(list(c(rep(0.0001, 10), rep(0.0001, 10))), 5)
accept_x3 = rep(list(rep(0, length(tau))), 5)
all_beta_x1 = list()
all_beta_x2 = list()
all_beta_x3 = list()
for (j in 1:reps){
  print(j)
  x1 = rexp(n, lambda)
  x2 = a0 + b0*x1 + rexp(n, lambda)
  x3 = a1 + b1*x1 + b2*x2 + rexp(n, lambda)
  vars = c("x1", "x2", "x3")
  #vars = c("x1", "x2")
  
  for (k in 1:length(vars)){
    syn_var = vars[k]
    print(syn_var)
    if (syn_var == "x1"){
      fml = "x1 ~ 1"
      data = as.data.frame(x1)
      all_beta = list()
      temp = originalKNG(data = data, total_eps = ep, tau = tau, nbatch = runs,
                         scale = 0.01, formula = fml)
      all_beta[[1]] = temp[[1]]
      accept_x1[[1]] = temp[[3]]
      scale_x1[[1]] = temp[[2]]
      
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 1/length(tau), 
                         tau = tau, scale = 0.03, nbatch = runs, method = "fixed", 
                         lb = 0, ub = 1000, formula = fml)
      all_beta[[2]] = temp[[1]]
      scale_x1[[2]] = temp[[2]]
      accept_x1[[2]] = temp[[3]]
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 1/length(tau), 
                         tau = tau, scale = 0.05, nbatch = runs, method = "varying", 
                         lb = 0, ub = 1000, formula = fml)
      
      all_beta[[3]] = temp[[1]]
      scale_x1[[3]] = temp[[2]]
      accept_x1[[3]] = temp[[3]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 1/length(main_tau),
                         main_tau_eps = length(main_tau)/length(tau), tau = tau, 
                         main_tau = main_tau, scale = 0.03, sw_scale = 0.02,
                         nbatch = runs, method = "fixed", 
                         lb = 0, ub = 1000, formula = fml)
      
      
      all_beta[[4]] = temp[[1]]
      scale_x1[[4]] = temp[[2]]
      accept_x1[[4]] = temp[[3]]
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 1/length(main_tau),
                         main_tau_eps = length(main_tau)/length(tau), tau = tau, 
                         main_tau = main_tau, scale = 0.05, sw_scale = 0.01,
                         nbatch = runs, method = "varying", lb = 0, ub = 1000, 
                         formula = fml)
      all_beta[[5]] = temp[[1]]
      scale_x1[[5]] = temp[[2]]
      accept_x1[[5]] = temp[[3]]
      all_beta[[6]] = coef(rq(x1 ~ 1, tau)) # change back to 6
      
      all_beta_x1[[j]] = all_beta
      
      X = rep(list(matrix(1, nrow = n)), 6) # change back to 6
      #X = rep(list(matrix(1, nrow = n)), 1)
      synx1 = mapply(syndata, beta_result = all_beta, x = X, 
                     MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
      synall = synx1
      synx1[[7]] = x1 # change back to 7
      names(synx1) = c("KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                       "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
                       "Raw Data")
      plotdata = melt(synx1)
      plotdata$L1 = factor(plotdata$L1, levels = c("Raw Data", "Non-Private", "KNG",
                                                   "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                                                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope"))
      filename = paste(paste0('seed', t), paste0('rep', j), paste0('e', ep*3), syn_var, sep = '_')
      png(paste("../plot/simulations/eps1/simulation_", filename, ".png", sep = ""), 
          width = 700, height = 400)
      print(ggplot(plotdata, aes(x=value, fill = L1)) + facet_wrap(~L1, ncol = 4) +
              stat_density(geom = "area", bw = 5, alpha = 0.5, size = 1) +
              scale_color_brewer(palette="Dark2") + theme_minimal() +
              theme(legend.position=c(0.875,0.2)) +
              labs(fill='Methods')  + coord_cartesian(xlim=c(-50, 100)) +
              ggtitle(paste("Density of variable", syn_var, "- Rep", j, "- Eps", round(ep, 2))))
      dev.off()
      
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
                         scale = 0.001, formula = mod)
      all_beta[[1]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[1]] = temp[[2]]
        accept_x2[[1]] = temp[[3]]
      } else {
        scale_x3[[1]] = temp[[2]]
        accept_x3[[1]] = temp[[3]]
      }
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.6, #0.003
                         tau = tau, scale = 0.03, nbatch = 10000, method = "fixed", 
                         lb = 0, ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod,
                         check_data = synall[[2]])
      all_beta[[2]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[2]] = temp[[2]]
        accept_x2[[2]] = temp[[3]]
      } else {
        scale_x3[[2]] = temp[[2]]
        accept_x3[[2]] = temp[[3]]
      }
      
      temp = stepwiseKNG(data = data, total_eps = ep, median_eps = 0.8, #0.003
                         tau = tau, scale = 5e-5, nbatch = 10000, method = "varying", 
                         lb = 0, ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[3]])
      all_beta[[3]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[3]] = temp[[2]]
        accept_x2[[3]] = temp[[3]]
      } else {
        scale_x3[[3]] = temp[[2]]
        accept_x3[[3]] = temp[[3]]
      }
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.6, main_tau_eps = 0.7,
                         tau = tau, main_tau = main_tau, scale = 0.05, sw_scale = 0.01,
                         nbatch = runs, method = "fixed", lb = 0, check_data = synall[[4]],
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod)
      
      all_beta[[4]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[4]] = temp[[2]]
        accept_x2[[4]] = temp[[3]]
      } else {
        scale_x3[[4]] = temp[[2]]
        accept_x3[[4]] = temp[[3]]
      }
      
      temp = sandwichKNG(data = data, total_eps = ep, median_eps = 0.8, main_tau_eps = 0.8,
                         tau = tau, main_tau = main_tau, scale = 0.00005, sw_scale = 0.00001,
                         nbatch = runs, method = "varying", lb = 0, 
                         ub = ifelse(syn_var == "x2", 1000, 2000), formula = mod, 
                         check_data = synall[[5]])
      all_beta[[5]] = temp[[1]]
      if (syn_var == "x2") {
        scale_x2[[5]] = temp[[2]]
        accept_x2[[5]] = temp[[3]]
      } else {
        scale_x3[[5]] = temp[[2]]
        accept_x3[[5]] = temp[[3]]
      }
      
      all_beta[[6]] = coef(rq(mod, tau)) # change back to 6
      
      if (syn_var == "x2"){
        all_beta_x2[[j]] = all_beta
      } else {
        all_beta_x3[[j]] = all_beta
      }
      
      X = mapply(cbind, rep(list(matrix(1, nrow = n)), 6), synall, SIMPLIFY = FALSE) # 6
      #X = mapply(cbind, synall, rep(list(matrix(1, nrow = n)), 1), SIMPLIFY = FALSE)
      syn = mapply(syndata, beta_result = all_beta, x = X, 
                   MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
      synall = mapply(cbind, synall, syn, SIMPLIFY = FALSE)
      syn[[7]] = data[, ncol(data)] #change to 7
      # names(syn) = c("Original KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope", 
      #                  "Non-Private", "Raw Data")
      names(syn) = c("KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                     "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private",
                     "Raw Data")
      # syn[[2]] = data[, ncol(data)]
      # names(syn) = c("Original", "Truth")
      
      plotdata = melt(syn)
      plotdata$L1 = factor(plotdata$L1, levels = c("Raw Data", "Non-Private", "KNG",
                                                   "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                                                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope"))
      
      filename = paste(paste0('seed', t), paste0('rep', j), paste0('e', ep*3), syn_var, sep = '_')
      png(paste("../plot/simulations/eps1/simulation_", filename, ".png", sep = ""), 
          width = 700, height = 400)
      print(ggplot(plotdata, aes(x=value, fill = L1)) + facet_wrap(~L1, ncol = 4) +
              stat_density(geom = "area", bw = ifelse(syn_var == "x2", 40, 80), alpha = 0.5, size = 1) + 
              scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(c(-250, 500))+
              #coord_cartesian(xlim=c(-500, 1000)) +
              theme(legend.position=c(0.875,0.2)) +
              ggtitle(paste("Density of variable", syn_var, "- Rep", j, "- Eps", round(ep, 2))) +
              labs(fill="Methods")) 
      dev.off()
      
    }
    
  }
  
  synall_name = lapply(synall, `colnames<-`, vars)
  inds = c(rep(0,n), rep(1,n))
  ut.data = lapply(synall_name, rbind, data)
  
  ut_logit[, j] = sapply(ut.data, utility_logit, inds)
  ut_logit_inter[, j] = t(sapply(ut.data, utility_logit_inter, inds))
}

# tab4 = cbind(ut_logit[,30], ut_logit_inter[,30])
# rownames(tab4) = c("Original KNG", "Stepwise-Fixed Slope", "Stepwise-Varying Slope", 
#                    "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
# colnames(tab4) = c("W/out interaction", "W/ interaction")
# tab4
# 
# 
# test = ut.data[[1]]
# test = cbind(test, inds)
# mod = glm(inds ~ ., data = test, family = binomial(link = "logit"))
# preds = predict(mod, type = "response")
# score = sum((preds-0.5)^2)/nrow(test)
# score
# 
# tab3 = cbind(ut_logit, ut_logit_inter)
# rownames(tab3) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
#   "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
# colnames(tab3) = c("W/out interaction", "W/ interaction")
# tab3
# 
# tab1 = apply(ut_logit, 1, quantile, na.rm = T)
# colnames(tab1) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
#                    "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
# tab1
# rownames(tab1) = c("Min", "Q1", "Median", "Q3", "Max")
# round(tab1, 4)
# 
# library(gridExtra)
# png(paste("utility_", t, "_edited.png", sep =""), height = 30*nrow(tab1), 
#     width = 150*ncol(tab1))
# grid.table(round(tab1, 5))
# dev.off()
# 
# 
# 

# tab2 = apply(ut_logit_inter, 1, quantile)
# colnames(tab2) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
#                    "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
# tab2
# rownames(tab2) = c("Min", "Q1", "Median", "Q3", "Max")
# png(paste("utility_inter_", t, "_edited.png", sep =""), height = 30*nrow(tab2), 
#     width = 150*ncol(tab2))
# grid.table(round(tab2, 5))
# dev.off()

filename = paste("../output/simulations_seed", t, "_e", ep*3, "_10reps.Rdata", sep = "")
save(list = c("ut_logit", "ut_logit_inter", "all_beta_x1", 
              "all_beta_x2", "all_beta_x3"), file = filename)
