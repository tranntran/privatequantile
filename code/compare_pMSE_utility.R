# Use this code with the output from compare_pMSE_data.R
# The output of this code is simulation_utility_results.Rdata

# This code compare the performance of KNG methods, pMSE mechanism,
# and non-private (using quantile regression) synthetic data. 
# The utility measures are: 
# 1. Wasseinstein Randomization Test (Arnold and Neunhoeffer, 2020)
# 2. pMSE score (Snoke et al, 2018)
# 3. K-marginals (Task et al, 2021), 
# 4. Standardized coefficient difference, 
# 5. RMSE on holdout sample (Arnold and Neunhoeffer, 2020).

rm(list = ls())
load("./output/data_eps1_50q_v2.Rdata")

methods = c("KNG", "StepF", "StepV", "SWF", "SWV", "NP","pMSE", "Raw")
set.seed(6789)

###############################################################################
## WRT
wasserstein_randomization_test <- function(a, b, n_rep = 1000){
  require(transport)
  rand_dist <- NULL
  if(class(a) == "factor"){
    sample_dist <- transport::wasserstein1d(table(a), table(b))
  } else {
    sample_dist <- transport::wasserstein1d(a, b)
  }
  for(i in 1:n_rep){  
    tmp <- c(a, b)
    sel <- sample(1:length(tmp), length(a))
    tmp_a <- tmp[sel]
    tmp_b <- tmp[-sel]
    if(class(a) == "factor"){
      w_dist <- transport::wasserstein1d(table(tmp_a), table(tmp_b))
    } else {
      w_dist <- transport::wasserstein1d(tmp_a, tmp_b)
    }
    
    rand_dist <- c(rand_dist, w_dist)
  }
  
  p <- sum(sample_dist < rand_dist)/n_rep
  
  return(list(dist = rand_dist, sample_distance = sample_dist, p = p))
}

wrt_training = matrix(NA, nrow = 100, ncol = 7)
wrt_testing = matrix(NA, nrow = 100, ncol = 7)
for (i in 1:100){
  print(i)
  training_data = oper[i, 8][[1]]
  testing_data = oper[i, 9][[1]]
  tmp = lapply(oper[i, c(1:7)], wasserstein_randomization_test, b = training_data)
  wrt_training[i, ] = sapply(tmp, function(x) x$sample_distance/median(x$dist))
  
  tmp = lapply(oper[i, c(1:7)], wasserstein_randomization_test, b = testing_data)
  wrt_testing[i, ] = sapply(tmp, function(x) x$sample_distance/median(x$dist))
  
}

apply(wrt_testing, 2, median)
apply(wrt_training, 2, median)

tab_wrt = rbind(apply(wrt_training, 2, mean), apply(wrt_testing, 2, mean))
colnames(tab_wrt) = methods[-8]
rownames(tab_wrt) = c("Training WRT", "Testing WRT")

tab_wrt_se = rbind(apply(wrt_training, 2, sd), apply(wrt_testing, 2, sd))/sqrt(100)
rownames(tab_wrt_se) = paste(c("Training WRT", "Testing WRT"), "SE")
tab_wrt = rbind(tab_wrt, tab_wrt_se)

saveRDS(tab_wrt, file="output/simulation_wrt.rds")
###############################################################################
## utility (pMSE score)

# without interaction
utility_logit = function(data, inds){
  data = cbind(data, inds)
  mod = glm(inds ~ ., data = as.data.frame(data), family = binomial(link = "logit"))
  preds = predict(mod, type = "response")
  score = sum((preds-0.5)^2)/nrow(data)
  return(score)
}


utility_training = matrix(NA, nrow = 100, ncol = 7)
utility_testing = matrix(NA, nrow = 100, ncol = 7)
inds = c(rep(0, 5000), rep(1, 5000))
for (i in 1:100){
  training_data = oper[i, 8][[1]]
  testing_data = oper[i, 9][[1]]
  
  tmp = lapply(oper[i, c(1:7)], rbind, training_data)
  utility_training[i, ] = sapply(tmp, utility_logit, inds)
  
  tmp = lapply(oper[i, c(1:7)], rbind, testing_data)
  utility_testing[i, ] = sapply(tmp, utility_logit, inds)
  
  
}

tab_ut = rbind(apply(utility_training, 2, mean), apply(utility_testing, 2, mean))
tab_ut_se = rbind(apply(utility_training, 2, sd), apply(utility_testing, 2, sd))/sqrt(100)
colnames(tab_ut) = methods[-8]
rownames(tab_ut) = paste(c("Training", "Testing"), "Mean pMSE")
colnames(tab_ut_se) = methods[-8]
rownames(tab_ut_se) = paste(c("Training", "Testing"), "SE")
tab_ut
tab_ut_se
tab_ut = rbind(tab_ut, tab_ut_se)


# with interaction

utility_logit_inter = function(data, inds){
  data = cbind(data, inds)
  mod = glm(inds ~ . ^2, data = as.data.frame(data), family = "binomial")
  preds = predict(mod, type = "response")
  score = sum((preds-0.5)^2)/nrow(data)
  return(score)
}


utility_training = matrix(NA, nrow = 100, ncol = 7)
utility_testing = matrix(NA, nrow = 100, ncol = 7)
inds = c(rep(0, 5000), rep(1, 5000))
for (i in 1:100){
  training_data = oper[i, 8][[1]]
  testing_data = oper[i, 9][[1]]
  
  tmp = lapply(oper[i, c(1:7)], rbind, training_data)
  utility_training[i, ] = sapply(tmp, utility_logit_inter, inds)
  
  tmp = lapply(oper[i, c(1:7)], rbind, testing_data)
  utility_testing[i, ] = sapply(tmp, utility_logit_inter, inds)
  
  
}

tab_ut_inter = rbind(apply(utility_training, 2, mean), apply(utility_testing, 2, mean))
colnames(tab_ut_inter) = methods[-8]
rownames(tab_ut_inter) = paste(c("Training", "Testing"), "Mean pMSE")
tab_ut_inter

tab_ut_inter_se = rbind(apply(utility_training, 2, sd), apply(utility_testing, 2, sd))/sqrt(100)
colnames(tab_ut_inter_se) = methods[-8]
rownames(tab_ut_inter_se) = paste(c("Training", "Testing"), "SE")
tab_ut_inter_se

tab_ut_inter = rbind(tab_ut_inter, tab_ut_inter_se)

###############################################################################
## k-marginals
library(data.table)

nist_score = NULL
for (j in 1:100){
  tab = CJ(x1=0:5,x2=0:5, x3=0:5)
  training_data = oper[j, 8][[1]]
  # x1_check = quantile(training_data[,1], c(0.1, 0.25, 0.5, 0.75, 0.9))
  # x2_check = quantile(training_data[,2], c(0.1, 0.25, 0.5, 0.75, 0.9))
  # x3_check = quantile(training_data[,3], c(0.1, 0.25, 0.5, 0.75, 0.9))
  
  x1_check = quantile(training_data[,1])
  x2_check = quantile(training_data[,2])
  x3_check = quantile(training_data[,3])
  
  methods = c("KNG", "StepF", "StepV", "SWF", "SWV", "NP","pMSE", "Raw")
  for (i in 1:8){
    # tmp = cbind(cut(oper[1, i][[1]][,1], breaks=x1_check, labels = F, include.lowest = T),
    #             cut(oper[1, i][[1]][,2], breaks=x2_check, labels = F, include.lowest = T),
    #             cut(oper[1, i][[1]][,3], breaks=x3_check, labels = F, include.lowest = T))
    tmp = cbind(findInterval(x = oper[j, i][[1]][,1], vec = x1_check, left.open = T),
                findInterval(x = oper[j, i][[1]][,2], vec = x2_check, left.open = T),
                findInterval(x = oper[j, i][[1]][,3], vec = x3_check, left.open = T))
    colnames(tmp) = c("x1", "x2", "x3")
    tmp = as.data.table(tmp)
    tmp_tab = tmp[,.N,by=c("x1", "x2", "x3")]
    setnames(tmp_tab, "N", methods[i])
    tab = merge.data.table(tab, tmp_tab, by = c("x1", "x2", "x3"), all.x = T)
  }
  tab[is.na(tab)] = 0
  colSums(tab)
  raw_score = rep(NA, 7)
  for (i in 1:7){
    var = methods[i]
    tmp = abs(tab[, ..var]-tab[,11])/5000
    raw_score[i] = sum(tmp)
  }
  nist_score = rbind(nist_score, (2-raw_score)/2*1000)
  
}

tab_km = apply(nist_score, 2, mean)

names(tab_km) = methods[-8]
tab_km = rbind(tab_km, apply(nist_score, 2, sd)/sqrt(100))
rownames(tab_km) = c('Mean', 'SE')
tab_km

###############################################################################
## standardized coef difference
# x2 ~ x1
int_training = matrix(NA, ncol = 7, nrow = 100)
sl_training = matrix(NA, ncol = 7, nrow = 100)
int_testing = matrix(NA, ncol = 7, nrow = 100)
sl_testing = matrix(NA, ncol = 7, nrow = 100)

for (i in 1:100){
  training_data = oper[i, 8][[1]]
  testing_data = oper[i, 9][[1]]
  syn = oper[i, c(1:7)]
  
  training_fit = lm(x2 ~ x1, data = as.data.frame(training_data))
  tmp = sapply(syn, function(x) coef(lm(x2 ~ x1, data = as.data.frame(x))))
  tmp_dif = abs(tmp - matrix(rep(coef(training_fit), 7), nrow = 2))/sqrt(diag(vcov(training_fit)))
  int_training[i, ] = tmp_dif[1, ]
  sl_training[i, ] = tmp_dif[2, ]
  
  holdout_fit = lm(x2 ~ x1, data = as.data.frame(testing_data))
  tmp_dif = abs(tmp - matrix(rep(coef(holdout_fit), 7), nrow = 2))/sqrt(diag(vcov(holdout_fit)))
  int_testing[i, ] = tmp_dif[1, ]
  sl_testing[i, ] = tmp_dif[2, ]
}


tab_coef = rbind(apply(int_training, 2, mean), apply(sl_training, 2, mean),
                 apply(int_testing, 2, mean), apply(sl_testing, 2, mean))

colnames(tab_coef) = methods[-8]
rownames(tab_coef) = c("Training Intercept", "Training Slope", "Testing Intercept",
                       "Testing Slope")

tab_coef_se = rbind(apply(int_training, 2, sd), apply(sl_training, 2, sd),
                    apply(int_testing, 2, sd), apply(sl_testing, 2, sd)) / sqrt(100)

colnames(tab_coef_se) = methods[-8]
rownames(tab_coef_se) = paste(c("Training Intercept", "Training Slope", "Testing Intercept",
                                "Testing Slope"), "SE")

tab_coef_x2 = rbind(tab_coef, tab_coef_se)
tab_coef_x2[c(1:2, 5:6),]


# x3 ~ x1 + x2
int_training = matrix(NA, ncol = 7, nrow = 100)
sl1_training = matrix(NA, ncol = 7, nrow = 100)
sl2_training = matrix(NA, ncol = 7, nrow = 100)
int_testing = matrix(NA, ncol = 7, nrow = 100)
sl1_testing = matrix(NA, ncol = 7, nrow = 100)
sl2_testing = matrix(NA, ncol = 7, nrow = 100)

for (i in 1:100){
  training_data = oper[i, 8][[1]]
  testing_data = oper[i, 9][[1]]
  syn = oper[i, c(1:7)]
  
  training_fit = lm(x3 ~ x1 + x2, data = as.data.frame(training_data))
  tmp = sapply(syn, function(x) coef(lm(x3 ~ x1 + x2, data = as.data.frame(x))))
  tmp_dif = abs(tmp - matrix(rep(coef(training_fit), 7), nrow = 3))/sqrt(diag(vcov(training_fit)))
  int_training[i, ] = tmp_dif[1, ]
  sl1_training[i, ] = tmp_dif[2, ]
  sl2_training[i, ] = tmp_dif[3, ]
  
  holdout_fit = lm(x3 ~ x1 + x2, data = as.data.frame(testing_data))
  tmp_dif = abs(tmp - matrix(rep(coef(holdout_fit), 7), nrow = 3))/sqrt(diag(vcov(holdout_fit)))
  int_testing[i, ] = tmp_dif[1, ]
  sl1_testing[i, ] = tmp_dif[2, ]
  sl2_testing[i, ] = tmp_dif[3, ]
}


tab_coef = rbind(apply(int_training, 2, mean), apply(sl1_training, 2, mean), 
                 apply(sl2_training, 2, mean), apply(int_testing, 2, mean), 
                 apply(sl1_testing, 2, mean), apply(sl2_testing, 2, mean))

colnames(tab_coef) = methods[-8]
rownames(tab_coef) = c("Training Intercept", "Training Slope x1", "Training Slope x2",
                       "Testing Intercept", "Testing Slope x1", "Testing Slope x2")

tab_coef_se = rbind(apply(int_training, 2, sd), apply(sl_training, 2, sd),
                    apply(int_testing, 2, sd), apply(sl_testing, 2, sd)) / sqrt(100)

colnames(tab_coef_se) = methods[-8]
rownames(tab_coef_se) = paste(c("Training Intercept", "Training Slope", "Testing Intercept",
                                "Testing Slope"), "SE")

tab_coef_x3 = rbind(tab_coef, tab_coef_se)
tab_coef_x3

###############################################################################
## RMSE
std_rmse = matrix(NA, ncol = 7, nrow = 100)
for (i in 1: 100){
  cur = oper[i, ]
  holdout_data = cur[[9]]
  
  tmp = lapply(cur[c(1:7)], function(x) lm(x2 ~ x1, data = as.data.frame(x)))
  tmp_pred = lapply(tmp, function(x) predict(x, newdata = as.data.frame(holdout_data)))
  tmp_rmse = sapply(tmp_pred, function(x) sqrt(mean((x - holdout_data[,2])^2))/sqrt(var(holdout_data[,2])))
  std_rmse[i, ] = tmp_rmse
}

tab_rmse = rbind(apply(std_rmse, 2, mean), apply(std_rmse, 2, sd)/sqrt(100))


std_rmse = matrix(NA, ncol = 7, nrow = 100)

for (i in 1: 100){
  cur = oper[i, ]
  holdout_data = cur[[9]]
  
  tmp = lapply(cur[c(1:7)], function(x) lm(x3 ~ x1 + x2, data = as.data.frame(x)))
  tmp_pred = lapply(tmp, function(x) predict(x, newdata = as.data.frame(holdout_data)))
  tmp_rmse = sapply(tmp_pred, function(x) sqrt(mean((x - holdout_data[,3])^2))/sd(holdout_data[,3]))
  std_rmse[i, ] = tmp_rmse
}

tab_rmse = rbind(tab_rmse, apply(std_rmse, 2, mean),apply(std_rmse, 2, sd)/sqrt(100))

colnames(tab_rmse) = methods[-8]
rownames(tab_rmse) = c("x2 ~ x1", "x2 ~ x1 SE", "x3 ~ x1 + x2", "x3 ~ x1 + x2 SE")

tab_rmse

###############################################################################
## save all simulation results in an Rdata file
save(tab_wrt, tab_ut, tab_ut_inter, tab_km, tab_coef_x2, tab_coef_x3, tab_rmse, 
     file = "./output/simulation_utility_results.Rdata")


###############################################################################
## plot data
library(ggplot2)
library(gridExtra)

plot_data = list()
for (i in 1:8){
  plot_data[[i]] = oper[96,][[i]][, 3]
}
m = c("KNG", "StepFixed", "StepVarying", "SWFixed", 
      "SWVarying", "NonPrivate", "pMSE", "RawData")
names(plot_data) = m
plot_output = list() 

j = 1
for (i in c(1:5, 7)){
  plotdata = reshape2::melt(plot_data[c(i, 6, 8)])
  plotdata$L1 = factor(plotdata$L1, levels = c(m[6], m[8], m[i]))
  plot_output[[j]] = ggplot(plotdata,aes(x=value, fill= L1)) +
    geom_density(alpha=0.51) + 
    labs(fill='Methods')
  j = j + 1
} 
do.call(grid.arrange, plot_output)


# plotdata$L1 = factor(plotdata$L1, levels = c("Raw Data", "Non-Private", "pMSE Mechanism", "KNG",
#                                              "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
#                                              "Sandwich-Fixed Slope", "Sandwich-Varying Slope"))
# filename = paste(paste0('seed', t), paste0('e', ep*2), syn_var, sep = '_')
# png(paste("./plot/simulations/pMSE/",paste0('eps', ep*2) , "/simulation_", filename, ".png", sep = ""),
#    width = 700, height = 400)
# 
# 
# print(ggplot(plotdata, aes(x=value, fill = L1)) + facet_wrap(~L1, ncol = 4) +
#         stat_density(geom = "area", bw = 5, alpha = 0.5, size = 1) +
#         scale_color_brewer(palette="Dark2") + theme_minimal() +
#         theme(legend.position='none') +
#         labs(fill='Methods')  + coord_cartesian(xlim=c(-50, 100)) +
#         ggtitle(paste("Density of variable", syn_var, "- Eps", round(ep, 2))))
# dev.off()

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












