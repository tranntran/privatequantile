# Use this code with the output from compare_pMSE_data.R
rm(list = ls())
load("./output/data_eps1_50q_unif.Rdata")

methods = c("KNG", "StepF", "StepV", "SWF", "SWV", "NP","pMSE", "Raw")
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

###############################################################################
## utility

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

ut = rbind(apply(utility_training, 2, mean), apply(utility_testing, 2, mean))
colnames(ut) = methods[-8]
rownames(ut) = c("Training", "Testing")
ut

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

tab_rmse = apply(std_rmse, 2, mean)


std_rmse = matrix(NA, ncol = 7, nrow = 100)

for (i in 1: 100){
  cur = oper[i, ]
  holdout_data = cur[[9]]
  
  tmp = lapply(cur[c(1:7)], function(x) lm(x3 ~ x1 + x2, data = as.data.frame(x)))
  tmp_pred = lapply(tmp, function(x) predict(x, newdata = as.data.frame(holdout_data)))
  tmp_rmse = sapply(tmp_pred, function(x) sqrt(mean((x - holdout_data[,3])^2))/sd(holdout_data[,3]))
  std_rmse[i, ] = tmp_rmse
}

tab_rmse = rbind(tab_rmse, apply(std_rmse, 2, mean))

colnames(tab_rmse) = methods[-8]
rownames(tab_rmse) = c("x2 ~ x1", "x3 ~ x1 + x2")

tab_rmse




###############################################################################
## plot data
library(ggplot2)
plot_data = list()
for (i in 1:8){
  plot_data[[i]] = oper[96,][[i]][, 2]
}
names(plot_data) = c("KNG", "StepFixed", "StepVarying", "SWFixed", 
                        "SWVarying", "NonPriv", "pMSE", "RawData")

for (i in 2:5){
  plotdata = reshape2::melt(plot_data[c(i, 6, 8)])
  print(ggplot(plotdata,aes(x=value, fill= L1)) + geom_density(alpha=0.51))
} 








# wrt_training_1 = matrix(NA, ncol = 7, nrow = 100)
# wrt_testing_1 = matrix(NA, ncol = 7, nrow = 100)
# for (k in 1:100){
#   training_data = oper[k, 8][[1]]
#   testing_data = oper[k, 9][[1]]
#   for (j in 1:7){
#     test = oper[1, j][[1]]
#     ans = rep(NA, 3)
#     for (i in 1:3){
#       tmp_wrt = wasserstein_randomization_test(a = test[,i], b = training_data[,i])
#       ans[i] = tmp_wrt$sample_distance/median(tmp_wrt$dist)
#     }
#     wrt_training_1[k, j] = mean(ans)
#     
#     ans = rep(NA, 3)
#     for (i in 1:3){
#       tmp_wrt = wasserstein_randomization_test(a = test[,i], b = testing_data[,i])
#       ans[i] = tmp_wrt$sample_distance/median(tmp_wrt$dist)
#     }
#     wrt_testing_1[k, j] = mean(ans)
#   }
# }
# 
# 
# apply(wrt_testing_1, 2, mean)
# apply(wrt_testing_1, 2, median)
# apply(wrt_training_1, 2, mean)
# apply(wrt_training_1, 2, median)

wrt_tab = rbind(apply(wrt_training, 2, mean), apply(wrt_testing, 2, mean))
colnames(wrt_tab) = methods[-8]
rownames(wrt_tab) = c("Training", "Testing")





# training_data = oper[1, 8][[1]]
# testing_data = oper[1, 9][[1]]
# ans = rep(NA, 7)
# for (i in 1:7){
#   tmp = wasserstein_randomization_test(oper[1, i][[1]][,3], training_data[,3])
#   ans[i] = tmp$sample_distance/median(tmp$dist)
#   
# }
# ansf = rbind(ansf, ans)
# colnames(ansf) = c("KNG", "StepFixed", "StepVarying", "SWFixed", 
#                         "SWVarying", "NonPriv","pMSE")
# ansf = rbind(ansf, apply(ansf, 2, mean))
# rownames(ansf) = c("x1", "x2", "x3", "Avg")
# 
# 
# sum(abs(sort(training_data[,2]) - sort(oper[1, ][[1]][,2]))/2)
# 
# 
# # Total variation distance
# ans = matrix(NA, ncol = 7, nrow = 2)
# for (i in 1:7){
#   ans[1, i] = sum(abs(sort(c(oper[1, i][[1]])) - sort(c(oper[1, 8][[1]])))/2)
#   ans[2, i] = (sum(abs(sort(c(oper[1, i][[1]][,1])) - sort(c(oper[1, 8][[1]][,1])))/2) + 
#                  sum(abs(sort(c(oper[1, i][[1]][,2])) - sort(c(oper[1, 8][[1]][,2])))/2) +
#                  sum(abs(sort(c(oper[1, i][[1]][,3])) - sort(c(oper[1, 8][[1]][,3])))/2))/3
# 
# }













