rm(list = ls())
load("./output/useful_syn_data_eps1.Rdata")

utility_logit = function(data, inds){
  data = cbind(data, inds)
  mod = glm(inds ~ ., data = data, family = binomial(link = "logit"))
  preds = predict(mod, type = "response")
  score = sum((preds-0.5)^2)/nrow(data)
  return(score)
}
tab_final = NULL

tab_10td = matrix(NA, nrow = 7, ncol = 10)
tab_10td_h = matrix(NA, nrow = 7, ncol = 10)
for (k in 1:10){
  tab_1td = matrix(NA, nrow = 7, ncol = 10)
  tab_1td_h = matrix(NA, nrow = 7, ncol = 10)
  for (j in 1:10){
    cur = oper[k, j][[1]]
    inds = as.factor(c(rep(0, n), rep(1, n)))
    training_data = cur[[8]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      ut_data = lapply(syn, rbind, training_data)
      ut_data = lapply(ut_data, as.data.frame)
      tab_1td[i, j] = mean(sapply(ut_data, utility_logit, inds = inds))
      
      ut_data_h = lapply(syn, rbind, holdout_data)
      ut_data_h = lapply(ut_data_h, as.data.frame)
      tab_1td_h[i, j] = mean(sapply(ut_data_h, utility_logit, inds = inds))
    }
  }
  tab_10td[,k] = apply(tab_1td, 1, mean)
  tab_10td_h[,k] = apply(tab_1td_h, 1, mean)
}

tab_final = rbind(tab_final, apply(tab_10td, 1, mean)) #training pmse
tab_final = rbind(tab_final, apply(tab_10td_h, 1, mean))


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


w_tab = matrix(NA, ncol = 8, nrow = 25)
for (j in 1: 25){
  test = oper[j,]
  for (i in 1:8){
    res = wasserstein_randomization_test(test[[i]], test[[8]], n_rep = 10000)
    null = median(res$dist)
    ratio = res$sample_distance/null
    w_tab[j, i] = ratio
  }
}

set.seed(123)
tab_10td = NULL
tab_10td_h = NULL
for (k in 1:10){
  tab_1td = matrix(NA, nrow = 7, ncol = 10)
  tab_1td_h = matrix(NA, nrow = 7, ncol = 10)
  for (j in 1:10){
    cur = oper[k, j][[1]]
    training_data = cur[[8]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      tmp = lapply(syn, wasserstein_randomization_test, b = training_data, n_rep = 1000)
      tab_1td[i, j] = mean(sapply(tmp, function(x) x$sample_distance/median(x$dist)))
      
      
      tmp_h = lapply(syn, wasserstein_randomization_test, b = holdout_data, n_rep = 1000)
      tab_1td_h[i, j] = mean(sapply(tmp_h, function(x) x$sample_distance/median(x$dist)))
    }
  }
  tab_10td = cbind(tab_10td, tab_1td)
  tab_10td_h = cbind(tab_10td_h, tab_1td_h)
}

tab_final = rbind(tab_final, apply(tab_10td, 1, mean)) #training wrt
tab_final = rbind(tab_final, apply(tab_10td_h, 1, mean))

tab_final
colnames(tab_final) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                     "Sandwich Fixed Slope", "Sandwich Varying Slope",
                     "NonPrivate", "pMSE Mechanism")
rownames(tab_final) = c("Training pMSE", "Generalization pMSE", 
                        "Training WRT", "Generalization WRT")
save(tab_final, file = "./output/useful_syn_quality.Rdata")


# set.seed(123)
# tab_10td = NULL
# tab_10td_h = NULL
# for (k in 1:10){
#   tab_1td = NULL
#   tab_1td_h = NULL
#   for (j in 1:10){
#     cur = oper[k, j][[1]]
#     training_data = cur[[8]]
#     holdout_data = cur[[9]]
#     for (i in 1:7){
#       syn = cur[[i]]
#       tmp = lapply(syn, wasserstein_randomization_test, b = training_data, n_rep = 100)
#       tab_1td = cbind(tab_1td, sapply(tmp, function(x) x$sample_distance/median(x$dist)))
#       
#       
#       tmp_h = lapply(syn, wasserstein_randomization_test, b = holdout_data, n_rep = 100)
#       tab_1td_h = cbind(tab_1td_h, sapply(tmp_h, function(x) x$sample_distance/median(x$dist)))
#     }
#   }
#   tab_10td = cbind(tab_10td, tab_1td)
#   tab_10td_h = cbind(tab_10td_h, tab_1td_h)
# }


tab_10td = rep(list(NULL), 7)
tab_10td_h = rep(list(NULL), 7)
for (k in 1:10){
  tab_1td = rep(list(NULL), 7)
  tab_1td_h = rep(list(NULL), 7)
  for (j in 1:10){
    cur = oper[k, j][[1]]
    training_data = cur[[8]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      training_coef = coef(lm(x2 ~ x1, data = as.data.frame(training_data)))
      tmp = sapply(syn, function(x) coef(lm(x2 ~ x1, data = as.data.frame(x))))
      tmp_dif = 100*abs(rowMeans(tmp - matrix(rep(training_coef, 10), nrow = 2))/training_coef)
      tab_1td[[i]] = cbind(tab_1td[[i]], tmp_dif)
      
      holdout_coef = coef(lm(x2 ~ x1, data = as.data.frame(holdout_data)))
      #holdout_coef = c(4, 3) #true values defining the data generating process
      tmp_dif = 100*abs(rowMeans(tmp - matrix(rep(holdout_coef, 10), nrow = 2))/holdout_coef)
      tab_1td_h[[i]] = cbind(tab_1td_h[[i]], tmp_dif)
      

    }
  }
  
  tmp = lapply(tab_1td, rowMeans)
  tmp_h = lapply(tab_1td_h, rowMeans)
  for (i in 1:7){
    tab_10td[[i]] = cbind(tab_10td[[i]], tmp[[i]])
    tab_10td_h[[i]] = cbind(tab_10td_h[[i]], tmp_h[[i]])
  }
  
}

tab_final = rbind(tab_final, sapply(lapply(tab_10td, rowMeans), mean))
tab_final = rbind(tab_final, sapply(lapply(tab_10td_h, rowMeans), mean))


rownames(tab_final) = c("Training pMSE", "Generalization pMSE", 
                        "Training WRT", "Generalization WRT", 
                        "Training Coef Bias %", "Generalization Coef Bias %")
save(tab_final, file = "./output/useful_syn_quality.Rdata")



tab_10td = rep(list(NULL), 7)
for (k in 1:10){
  tab_1td = rep(list(NULL), 7)
  for (j in 1:10){
    cur = oper[k, j][[1]]
    training_data = cur[[8]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      training_vcov = vcov(lm(x2 ~ x1, data = as.data.frame(training_data)))
      tmp = lapply(syn, function(x) vcov(lm(x2 ~ x1, data = as.data.frame(x))))
      tmp_r = lapply(tmp, function(x) abs((x - training_vcov) / training_vcov))
      tmp_r = sapply(tmp_r, function(x) x[lower.tri(x, diag = T)])
      tab_1td[[i]] = cbind(tab_1td[[i]], rowMeans(tmp_r))
    }
  }
  
  tmp = lapply(tab_1td, rowMeans)
  for (i in 1:7){
    tab_10td[[i]] = cbind(tab_10td[[i]], tmp[[i]])
  }
  
}

tab_final = rbind(tab_final, sapply(lapply(tab_10td, rowMeans), mean))

rownames(tab_final) = c("Training pMSE", "Generalization pMSE", 
                        "Training WRT", "Generalization WRT", 
                        "Training Coef Bias %", "Generalization Coef Bias %",
                        "Training Cov Ratio")
save(tab_final, file = "./output/useful_syn_quality.Rdata")



tab_10td = rep(list(NULL), 7)
tab_10td_h = rep(list(NULL), 7)
for (k in 1:10){
  tab_1td = rep(list(NULL), 7)
  tab_1td_h = rep(list(NULL), 7)
  for (j in 1:10){
    cur = oper[k, j][[1]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      true_coef = c(4, 3)
      holdout_coef = coef(lm(x2 ~ x1, data = as.data.frame(holdout_data)))
      tmp = lapply(syn, function(x) confint(lm(x2 ~ x1, data = as.data.frame(x)), level = 0.9))
      tmp_p = sapply(tmp, function(x) x[,1] < holdout_coef & holdout_coef < x[,2])
      tab_1td[[i]] = cbind(tab_1td[[i]], rowMeans(tmp_p))
      
      tmp_p = sapply(tmp, function(x) x[,1] < true_coef & true_coef < x[,2])
      tab_1td_h[[i]] = cbind(tab_1td_h[[i]], rowMeans(tmp_p))
      
    }
  }
  
  tmp = lapply(tab_1td, rowMeans)
  tmp_h = lapply(tab_1td_h, rowMeans)
  for (i in 1:7){
    tab_10td[[i]] = cbind(tab_10td[[i]], tmp[[i]])
    tab_10td_h[[i]] = cbind(tab_10td_h[[i]], tmp_h[[i]])
  }
  
}
tab_tmp = rbind(sapply(lapply(tab_10td, rowMeans), mean), sapply(lapply(tab_10td_h, rowMeans), mean))
rownames(tab_tmp) = c("Gen. Coverage Rate Holdout Data", "Gen. Coverage Rate True Coef" )

tab_final = rbind(tab_final, tab_tmp)


tab_10td = matrix(NA, nrow = 7, ncol = 10)
for (k in 1:10){
  tab_1td = matrix(NA, nrow = 7, ncol = 10)
  for (j in 1:10){
    cur = oper[k, j][[1]]
    training_data = cur[[8]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      tmp = lapply(syn, function(x) lm(x2 ~ x1, data = as.data.frame(x)))
      tmp_pred = lapply(tmp, function(x) predict(x, newdata = as.data.frame(holdout_data)))
      tmp_rmse = sapply(tmp_pred, function(x) sqrt(mean((x - holdout_data[,2])^2))/sqrt(var(holdout_data[,2])))
      tab_1td[i, j] = mean(tmp_rmse)
    }
  }
  
  tab_10td[,k] = apply(tab_1td, 1, mean)
  
}
tab_final = rbind(tab_final, rowMeans(tab_10td))
row.names(tab_final)[10] = "Gen. Prediction RMSE"


row.names(tab_final)[11] = "Gen. Coef Bias % Data"
save(tab_final, file = "./output/useful_syn_quality.Rdata")

round(tab_final, 4)


tab_10td = rep(list(NULL), 7)
tab_10td_h = rep(list(NULL), 7)
for (k in 1:10){
  for (j in 1:10){
    cur = oper[k, j][[1]]
    training_data = cur[[8]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      training_fit = lm(x2 ~ x1, data = as.data.frame(training_data))
      tmp = sapply(syn, function(x) coef(lm(x2 ~ x1, data = as.data.frame(x))))
      tmp_dif = abs(tmp - matrix(rep(coef(training_fit), 10), nrow = 2))/sqrt(diag(vcov(training_fit)))
      tab_10td[[i]] = cbind(tab_10td[[i]], tmp_dif)
      
      holdout_fit = lm(x2 ~ x1, data = as.data.frame(holdout_data))
      tmp_dif = abs(tmp - matrix(rep(coef(holdout_fit), 10), nrow = 2))/sqrt(diag(vcov(holdout_fit)))
      tab_10td_h[[i]] = cbind(tab_10td_h[[i]], tmp_dif)
      
      
    }
  }
  
}


tab1 = sapply(tab_10td, rowMeans)
colnames(tab1) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                                         "Sandwich Fixed Slope", "Sandwich Varying Slope",
                                         "NonPrivate", "pMSE Mechanism")

tab1_h = sapply(tab_10td_h, rowMeans)
colnames(tab1_h) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                   "Sandwich Fixed Slope", "Sandwich Varying Slope",
                   "NonPrivate", "pMSE Mechanism")


calculate_io = function(syn_int, orig_int){
  io = rep(NA, 2)
  for (i in 1:2){
    val = min(orig_int[i, 2], syn_int[i, 2]) - max(orig_int[i, 1], syn_int[i, 1])
    io[i] = (val/(orig_int[i, 2] - orig_int[i, 1]) + val/(syn_int[i, 2] - orig_int[i, 1]))/2
  }
  return(io)
}

tab_10td = rep(list(NULL), 7)
tab_10td_h = rep(list(NULL), 7)
for (k in 1:10){
  for (j in 1:10){
    cur = oper[k, j][[1]]
    training_data = cur[[8]]
    holdout_data = cur[[9]]
    for (i in 1:7){
      syn = cur[[i]]
      training_fit = confint(lm(x2 ~ x1, data = as.data.frame(training_data)), level = 0.95)
      tmp = lapply(syn, function(x) confint(lm(x2 ~ x1, data = as.data.frame(x)), level = 0.95))
      tmp_dif = sapply(tmp, calculate_io, orig_int = training_fit)
      tab_10td[[i]] = cbind(tab_10td[[i]], tmp_dif)
      
      holdout_fit = confint(lm(x2 ~ x1, data = as.data.frame(holdout_data)), level = 0.95)
      tmp_dif = sapply(tmp, calculate_io, orig_int = holdout_fit)
      tab_10td_h[[i]] = cbind(tab_10td_h[[i]], tmp_dif)
      
      
    }
  }
  
}

tab2 = sapply(tab_10td, rowMeans)
rownames(tab2) = c("Intercept", "Slope")
colnames(tab2) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                   "Sandwich Fixed Slope", "Sandwich Varying Slope",
                   "NonPrivate", "pMSE Mechanism")

tab2 = sapply(tab_10td_h, rowMeans)
rownames(tab2) = c("Intercept", "Slope")
colnames(tab2) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                   "Sandwich Fixed Slope", "Sandwich Varying Slope",
                   "NonPrivate", "pMSE Mechanism")

