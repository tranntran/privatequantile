# Use this code with the output from compare_pMSE_data.R
rm(list = ls())
load("./output/data_eps1.Rdata")

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
  training_data = oper[i, 8][[1]]
  testing_data = oper[i, 9][[1]]
  tmp = lapply(oper[i, c(1:7)], wasserstein_randomization_test, b = training_data)
  wrt_training[i, ] = sapply(tmp, function(x) x$sample_distance/median(x$dist))
  
  tmp = lapply(oper[i, c(1:7)], wasserstein_randomization_test, b = testing_data)
  wrt_testing[i, ] = sapply(tmp, function(x) x$sample_distance/median(x$dist))
  
}
apply(wrt_testing, 2, mean)
apply(wrt_testing, 2, median)
apply(wrt_training, 2, mean)
apply(wrt_training, 2, median)


wrt_training_1 = matrix(NA, ncol = 7, nrow = 100)
wrt_testing_1 = matrix(NA, ncol = 7, nrow = 100)
for (k in 1:100){
  training_data = oper[k, 8][[1]]
  testing_data = oper[k, 9][[1]]
  for (j in 1:7){
    test = oper[1, j][[1]]
    ans = rep(NA, 3)
    for (i in 1:3){
      tmp_wrt = wasserstein_randomization_test(a = test[,i], b = training_data[,i])
      ans[i] = tmp_wrt$sample_distance/median(tmp_wrt$dist)
    }
    wrt_training_1[k, j] = mean(ans)
    
    ans = rep(NA, 3)
    for (i in 1:3){
      tmp_wrt = wasserstein_randomization_test(a = test[,i], b = testing_data[,i])
      ans[i] = tmp_wrt$sample_distance/median(tmp_wrt$dist)
    }
    wrt_testing_1[k, j] = mean(ans)
  }
}


apply(wrt_testing_1, 2, mean)
apply(wrt_testing_1, 2, median)
apply(wrt_training_1, 2, mean)
apply(wrt_training_1, 2, median)

wrt_tab = rbind(apply(wrt_training, 2, mean), apply(wrt_testing, 2, mean))
colnames(wrt_tab) = methods[-8]
rownames(wrt_tab) = c("Training", "Testing")

plot_data = NULL
for (i in 1:8){
  plot_data = cbind(plot_data, oper[1, i][[1]][, 2])
}
colnames(plot_data) = c("KNG", "StepFixed", "StepVarying", "SWFixed", 
                        "SWVarying", "NonPriv", "pMSE", "RawData")

plot_data_m = melt(plot_data[, -1])
ggplot(plot_data_m,aes(x=value, fill= Var2)) + geom_density(alpha=0.25)

ggplot(plot_data_m,aes(x=value, fill= Var2)) + geom_density(alpha=0.25) +
  coord_cartesian(xlim = c(200, 1000), ylim = c(0, 0.0025))


plot_data1 = plot_data[, c(8, 5, 6)]
plot_data_m = melt(plot_data1)
ggplot(plot_data_m,aes(x=value, fill= Var2)) + geom_density(alpha=0.25)
ggplot(plot_data_m,aes(x=value, fill= Var2)) + geom_density(alpha=0.25) +
  coord_cartesian(xlim = c(100, 250), ylim = c(0, 0.005))

training_data = oper[1, 8][[1]]
testing_data = oper[1, 9][[1]]
ans = rep(NA, 7)
for (i in 1:7){
  tmp = wasserstein_randomization_test(oper[1, i][[1]][,3], training_data[,3])
  ans[i] = tmp$sample_distance/median(tmp$dist)
  
}
ansf = rbind(ansf, ans)
colnames(ansf) = c("KNG", "StepFixed", "StepVarying", "SWFixed", 
                        "SWVarying", "NonPriv","pMSE")
ansf = rbind(ansf, apply(ansf, 2, mean))
rownames(ansf) = c("x1", "x2", "x3", "Avg")


sum(abs(sort(training_data[,2]) - sort(oper[1, ][[1]][,2]))/2)


# Total variation distance
ans = matrix(NA, ncol = 7, nrow = 2)
for (i in 1:7){
  ans[1, i] = sum(abs(sort(c(oper[1, i][[1]])) - sort(c(oper[1, 8][[1]])))/2)
  ans[2, i] = (sum(abs(sort(c(oper[1, i][[1]][,1])) - sort(c(oper[1, 8][[1]][,1])))/2) + 
                 sum(abs(sort(c(oper[1, i][[1]][,2])) - sort(c(oper[1, 8][[1]][,2])))/2) +
                 sum(abs(sort(c(oper[1, i][[1]][,3])) - sort(c(oper[1, 8][[1]][,3])))/2))/3

}












###############################################################################
# Draft implementation of k-marginals
library(data.table)
tab = CJ(x1=0:5,x2=0:5, x3=0:5)

training_data = oper[1, 8][[1]]
x1_check = quantile(training_data[,1])
x2_check = quantile(training_data[,2])
x3_check = quantile(training_data[,3])

methods = c("KNG", "StepF", "StepV", "SWF", "SWV", "NP","pMSE", "Raw")
for (i in 1:8){
  # tmp = cbind(cut(oper[1, i][[1]][,1], breaks=x1_check, labels = F, include.lowest = T),
  #             cut(oper[1, i][[1]][,2], breaks=x2_check, labels = F, include.lowest = T),
  #             cut(oper[1, i][[1]][,3], breaks=x3_check, labels = F, include.lowest = T))
  tmp = cbind(findInterval(x = oper[1, i][[1]][,1], vec = x1_check, left.open = T),
              findInterval(x = oper[1, i][[1]][,2], vec = x2_check, left.open = T),
              findInterval(x = oper[1, i][[1]][,3], vec = x3_check, left.open = T))
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
nist_score = (2-raw_score)/2*1000
names(nist_score) = methods[-8]
