library(doParallel)
library(tictoc)

num_cores=detectCores()-1 #use all available core but 1

workers=makeCluster(num_cores,type="SOCK",outfile="log.txt")
registerDoParallel(workers)


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


calculate_utility = function(truth, data, inds){
  ut.data = lapply(data, rbind, truth)
  ut_logit= sapply(ut.data, utility_logit, inds)
  ut_logit_inter = t(sapply(ut.data, utility_logit_inter, inds))
  return(list(ut_logit, ut_logit_inter))
}

tic()
Out = foreach(seed = c(1:5),.combine=rbind, .errorhandling='stop') %dopar% {
  library(data.table)
  set.seed(seed)
  inds = c(rep(1, 10), rep(0, 10))
  truth = data.table(sample(1:10, 30, replace = T), ncol = 3)
  data = list( data.table(sample(1:50, 30, replace = T), ncol = 3),
               data.table(sample(1:20, 30, replace = T), ncol = 3),
               truth)
  #calculate_utility(truth, data, inds)
  ut.data = lapply(data, rbind, truth)
  sapply(ut.data, utility_logit, inds)
  
}
filename = paste("SIC178_utility_logitstic.Rdata")

save(list = c("Out"), file = filename)

Out1 = foreach(seed = c(1:5),.combine=rbind, .errorhandling='stop') %dopar% {
  library(data.table)
  set.seed(seed)
  inds = c(rep(1, 10), rep(0, 10))
  truth = data.table(sample(1:10, 30, replace = T), ncol = 3)
  data = list( data.table(sample(1:50, 30, replace = T), ncol = 3),
               data.table(sample(1:20, 30, replace = T), ncol = 3),
               truth)
  #calculate_utility(truth, data, inds)
  ut.data = lapply(data, rbind, truth)
  t(sapply(ut.data, utility_logit_inter, inds))
  
}

colnames(Out1) = paste0("var", 1:3)

filename = paste("SIC178_utility_logitstic_inter.csv")

fwrite(Out1, filename)
#save(list(Out))
stopCluster(workers)

Out2 = fread("SIC178_utility_logitstic_inter.csv")
