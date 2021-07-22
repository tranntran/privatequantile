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

wf_tab = apply(w_tab, 2, mean)
wf_tab = rbind(wf_tab, apply(w_tab, 2, sd)/5)
wf_tab = wf_tab[,-8]
colnames(wf_tab) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                     "Sandwich Fixed Slope", "Sandwich Varying Slope",
                     "NonPrivate", "pMSE Mechanism")
rownames(wf_tab) = c("WD Ratio Mean", "WD Ratio SE")
wf_tab = round(t(wf_tab), 5)
library(gridExtra)
png(paste("./plot/WD_pMSE_seed1to25_eps1.png", sep =""), height = 35*nrow(wf_tab),
    width = 200*ncol(wf_tab))
grid.table(wf_tab)
dev.off()



rm(list = ls())
load("./output/pMSE/data_eps1.Rdata")

dif_int = matrix(NA, ncol = 7, nrow = 25)
dif_sl = matrix(NA, ncol = 7, nrow = 25)

for (j in 1:25){
  test = oper[j, ]
  truth = coef(lm(test[[8]][,2] ~ test[[8]][,1]))
  for (i in 1:7){
    syn = coef(lm(test[[i]][,2] ~ test[[i]][,1]))
    dif = abs(syn-truth)
    dif_int[j, i] = dif[1]
    dif_sl[j, i] = dif[2]
  }
}

tab_ana = rbind(apply(dif_int, 2, mean), apply(dif_int, 2, sd)/5, 
                apply(dif_sl, 2, mean), apply(dif_sl, 2, sd)/5)
tab_ana
colnames(tab_ana) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                      "Sandwich Fixed Slope", "Sandwich Varying Slope",
                      "NonPrivate", "pMSE Mechanism")
rownames(tab_ana) = c("Mean Difference Intercept", "SE Difference Intercept", 
                      "Mean Difference Slope", "SE Difference Slope")
tab_anaf = t(round(tab_ana, 5))

library(gridExtra)
png(paste("./plot/analysis_pMSE_seed1to25_eps1.png", sep =""), height = 35*nrow(tab_anaf),
    width = 200*ncol(tab_anaf))
grid.table(tab_anaf)
dev.off()


rm(list = ls())
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

utility <- function(data,inds){
  mod=rpart::rpart(inds~., data=data, method="class",minbucket=5,cp=1e-4)
  preds=predict(mod)[,2]
  score=sum((preds-0.5)^2)/nrow(data)
  return(score)
}


eps = c("01", "025", "05", "1")
#eps = c("025", "05", "1")
tab = NULL
for (i in 1:length(eps)){
  load(paste0("./output/pMSE/data_eps", eps[i], ".Rdata"))
  
  n = 5000
  inds = as.factor(c(rep(0,n), rep(1,n)))
  tab_logit = NULL
  #tab_logit_inter = NULL
  for (i in 1:25){
    df = oper[i,]
    data = df[[8]]
    ut.data = lapply(oper[2,], rbind, data)
    ut.data = lapply(ut.data, as.data.frame)
    tab_logit = cbind(tab_logit, sapply(ut.data, utility_logit, inds))
    #tab_logit_inter = cbind(tab_logit_inter, sapply(ut.data, utility_logit_inter, inds))
  }
  
  tab = rbind(tab, apply(tab_logit, 1, mean))
  
}

tab = tab[,-8]
rownames(tab) = c(0.1, 0.25, 0.5, 1)
colnames(tab) = c("KNG", "Stepwise Fixed Slope", "Stepwise Varying Slope",
                  "Sandwich Fixed Slope", "Sandwich Varying Slope",
                  "NonPrivate", "pMSE Mechanism")
tab = 1/tab
tab = tab[,-6]
tab_melt = reshape2::melt(tab)
tab_melt$Var2 = as.factor(tab_melt$Var2)
tab_melt$Eps <- ordered(tab_melt$Var1, levels=c(0.1, 0.25, 0.5, 1))
ggplot(as.data.frame(tab_melt), aes(x=Eps, y=value, col = Var2, group = Var2)) +
  geom_line() +
  geom_point() + 
  labs(x = "Epsilon", y = "1/pMSE", col = "Methods")


# library(ggplot2)
# ggplot(data = as.data.frame(tab_melt), aes(x = Var1, y = value, col = Var2)) +
#   geom_point() + geom_line()
# 
# 
# ggplot(data = as.data.frame(tab_melt1), aes(x = Var1, y = value, c, group = 1)) +
#   geom_point() + geom_line() +
#   scale_x_continuous(breaks = 1:3, labels=c("0.25", '0.5', '1')) + 
#   labs(x = "Epsilon", y = "1/pMSE", col = "Methods")

tab1 = tab[, -c(1, 6, 7)]
#tab1 = 1/tab1
tab_melt1 = reshape2::melt(tab1)
tab_melt1$Var2 = as.factor(tab_melt1$Var2)
tab_melt1$Eps <- ordered(tab_melt1$Var1, levels=c(0.1, 0.25, 0.5, 1))
ggplot(as.data.frame(tab_melt1), aes(x=Eps, y=value, col = Var2, group = Var2)) +
  geom_line() +
  geom_point() + 
  labs(x = "Epsilon", y = "1/pMSE", col = "Methods")



ydensity <- ggplot(df, aes(y, fill=group)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none")
ydensity


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