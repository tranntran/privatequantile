# This code is used to calculate the utility of generated synthetic data based on Task et al (2021).

rm(list = ls())
# change working directory.  Note that this is not an absolute pathname. 
# It is the pathname after one logs into the appropriate project space on the IRE.
wd = "<project directory>/programs"
library(ggplot2)
library(dplyr)
library(utils)

set.seed(123)
tab_nist = matrix(NA, nrow = 4, ncol = 6)
sics = c(178, 239, 542, 829) # These are SIC codes. 
methods = c("KNG", "StepFixed", "StepVarying", "SWFixed", 
            "SWVarying", "NonPriv", "Truth")
for (j in 1:length(sics)){
  nist_score = NULL
  sic = sics[j]
  for (seed in 1:5){
    data = list()
    # randomly select 3 columns to calculate the utility score
    samp = sample(1:25, 3)
    for (i in 1:length(methods)){
      df = read.csv(paste0(wd, "/amh1/", sic, "/", sic, "/output/SIC", sic, "_seed", seed, "_", methods[i], ".csv"))
      df_yr = df[, samp]
      data[[i]] = df_yr  
    }
    
    data = lapply(data, na.omit)
    
    # create 5 bins to sort data into
    tab = merge(merge(c(0:5), c(0:5), by = NULL), c(0:5), by = NULL)
    
    training_data = data[[7]]
    bins = apply(training_data, 2, quantile)
    vars = colnames(training_data)
    colnames(tab) = vars
    
  
    for (i in 1:7){
      tmp = cbind(findInterval(x = unlist(data[[i]][,1]), vec = bins[,1], left.open = T),
                  findInterval(x = unlist(data[[i]][,2]), vec = bins[,2], left.open = T),
                  findInterval(x = unlist(data[[i]][,3]), vec = bins[,3], left.open = T))
      colnames(tmp) = vars
      tmp = as.data.frame(tmp)
      tmp_tab = tmp %>% group_by_at(vars) %>% tally()
      names(tmp_tab)[colnames(tmp_tab) == "n"] = methods[i]
      tmp_tab = as.data.frame(tmp_tab)
      tab = merge(tab, tmp_tab, by = vars, all.x = T)
    }
    tab[is.na(tab)] = 0
    colSums(tab)
    raw_score = rep(NA, 6)
    for (i in 1:6){
      var = methods[i]
      tmp = abs(tab[, var]-tab[,ncol(tab)])/nrow(training_data)
      raw_score[i] = sum(tmp)
    }
  nist_score = rbind(nist_score, (2-raw_score)/2*1000) # hard-coded values are due to the formula
  
  }
  
  tab_nist[j, ] = apply(nist_score, 2, mean)
}


methods = c("KNG", "StepF", "StepV", "SWF", "SWV", "NP")
colnames(tab_nist) = methods
rownames(tab_nist) = paste("SIC", sics)
tab_nist = as.data.frame(tab_nist)
write.csv(tab_nist, paste0(wd, "/amh1/kmarginals_results.csv"))
