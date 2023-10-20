# This code generate synthetic datasets under SIC 178 using KNG, 
# stepwise fixed slope KNG,stepwise varying slope KNG, sandwich 
# fixed slope KNG, and sandwich varying slope KNG. It synthesizes
# variable emp from year 1976 to year 2000 sequentially. It can
# be used to synthesize both variable emp and pay over the 25 years
# if endcol (in line xx) is increased to 54.  

# Input: Real data (longitudinalSynLBD.csv)
# Output: 5 different synthetic datasets generated under 5 different seeds for
# each method, save under csv form.

# This function generates synthetic data using KNG
# Note: The hard-coded numbers in this function is the default parameters for MCMC.
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

# This function subsets the full dataset to only include variables
# needed for the synthesis procedure. The variable being synthesized
# is set to be in the last column. This modeling procedure is adapted
# from Kinney et al (2011) and Pistner et al (2018).
# For ex: 
# emp_1976 is generated with private quantiles (emp_1976 ~ 1), 
# so the data should only include 1 column of emp_1976.
# emp_1977 is generated using emp_1976, so the resulting dataset will
# include emp_1976 and emp_1977.
# pay_1976 is generated using emp_1976, so the data will have one column
# of emp_1976 and another column of pay_1976.
# pay_1977 is generated using emp_1976 and pay_1976, so the first column 
# will be emp_1976, then pay_1976 and pay_1977.
# This model assumes that firstyear, lastyear, and mu to be publicly 
# available information and we do not use or synthesize them.
select_var = function(syn_var){
  check = substr(syn_var, 1, 3)
  year = as.numeric(substr(syn_var, 5,8))
  
  if (check == "emp"){
    if (year == 1976){
      return(syn_var)
    } else {
      #ans = c("firstyear", "lastyear", "mu", paste("emp_", year-1, sep = ""))
      ans = paste("emp_", year-1, sep = "")
    }
  } else if (check == "pay"){
    if (year == 1976){
      ans = c("emp_1976")
    } else {
      ans = c(paste("emp_", year, sep = ""), paste("pay_", year-1, sep = ""))
    }
    
  }
  
  return(c(ans, syn_var))
}


# function for subsetting establishments based on  the year of synthesizing
# variable. For example, if we want to synthesize variable emp_1978,
# then the year is 1978. Any establishments founded before this year and closed
# after this year are put into df_model. Any establishments found on this year
# are put into df_birth. Finally, any establishments not existed during this year
# is put into df_na. It is because df_model and df_birth are generated differently
# based on Kinney et al, 2011. emp_1978 for any establishments in df_na is filled
# with NA.
subset_data = function(data, syn_var){
  check = substr(syn_var, 1, 3)
  if(check == "emp" | check == "pay"){
    year=substr(syn_var,5,8)
    
    df_model=data[which(data$firstyear < year & data$lastyear > year),]
    df_na= data[which(data$firstyear > year | data$lastyear <= year),]
    df_birth=data[which(data$firstyear== year),]
    
  }
  return(list(df_model, df_birth, df_na))
}

# This function is used to create synthetic data for a variable based
# on all possible predicted values on all quantiles and a vector of randomly
# selected quantiles. For example, if Y is model with 1 variable, then X is 
# an nx2 matrix, with the first column of 1, and beta_result is a 2xp matrix, 
# where p is the number of quantiles used to generate synthetic data. X*B will
# result in a nxp matrix. For each of the row from 1 to n, we choose a quantile 
# to be the synthetic value, based on the previously sample vector select_quantile.
# The function will return a nx1 matrix.
syndata = function(beta_result, x, select_quantile){
  allsyn = x%*%beta_result
  coord = cbind(c(1:nrow(x)), select_quantile)
  ans = allsyn[coord]
  return(ans)
}

# function to jitter the data
# to avoid singular design issue when fitting quantile regression
data_jitter <- function(data){
  data = as.data.frame(data)
  for(i in 1:ncol(data)){
    names=names(data)
    tmp=data[,names[i]]
    tmp=tmp + runif(length(tmp),-1,1)
    data[,names[i]]=tmp
  }
  return(data)
}


# This function generate synthetic data using all methods and saved them as csv files.
# totaleps is the total privacy budget used for the whole process.
# If endcol = 29, then we synthesize emp from 1976 to 2000.
# If endcol = 54, then we synthesize emp and pay from 1976 to 2000.

# Notes: 
# tau and main_tau are the quantiles to synthesize. 
# runs is the number of steps in MCMC.
# totaleps is the total privacy-loss budget.
# startcol, endcol, and start indicate which variables to synthesize.

syn_all_methods = function(sic, startcol = 2, endcol = 54, start = 4, totaleps = 1, 
                           runs = 1000, main_tau = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99),
                           tau = c(seq(0.05, 0.95, 0.05), 0.99)) { 
  data = read.csv(paste0(wd, "/LBD_178.csv"), na.strings="NA") #change this based on SIC
  data = data[, c(startcol: endcol)]
  
  nsyn = ncol(data) - start + 1
  vars = names(data)
  syn_final = rep(list(data[, c(1:(start-1))]), 6)
  ep = totaleps/nsyn #need to account for ep*2*25 
  
  # scale (mcmc step sizes) for variable emp
  # Notes: These are the default step sizes based on simulation study
  scale_emp_fy = list(0.05, 0.05, 0.05, c(0.1, 0.05), c(0.1, 0.05))
  scale_emp_oy = list(5e-4, 5e-4, 5e-4, c(1e-3, 5e-4), c(1e-3, 5e-4))
  
  for (j in start:ncol(data)){
    syn_var = vars[j]
    print(syn_var)
    
    df_subset = subset_data(data, syn_var)
    df_mod = data_jitter(df_subset[[1]]) 
    df_birth = data_jitter(df_subset[[2]])
    df_na = df_subset[[3]]
    
    original = subset_data(syn_final[[1]], syn_var)
    step_k = subset_data(syn_final[[2]], syn_var)
    step_d = subset_data(syn_final[[3]], syn_var)
    sw_k = subset_data(syn_final[[4]], syn_var)
    sw_d = subset_data(syn_final[[5]], syn_var)
    syn_rq = subset_data(syn_final[[6]], syn_var)
    
    syn_mod = list(original[[1]], step_k[[1]], step_d[[1]], sw_k[[1]], sw_d[[1]], syn_rq[[1]])
    syn_birth = list(original[[2]], step_k[[2]], step_d[[2]], sw_k[[2]], sw_d[[2]], syn_rq[[2]])
    syn_na = list(original[[3]], step_k[[3]], step_d[[3]], sw_k[[3]], sw_d[[3]], syn_rq[[3]])
    
    pred_var_mod = select_var(syn_var)
    pred_var_mod_ns = pred_var_mod[-length(pred_var_mod)]
    df_mod = as.data.frame(df_mod[, pred_var_mod])
    colnames(df_mod) = c(pred_var_mod)
    n = nrow(df_mod)
    fit = paste(syn_var, "~", 
                ifelse(length(pred_var_mod_ns) > 0, paste(pred_var_mod_ns, collapse = "+"), "1"))
    print(fit)
    fit_rq = rq(fit, data = df_mod, tau = tau)
    
    if (substr(syn_var, 5, 8) == "1976"){
      if(substr(syn_var, 1, 3) == "emp"){
        scale = scale_emp_fy
      } else {
        scale = scale_pay_fy
      }
      
      all_beta = list()
      
      
      # due to the characteristics of the KNG function, sometimes MCMC can move around too much
      # adding upper and lower bound helps with convergence.
      # upper bound is [arbitrarily] selected to be a high enough bound
      # lower bound is set as 0, because the number of employees cannot be negative
      ub = 10000000
      lb = 0
      
      all_x = lapply(syn_mod,  `[`, i =, j = pred_var_mod_ns)
      all_x = lapply(all_x, as.matrix)
      
      # median_eps, total_eps, and main_tau_eps indicate the percentage of privacy-loss
      # budget allocated to the median, the main quantiles, and the total procedure.
      
      temp = originalKNG(data = df_mod, total_eps = ep*0.75, tau = tau, nbatch = runs,
                         scale = scale[[1]])
      all_beta[[1]] = temp[[1]]
      
      temp = stepwiseKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.6, 
                         tau = tau, scale = scale[[2]], nbatch = runs, method = "fixed", 
                         lb = lb, ub = ub, formula = fit)
      all_beta[[2]] = temp[[1]]
      
      temp = stepwiseKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.8, 
                         tau = tau, scale = scale[[3]], nbatch = runs, method = "varying", 
                         lb = lb, ub = ub, formula = fit)
      
      all_beta[[3]] = temp[[1]]
      
      temp = sandwichKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.6,
                         main_tau_eps = 0.7, tau = tau,
                         main_tau = main_tau, scale = scale[[4]][1], sw_scale = scale[[4]][2], nbatch = runs, 
                         method = "fixed", lb = lb, ub = ub, formula = fit)
      
      all_beta[[4]] = temp[[1]]
      
      temp = sandwichKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.8,
                         main_tau_eps = 0.8, tau = tau,
                         main_tau = main_tau, scale = scale[[5]][1], sw_scale = scale[[5]][2], nbatch = runs, 
                         method = "varying", lb = lb, ub = ub, formula = fit)
      all_beta[[5]] = temp[[1]]
      
      all_beta[[6]] = coef(fit_rq)
      
      X = mapply(cbind, rep(list(matrix(1, nrow = n)), 6), all_x, SIMPLIFY = FALSE)
      synall = mapply(syndata, beta_result = all_beta, x = X, 
                      MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
      syn_mod = mapply(cbind, syn_mod, synall, SIMPLIFY = FALSE)
      
    } else {
      if(substr(syn_var, 1, 3) == "emp"){
        scale = scale_emp_oy
      } else {
        scale = scale_pay_oy
      }
      
      all_beta = list()
      
      ub = 10*max(df_mod[, ncol(df_mod)])
      lb = 0
      all_x = lapply(syn_mod,  `[`, i =, j = pred_var_mod_ns)
      all_x = lapply(all_x, as.matrix)
      
      temp = originalKNG(data = df_mod, total_eps = ep*0.75, tau = tau, nbatch = runs,
                         scale = scale[[1]])
      all_beta[[1]] = temp[[1]]
      
      temp = stepwiseKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.6, 
                         tau = tau, scale = scale[[2]], nbatch = runs, method = "fixed", 
                         lb = lb, ub = ub, formula = fit, check_data = all_x[[2]])
      all_beta[[2]] = temp[[1]]
      
      temp = stepwiseKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.8, 
                         tau = tau, scale = scale[[3]], nbatch = runs, method = "varying", 
                         lb = lb, ub = ub, formula = fit, check_data = all_x[[3]])
      
      all_beta[[3]] = temp[[1]]
      
      
      temp = sandwichKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.6,
                         main_tau_eps = 0.7, tau = tau, check_data = all_x[[4]],
                         main_tau = main_tau, scale = scale[[4]][1], sw_scale = scale[[4]][2], 
                         nbatch = runs, 
                         method = "fixed", lb = lb, ub = ub, formula = fit)
      
      all_beta[[4]] = temp[[1]]
      
      temp = sandwichKNG(data = df_mod, total_eps = ep*0.75, median_eps = 0.8,
                         main_tau_eps = 0.8, tau = tau, check_data = all_x[[5]],
                         main_tau = main_tau, scale = scale[[5]][1], sw_scale = scale[[5]][2], 
                         nbatch = runs, 
                         method = "varying", lb = lb, ub = ub, formula = fit)
      all_beta[[5]] = temp[[1]]
      
      all_beta[[6]] = coef(fit_rq)
      X = mapply(cbind, rep(list(matrix(1, nrow = n)), 6), all_x, SIMPLIFY = FALSE)
      synall = mapply(syndata, beta_result = all_beta, x = X, 
                      MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
      syn_mod = mapply(cbind, syn_mod, synall, SIMPLIFY = FALSE)
    }
    
    name_var = c(names(syn_mod[[1]])[-j], syn_var)
    syn_mod = lapply(syn_mod, `colnames<-`, name_var)
    
    df_birth = df_birth[, c(1:j)]
    pred_var_birth = syn_var
    pred_var_birth_ns = pred_var_birth[-length(pred_var_birth)]
    df_birth = as.data.frame(df_birth[, pred_var_birth])
    colnames(df_birth) = c(pred_var_birth)
    n = nrow(df_birth)
    fit = paste(syn_var, "~", 
                ifelse(length(pred_var_birth_ns) > 0, paste(pred_var_birth_ns, collapse = "+"), "1"))
    print(fit)
    fit_rq = rq(fit, data = df_birth, tau = tau)
    all_beta_birth = list()
    
    if(substr(syn_var, 1, 3) == "emp"){
      scale = scale_emp_fy
    } else {
      scale = scale_pay_fy
    }
    
    # due to the characteristics of the KNG function, sometimes MCMC can move around too much
    # adding upper and lower bound helps with convergence.
    # upper bound is [arbitrarily] selected to be a high enough bound
    # lower bound is set as 0, because the number of employees cannot be negative
    
    ub = 1000000
    lb = 0
    
    temp = originalKNG(data = df_birth, total_eps = ep*0.25, tau = tau, nbatch = runs,
                       scale = scale[[1]])
    all_beta_birth[[1]] = temp[[1]]
    
    temp = stepwiseKNG(data = df_birth, total_eps = ep*0.25, median_eps = 0.6, 
                       tau = tau, scale = scale[[2]], nbatch = runs, method = "fixed", 
                       lb = lb, ub = ub, formula = fit)
    all_beta_birth[[2]] = temp[[1]]
    
    temp = stepwiseKNG(data = df_birth, total_eps = ep*0.25, median_eps = 0.8, 
                       tau = tau, scale = scale[[3]], nbatch = runs, method = "varying", 
                       lb = lb, ub = ub, formula = fit)
    
    all_beta_birth[[3]] = temp[[1]]
    
    temp = sandwichKNG(data = df_birth, total_eps = ep*0.25, median_eps = 0.6,
                       main_tau_eps = 0.7, tau = tau, nbatch = runs, 
                       main_tau = main_tau, scale = scale[[4]][1], sw_scale = scale[[4]][2],  
                       method = "fixed", lb = lb, ub = ub, formula = fit)
    
    all_beta_birth[[4]] = temp[[1]]
    
    temp = sandwichKNG(data = df_birth, total_eps = ep*0.25, median_eps = 0.8,
                       main_tau_eps = 0.8, tau = tau, nbatch = runs,
                       main_tau = main_tau, scale = scale[[5]][1], sw_scale = scale[[5]][2],  
                       method = "varying", lb = lb, ub = ub, formula = fit)
    all_beta_birth[[5]] = temp[[1]]
    
    all_beta_birth[[6]] = coef(fit_rq)
    
    X = rep(list(matrix(1, nrow = n)), 6)
    synallbirth = mapply(syndata, beta_result = all_beta_birth, x = X, 
                         MoreArgs = list(sample(1:length(tau), n, replace = TRUE)), SIMPLIFY = FALSE)
    syn_birth = mapply(cbind, syn_birth, synallbirth, SIMPLIFY = FALSE)
    syn_birth = lapply(syn_birth, `colnames<-`, name_var)
    
    n = nrow(syn_na[[1]])
    na_data = rep(list(rep(NA, n)), 6)
    syn_na = mapply(cbind, syn_na, na_data, SIMPLIFY = FALSE)
    syn_na = lapply(syn_na, `colnames<-`, name_var)
    
    syn_final = mapply(rbind, syn_mod, syn_birth, syn_na, SIMPLIFY = FALSE) 
  }
  
  file_name=paste0(wd, "/amh1/", sic, "/SIC", sic, "_seed", seed, "_KNG.csv")
  write.csv(as.data.frame(syn_final[[1]]),file=file_name,append=FALSE,row.names = FALSE)
  file_name=paste0(wd, "/amh1/", sic, "/SIC", sic, "_seed", seed, "_StepFixed.csv")
  write.csv(as.data.frame(syn_final[[2]]),file=file_name,append=FALSE,row.names = FALSE)
  file_name=paste0(wd, "/amh1/", sic, "/SIC", sic, "_seed", seed, "_StepVarying.csv")
  write.csv(as.data.frame(syn_final[[3]]),file=file_name,append=FALSE,row.names = FALSE)
  file_name=paste0(wd, "/amh1/", sic, "/SIC", sic, "_seed", seed, "_SWFixed.csv")
  write.csv(as.data.frame(syn_final[[4]]),file=file_name,append=FALSE,row.names = FALSE)
  file_name=paste0(wd, "/amh1/", sic, "/SIC", sic, "_seed", seed, "_SWVarying.csv")
  write.csv(as.data.frame(syn_final[[5]]),file=file_name,append=FALSE,row.names = FALSE)
  file_name=paste0(wd, "/amh1/", sic, "/SIC", sic, "_seed", seed, "_NonPriv.csv")
  write.csv(as.data.frame(syn_final[[6]]),file=file_name,append=FALSE,row.names = FALSE)
  file_name=paste0(wd, "/amh1/", sic, "/SIC", sic, "_seed", seed, "_Truth.csv")
  write.csv(as.data.frame(data),file=file_name,append=FALSE,row.names = FALSE)
}

# run it in parallel to speed up runtime
source("../../config.R", echo = TRUE)
library(doParallel, lib = dir.Rpackages)
num_cores= detectCores()-1 #use all available core but 1
workers=makeCluster(num_cores,type="SOCK")
registerDoParallel(workers)

Out = foreach(seed=c(1:5), .errorhandling='stop') %dopar% {
  source(paste0(wd, "/config.R"))
  source(paste0(wd, "/amh1/functions_final.R"))
  
	library(MASS, lib = dir.Rpackages)
	library(tictoc, lib = dir.Rpackages)
	library(FNN, lib = dir.Rpackages)
	library(SparseM, lib = dir.Rpackages)
	library(ggplot2, lib = dir.Rpackages)
	library(quantreg, lib = dir.Rpackages)
	library(rpart,lib=dir.Rpackages)
	library(mvtnorm, lib = dir.Rpackages)
  library(dplyr, lib = dir.Rpackages)

	set.seed(seed)
	syn_all_methods(sic = 178, runs = 1000, endcol = 29)
}
stopCluster(workers)

