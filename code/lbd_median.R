# This code is used to plot the trend of gross employment, job creation rate,
# and net job creation rate using the generated synthetic data.sets

rm(list = ls())
wd = "<project directory>/amh1/"  # This was the path on the Cornell Synthetic Data Server
											# It is not a valid path on a Census Bureau server.
setwd(wd)
source("../config.R", echo = TRUE)
library(ggplot2, lib = dir.Rpackages)
library(digest, lib = dir.Rpackages)
library(data.table, lib = dir.Rpackages)
set.seed(1)
median_ci = function(x) {
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median)
  return(quantile(bootmed, c(.025, 0.975)))
}

#sic = c(178, 239, 542, 829)  # These are SIC codes.
sic = 542  					  # This is an SIC code.
methods = c("StepVarying", "SWVarying")

setwd(paste0(wd,'/amh1/', sic, "/"))
df = fread(paste0("SIC", sic, "_seed1_Truth.csv"))
drops = c("firstyear", "lastyear", "mu")
df_yr = df[, !drops, with = F]
df_yr = apply(df_yr, 2, sum, na.rm = T)
df_yr
df_p = NULL
df_ub = NULL
df_lb = NULL
for (i in 1:2){
  tab = NULL
  for (seed in 1:5){
    syn = fread(paste0("SIC", sic, "_seed", seed, "_", methods[i], ".csv"))
    syn_yr = syn[, !drops, with = F]
    syn_yr = apply(syn_yr, 2, sum, na.rm = T)
    tab = rbind(tab, syn_yr)
  }
  tab_p = apply(tab, 2, median)
  tmp = apply(tab, 2, median_ci)
  df_lb = rbind(df_lb, tmp[1,])
  df_ub = rbind(df_ub, tmp[2,])
  df_p = rbind(df_p, tab_p) 
}
rownames(df_p) = methods
rownames(df_ub) = methods
rownames(df_lb) = methods

tmp = NULL
for (j in 1:5){
  p = fread(paste0("<project directory>/spec989/Desktop/final/Pistner/QR", j, "/QR_SynLBD_", sic, ".csv")) # This was the path on the Cornell Synthetic Data Server
                                                                                            # It is not a valid path on a Census Bureau server.
  p = p[, !c(drops, "sic3"), with = F]
  tmp = rbind(tmp, apply(p, 2, sum, na.rm = T)[1:25])
  
}
df_p = as.data.table(rbind(df_p, apply(tmp, 2, median), df_yr))
tmp1 = apply(tmp, 2, median_ci)
df_lb = as.data.table(rbind(df_lb, tmp1[1,], df_yr))
df_ub = as.data.table(rbind(df_ub, tmp1[2,], df_yr))

df_p$methods = c("StepVarying", "SWVarying", "Pistner et al (2018)", "Raw Data")
df_lb$methods = c("StepVarying", "SWVarying", "Pistner et al (2018)", "Raw Data")
df_ub$methods = c("StepVarying", "SWVarying", "Pistner et al (2018)", "Raw Data")
measure_variables = paste0('emp_', c(1976:2000))
df_p = melt(df_p, id.vars = c("methods"), measure.vars = measure_variables)
df_lb = melt(df_lb, id.vars = c("methods"), measure.vars = measure_variables)
df_ub = melt(df_ub, id.vars = c("methods"), measure.vars = measure_variables)

df_final = merge(df_p, df_lb, by = c('methods', 'variable'), all_x = T)
df_final = merge(df_final, df_ub, by = c('methods', 'variable'), all_x = T)
colnames(df_final) = c('methods', 'variable', 'median', 'min_val', 'max_val')
df_final$year = as.numeric(substr(df_final$variable, 5,8))


gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
# Dynamically generate default color values, but have Raw Data = "black"
adj_names = sort(setdiff(unique(df_final$methods), "Raw Data"))
values = gg_color_hue(length(adj_names))
names(values) = adj_names
values = c(values, c('Raw Data'="black"))
thickness = c(1, 1, 1, 1.5)
names(thickness) = c(adj_names, 'Raw Data')

palette = c("#000000", "#CC79A7", "#E69F00", "#0072B2")
names(palette) = c('Raw Data', adj_names)


ggplot(data = df_final, aes(x = year, group = methods)) + 
  geom_line(aes(y = median, color = methods, size = methods)) + 
  geom_ribbon(aes(y = median, ymin = min_val, ymax = max_val, fill = methods), alpha = .3) +
  labs(title=paste("Gross Employment by Year for Industry Code", sic), 
       x ="Year", y = "Gross Employment (in thousands)", 
       color = 'Method', fill = 'Method', size = 'Method') +
  scale_x_continuous(breaks=seq(1980, 2000, by = 5)) + 
  scale_colour_manual(values=palette) +
  scale_fill_manual(values=palette) +
  scale_size_manual(values = thickness) +
  theme(axis.text = element_text(size=13),
        legend.position = 'bottom',
        plot.title = element_text(size = 14, face = "bold"),
        legend.title=element_text(size=13), 
        legend.text=element_text(size=12),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'),
        axis.title=element_text(size=13))



###############################################################################
#job creation rate

# total emp across 4 industries plot
rm(list = ls())

setwd("<project directory>/spec989/Desktop/project/") # This was the path on the Cornell Synthetic Data Server
                                        # It is not a valid path on a Census Bureau server.
source("config.R", echo = TRUE)
library(ggplot2, lib = dir.Rpackages)
library(data.table, lib = dir.Rpackages)
drops = c("firstyear", "lastyear", "mu")
sic = c(178, 239, 542, 829)  # These are SIC codes.


jc = list()
njc = list()

methods = c("StepVarying", "SWVarying")
for (seed in 1:5){
  tab_jc = matrix(NA, nrow = 4, ncol = 24)
  tab_jd = matrix(NA, nrow = 4, ncol = 24)
  tab_n = matrix(NA, nrow = 4, ncol = 24)
  df = NULL
  for (k in 1:length(sic)){
    # The line below was for a path on the Cornell Synthetic Data Server.
    # It is not a valid path on a Census Bureau server.
    df = rbind(df, fread(paste0("<project directory>/spec989/Desktop/project/amh1/", sic[k], "/SIC", sic[k], "_seed", seed,"_Truth.csv")))
  }
  df[df == ""] = NA
  df = df[, !drops, with = FALSE]
  
  for (t in 2:25){
    keep = c((t-1):t)
    df_t = na.omit(df[, ..keep])
    z_et = rowSums(df_t)/2
    z = sum(z_et)
    diff = df_t[,2] - df_t[,1]
    diff = cbind(diff, rep(0, length(diff)))
    tab_jc[1, (t-1)] = sum(abs(apply(diff, 1, max))/z)
    tab_jd[1, (t-1)] = sum(abs(apply(diff, 1, min))/z)
    tab_n[1, (t-1)] = tab_jc[1, (t-1)] - tab_jd[1, (t-1)]
    
  }
  for (i in 1:length(methods)){
    df = NULL
    for (k in 1:length(sic)){
      # The line below was for a path on the Cornell Synthetic Data Server.
      # It is not a valid path on a Census Bureau server.
      tab = fread(paste0("<project directory>/spec989/Desktop/project/amh1/", sic[k], "/SIC", sic[k], "_seed", seed, "_", methods[i], ".csv"))
      tab = tab[, !drops, with = F]
      df = rbind(df, tab)
    }
    
    for (t in 2:25){
      keep = c((t-1):t)
      df_t = na.omit(df[, ..keep])
      z_et = rowSums(df_t)/2
      z = sum(z_et)
      diff = df_t[,2] - df_t[,1]
      diff = cbind(diff, rep(0, length(diff)))
      tab_jc[(i+1), (t-1)] = sum(abs(apply(diff, 1, max))/z)
      tab_jd[(i+1), (t-1)] = sum(abs(apply(diff, 1, min))/z)
      tab_n[(i+1), (t-1)] = tab_jc[(i+1), (t-1)] - tab_jd[(i+1), (t-1)]
      
    }
    
  }
  
  
  tab = NULL
  for (k in 1:length(sic)){
    # The line below was for a path on the Cornell Synthetic Data Server.
    # It is not a valid path on a Census Bureau server.
    p = fread(paste0("<project directory>/spec989/Desktop/project/amh1/Pistner/QR", seed,"/QR_SynLBD_", sic[k], ".csv"))
    p = p[, !c(drops, "sic3"), with = F]
    tab = rbind(tab, p[, 1:25])
  }
  df = tab
  for (t in 2:25){
    keep = c((t-1):t)
    df_t = na.omit(df[, ..keep])
    z_et = rowSums(df_t)/2
    z = sum(z_et)
    diff = df_t[,2] - df_t[,1]
    diff = cbind(diff, rep(0, length(diff)))
    tab_jc[4, (t-1)] = sum(abs(apply(diff, 1, max))/z)
    tab_jd[4, (t-1)] = sum(abs(apply(diff, 1, min))/z)
    tab_n[4, (t-1)] = tab_jc[4, (t-1)] - tab_jd[4, (t-1)]
    
  }
  
  jc[[seed]] = tab_jc*100
  njc[[seed]] = tab_n*100
}

######################## 
# plot job creation rate by year

plot_jc = matrix(NA, ncol = 24, nrow = 4)
plot_jc_sd = matrix(NA, ncol = 24, nrow = 4)
for (i in 1:4){
  tmp = rbind(jc[[1]][i,], jc[[2]][i,], jc[[3]][i,], jc[[4]][i,], jc[[5]][i,])
  plot_jc[i,] = apply(tmp, 2, mean)
  plot_jc_sd[i,] = apply(tmp, 2, sd)
}

vars = paste0("emp_", c(1977:2000))
plot_jc = as.data.table(plot_jc)
plot_jc_sd = as.data.table(plot_jc_sd)
colnames(plot_jc) = vars
colnames(plot_jc_sd) = vars
plot_jc$methods= c("Raw Data", "Stepwise Varying", 
                   "SW Varying", "Pistner et al (2018)") 
plot_jc_sd$methods = c("Raw Data", "Stepwise Varying", 
                       "SW Varying", "Pistner et al (2018)") 

plot_jc= melt(plot_jc, id.vars = c("methods"), measure.vars = vars)
plot_jc_sd = melt(plot_jc_sd, id.vars = c("methods"), measure.vars = vars)
plot_final = merge(plot_jc, plot_jc_sd, by = c('methods', 'variable'), all_x = T)
colnames(plot_final) = c('methods', 'variable', 'mean', 'sd')
plot_final$min_val = plot_final$mean - plot_final$sd
plot_final$max_val = plot_final$mean + plot_final$sd
plot_final$year = as.numeric(substr(plot_final$variable, 5,8))


# Dynamically generate default color values, but have Raw Data = "black"
adj_names = sort(setdiff(unique(plot_final$methods), "Raw Data"))
thickness = c(1, 1, 1, 1.5)
names(thickness) = c(adj_names, 'Raw Data')

palette = c("#000000", "#CC79A7", "#E69F00", "#0072B2")
names(palette) = c('Raw Data', adj_names)

ggplot(data = plot_final, aes(x = year, group = methods)) + 
  geom_line(aes(y = mean, color = methods, size = methods)) + 
  geom_ribbon(aes(y = mean, ymin = min_val, ymax = max_val, fill = methods), alpha = .3) +
  labs(title="Job Creation Rate by Year", 
       x ="Year", y = "Job Creation Rate (Percent)", 
       color = 'Method', fill = 'Method', size = 'Method') +
  scale_x_continuous(breaks=seq(1980, 2000, by = 5)) + 
  scale_colour_manual(values=palette) +
  scale_fill_manual(values=palette) +
  scale_size_manual(values = thickness) +
  theme(axis.text = element_text(size=13),
        legend.position = 'bottom',
        plot.title = element_text(size = 14, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'),
        axis.title=element_text(size=13))

######################## 
# plot net job creation rate by year

plot_njc = matrix(NA, ncol = 24, nrow = 4)
plot_njc_sd = matrix(NA, ncol = 24, nrow = 4)
for (i in 1:4){
  tmp = rbind(njc[[1]][i,], njc[[2]][i,], njc[[3]][i,], njc[[4]][i,], njc[[5]][i,])
  plot_njc[i,] = apply(tmp, 2, mean)
  plot_njc_sd[i,] = apply(tmp, 2, sd)
}

vars = paste0("emp_", c(1977:2000))
plot_njc = as.data.table(plot_njc)
plot_njc_sd = as.data.table(plot_njc_sd)
colnames(plot_njc) = vars
colnames(plot_njc_sd) = vars
plot_njc$methods= c("Raw Data", "Stepwise Varying", 
                    "SW Varying", "Pistner et al (2018)") 
plot_njc_sd$methods = c("Raw Data", "Stepwise Varying", 
                        "SW Varying", "Pistner et al (2018)") 

plot_njc= melt(plot_njc, id.vars = c("methods"), measure.vars = vars)
plot_njc_sd = melt(plot_njc_sd, id.vars = c("methods"), measure.vars = vars)
plot_final = merge(plot_njc, plot_njc_sd, by = c('methods', 'variable'), all_x = T)
colnames(plot_final) = c('methods', 'variable', 'mean', 'sd')
plot_final$min_val = plot_final$mean - plot_final$sd
plot_final$max_val = plot_final$mean + plot_final$sd
plot_final$year = as.numeric(substr(plot_final$variable, 5,8))


palette = c("#000000", "#CC79A7", "#E69F00", "#0072B2")
names(palette) = c('Raw Data', adj_names)

ggplot(data = plot_final, aes(x = year, group = methods)) + 
  geom_line(aes(y = mean, color = methods, size = methods)) + 
  geom_ribbon(aes(y = mean, ymin = min_val, ymax = max_val, fill = methods), alpha = .3) +
  labs(title="Net Job Creation Rate by Year", 
       x ="Year", y = "Net Job Creation Rate (Percent)", 
       color = 'Method', fill = 'Method', size = 'Method') +
  scale_x_continuous(breaks=seq(1980, 2000, by = 5)) + 
  scale_colour_manual(values=palette) +
  scale_fill_manual(values=palette) +
  scale_size_manual(values = thickness) +
  theme(axis.text = element_text(size=13),
        legend.position = 'bottom',
        plot.title = element_text(size = 14, face = "bold"),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=13),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'),
        axis.title=element_text(size=13))
