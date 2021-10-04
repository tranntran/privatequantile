# This file can be used to plot any results from icds_distance_sandwich_fixed.R.
# Change total eps and median eps in title and variable reps
library(ggplot2)
library(data.table)
rm(list = ls())
filename = "sandwich_fixed_sd_123_e0.5"
load(paste("output/", filename, ".Rdata", sep = ""))

tau = c(seq(0.05, 0.95, 0.05), 0.99)
allocate_eps = c(1/length(tau), seq(0.1, 0.9, 0.1))
reps = 100

# this section plots the distance to the truth for intercept
png(paste("plot/distance/", filename, "_distance_intercept.png", sep = ""), 
    width = 600, height = 400)
colnames(sandwich_fixed_int) = tau
sandwich_fixed_int = as.data.table(sandwich_fixed_int)
sandwich_fixed_int$Eps = paste0(allocate_eps*100, '%')
sandwich_fixed_int = melt(id.vars = 21, sandwich_fixed_int)
colnames(sandwich_fixed_int) = c("Eps", "Quantile", "Value")
sandwich_fixed_int$Eps = factor(sandwich_fixed_int$Eps, paste0(allocate_eps*100, '%'))
sandwich_fixed_int$Quantile = as.numeric(as.character(sandwich_fixed_int$Quantile))
ggplot(sandwich_fixed_int, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau)  +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation", 
       title = "L1 Distance to the True Intercept by Epsilon Allocation to Main Quantiles - Total Epsilon = 0.5",
       subtitle = "Median Epsilon = 60% of Main Quantiles Epsilon")
dev.off()

# this section plots the distance to the truth for slope
png(paste("plot/distance/", filename, "_distance_slope.png", sep = ""), 
    width = 600, height = 400)
colnames(sandwich_fixed_slope) = tau
sandwich_fixed_slope = as.data.table(sandwich_fixed_slope)
sandwich_fixed_slope$Eps = paste0(allocate_eps*100, '%')
sandwich_fixed_slope = melt(id.vars = 21, sandwich_fixed_slope)
colnames(sandwich_fixed_slope) = c("Eps", "Quantile", "Value")
sandwich_fixed_slope$Eps = factor(sandwich_fixed_slope$Eps, paste0(allocate_eps*100, '%'))
sandwich_fixed_slope$Quantile = as.numeric(as.character(sandwich_fixed_slope$Quantile))
ggplot(sandwich_fixed_slope, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation", 
       title = "L1 Distance to the True Slope by Epsilon Allocation to the Main Quantiles - Total Epsilon = 0.5",
       subtitle = "Median Epsilon = 60% of Main Quantiles Epsilon")
dev.off()


# this section plots the l2 distance to the truth
colnames(sandwich_fixed_l2) = tau
sandwich_fixed_l2 = as.data.table(sandwich_fixed_l2)
sandwich_fixed_l2$Eps = paste0(allocate_eps*100, '%')
sandwich_fixed_l2 = melt(id.vars = c(21), sandwich_fixed_l2)
sandwich_fixed_l2 = sandwich_fixed_l2[order(sandwich_fixed_l2$Eps)]
colnames(sandwich_fixed_l2) = c("Eps", "Quantile", "Value")
sandwich_fixed_l2$Eps = factor(sandwich_fixed_l2$Eps, paste0(allocate_eps*100, '%'))
sandwich_fixed_l2$Quantile = as.numeric(as.character(sandwich_fixed_l2$Quantile))

png(paste("plot/distance/", filename, "_distance_meanl2.png", sep = ""), 
    width = 600, height = 400)
ggplot(sandwich_fixed_l2, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation",
       title = "L2 Distance to the Truth by Epsilon Allocation to the Main Quantiles - Total Epsilon = 0.5",
       subtitle = "Median Epsilon = 60% of Main Quantiles Epsilon")
dev.off()

# this section plots the l2 distance to the truth with se
# the first part of the code should only be run once
# i.e. if l2 distance to the truth has been plot, skip the first block of the code
colnames(sandwich_fixed_l2) = tau
sandwich_fixed_l2 = as.data.table(sandwich_fixed_l2)
sandwich_fixed_l2$Eps = paste0(allocate_eps*100, '%')
sandwich_fixed_l2 = melt(id.vars = c(21), sandwich_fixed_l2)
sandwich_fixed_l2 = sandwich_fixed_l2[order(sandwich_fixed_l2$Eps)]
colnames(sandwich_fixed_l2) = c("Eps", "Quantile", "Value")
sandwich_fixed_l2$Eps = factor(sandwich_fixed_l2$Eps, paste0(allocate_eps*100, '%'))
sandwich_fixed_l2$Quantile = as.numeric(as.character(sandwich_fixed_l2$Quantile))

colnames(sd_l2) = tau
sd_l2 = as.data.table(sd_l2)
sd_l2$Eps = paste0(allocate_eps*100, '%')
sd_l2_m = melt(sd_l2, id.vars = "Eps")
colnames(sd_l2_m) = c("Eps", "Quantile", "Sd")
sd_l2_m$Eps = factor(sd_l2_m$Eps, paste0(allocate_eps*100, '%'))
sd_l2_m$Quantile = as.numeric(as.character(sd_l2_m$Quantile))
setkey(sandwich_fixed_l2, Eps, Quantile)
setkey(sd_l2_m, Eps, Quantile)
sandwich_fixed_l2 = merge(sandwich_fixed_l2, sd_l2_m, all.x = T)
sandwich_fixed_l2$Sd_n = sandwich_fixed_l2$Sd/sqrt(reps)

png(paste("plot/", filename, "_distance_l2se.png", sep = ""), 
    width = 600, height = 400)
ggplot(sandwich_fixed_l2, aes(x = Quantile, log(Value), ymin=log(Value-2*Sd_n), 
                             ymax=log(Value+2*Sd_n), fill = Eps)) + 
  facet_wrap(~Eps) +
  geom_line(aes(color = Eps), size = 1) +
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(breaks= tau[seq(1, 20, 3)]) +
  ggtitle("L2 Distance to the Truth by Epsilon Allocation to the Main Quantiles") + 
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation")
dev.off()



