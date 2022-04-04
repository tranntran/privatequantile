# This file can be used to plot any results from icds_distance_stepwise_varying.R.

library(ggplot2)
library(data.table)
rm(list = ls())
filename = "stepwise_varying_sd_123_e0.5"
load(paste("output/", filename, ".Rdata", sep = ""))

tau = c(seq(0.05, 0.95, 0.05), 0.99)
allocate_eps = c(1/length(tau), seq(0.1, 0.9, 0.1))
reps = 100


png(paste("plot/distance/", filename, "_distance_intercept.png", sep = ""), 
    width = 600, height = 400)
colnames(stepwise_varying_int) = tau
stepwise_varying_int = as.data.table(stepwise_varying_int)
stepwise_varying_int$Eps = paste0(allocate_eps*100, '%')
stepwise_varying_int = melt(id.vars = 21, stepwise_varying_int)
colnames(stepwise_varying_int) = c("Eps", "Quantile", "Value")
stepwise_varying_int$Eps = factor(stepwise_varying_int$Eps, paste0(allocate_eps*100, '%'))
stepwise_varying_int$Quantile = as.numeric(as.character(stepwise_varying_int$Quantile))
ggplot(stepwise_varying_int, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau)  +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation",
       title = "L1 Distance to the True Intercept by Epsilon Allocation to the Median - Total Eps = 0.5")
dev.off()

png(paste("plot/distance/", filename, "_distance_slope.png", sep = ""), 
    width = 600, height = 400)
colnames(stepwise_varying_slope) = tau
stepwise_varying_slope = as.data.table(stepwise_varying_slope)
stepwise_varying_slope$Eps = paste0(allocate_eps*100, '%')
stepwise_varying_slope = melt(id.vars = 21, stepwise_varying_slope)
colnames(stepwise_varying_slope) = c("Eps", "Quantile", "Value")
stepwise_varying_slope$Eps = factor(stepwise_varying_slope$Eps, paste0(allocate_eps*100, '%'))
stepwise_varying_slope$Quantile = as.numeric(as.character(stepwise_varying_slope$Quantile))
ggplot(stepwise_varying_slope, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation", 
       title = "L1 Distance to the True Slope by Epsilon Allocation to the Median - Total Epsilon = 0.5")
dev.off()



colnames(stepwise_varying_l2) = tau
stepwise_varying_l2 = as.data.table(stepwise_varying_l2)
stepwise_varying_l2$Eps = paste0(allocate_eps*100, '%')
stepwise_varying_l2 = melt(id.vars = c(21), stepwise_varying_l2)
stepwise_varying_l2 = stepwise_varying_l2[order(stepwise_varying_l2$Eps)]
colnames(stepwise_varying_l2) = c("Eps", "Quantile", "Value")
stepwise_varying_l2$Eps = factor(stepwise_varying_l2$Eps, paste0(allocate_eps*100, '%'))
stepwise_varying_l2$Quantile = as.numeric(as.character(stepwise_varying_l2$Quantile))

png(paste("plot/distance/", filename, "_distance_meanl2.png", sep = ""), 
    width = 600, height = 400)
ggplot(stepwise_varying_l2, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ggtitle("L2 Distance to the Truth by Epsilon Allocation to the Median") +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation", 
       title = "L2 Distance to the Truth by Epsilon Allocation to the Median - Total Epsilon = 0.5")
dev.off()


colnames(sd_l2) = tau
sd_l2 = as.data.table(sd_l2)
sd_l2$Eps = paste0(allocate_eps*100, '%')
sd_l2_m = melt(sd_l2, id.vars = "Eps")
colnames(sd_l2_m) = c("Eps", "Quantile", "Sd")
sd_l2_m$Eps = factor(sd_l2_m$Eps, paste0(allocate_eps*100, '%'))
sd_l2_m$Quantile = as.numeric(as.character(sd_l2_m$Quantile))
setkey(stepwise_varying_l2, Eps, Quantile)
setkey(sd_l2_m, Eps, Quantile)
stepwise_varying_l2 = merge(stepwise_varying_l2, sd_l2_m, all.x = T)
stepwise_varying_l2$Sd_n = stepwise_varying_l2$Sd/sqrt(reps)

png(paste("plot/", filename, "_distance_l2se.png", sep = ""), 
    width = 600, height = 400)
ggplot(stepwise_varying_l2, aes(x = Quantile, Value, ymin=Value-2*Sd_n, 
                                ymax=Value+2*Sd_n, fill = Eps)) + 
  facet_wrap(~Eps) +
  geom_line(aes(color = Eps), size = 1) +
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(breaks= tau[seq(1, 20, 3)]) +
  ggtitle("L2 Distance to the Truth by Epsilon Allocation to the Median") + 
  ylab("Mean Distance") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation") +
   ylim(c(0, 50))
dev.off()



