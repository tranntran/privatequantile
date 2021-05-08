# This file can be used to plot any results from icds_distance_sandwich_data.R.

library(ggplot2)
library(data.table)
rm(list = ls())
filename = "sandwich_data_sd_1_e0.5"
load(paste("output/", filename, ".Rdata", sep = ""))

tau = c(seq(0.05, 0.95, 0.05), 0.99)
allocate_eps = seq(0.1, 0.9, 0.1)
reps = 100


png(paste("plot/distance/", filename, "_distance_intercept.png", sep = ""), 
    width = 600, height = 400)
colnames(sandwich_data_int) = tau
sandwich_data_int = as.data.table(sandwich_data_int)
sandwich_data_int$Eps = paste0(allocate_eps*100, '%')
sandwich_data_int = melt(id.vars = 21, sandwich_data_int)
colnames(sandwich_data_int) = c("Eps", "Quantile", "Value")
sandwich_data_int$Eps = as.factor(sandwich_data_int$Eps)
sandwich_data_int$Quantile = as.numeric(as.character(sandwich_data_int$Quantile))
ggplot(sandwich_data_int, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau)  +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation",
       title = "L1 Distance to the True Intercept by Epsilon Allocation to Main Quantiles - Total Epsilon = 0.5",
       subtitle = "Median Epsilon = 80% of Main Quantiles Epsilon")
dev.off()

png(paste("plot/distance/", filename, "_distance_slope.png", sep = ""), 
    width = 600, height = 400)
colnames(sandwich_data_slope) = tau
sandwich_data_slope = as.data.table(sandwich_data_slope)
sandwich_data_slope$Eps = paste0(allocate_eps*100, '%')
sandwich_data_slope = melt(id.vars = 21, sandwich_data_slope)
colnames(sandwich_data_slope) = c("Eps", "Quantile", "Value")
sandwich_data_slope$Eps = as.factor(sandwich_data_slope$Eps)
sandwich_data_slope$Quantile = as.numeric(as.character(sandwich_data_slope$Quantile))
ggplot(sandwich_data_slope, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation",
       title = "L1 Distance to the True Slope by Epsilon Allocation to the Main Quantiles - Total Epsilon = 0.5",
       subtitle = "Median Epsilon = 80% of Main Quantiles Epsilon")
dev.off()



colnames(sandwich_data_l2) = tau
sandwich_data_l2 = as.data.table(sandwich_data_l2)
sandwich_data_l2$Eps = paste0(allocate_eps*100, '%')
sandwich_data_l2 = melt(id.vars = c(21), sandwich_data_l2)
sandwich_data_l2 = sandwich_data_l2[order(sandwich_data_l2$Eps)]
colnames(sandwich_data_l2) = c("Eps", "Quantile", "Value")
sandwich_data_l2$Eps = as.factor(sandwich_data_l2$Eps)
sandwich_data_l2$Quantile = as.numeric(as.character(sandwich_data_l2$Quantile))

png(paste("plot/distance/", filename, "_distance_meanl2.png", sep = ""), 
    width = 600, height = 400)
ggplot(sandwich_data_l2, aes(x = Quantile, log10(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ylab("log10(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation",
       title = "L2 Distance to the Truth by Epsilon Allocation to the Main Quantiles - Total Epsilon - 0.5",
       subtitle = "Median Epsilon = 80% of Main Quantiles Epsilon")
dev.off()


colnames(sd_l2) = tau
sd_l2 = as.data.table(sd_l2)
sd_l2$Eps = paste0(allocate_eps*100, '%')
sd_l2_m = melt(sd_l2, id.vars = "Eps")
colnames(sd_l2_m) = c("Eps", "Quantile", "Sd")
sd_l2_m$Eps = as.factor(sd_l2_m$Eps)
sd_l2_m$Quantile = as.numeric(as.character(sd_l2_m$Quantile))
setkey(sandwich_data_l2, Eps, Quantile)
setkey(sd_l2_m, Eps, Quantile)
sandwich_data_l2 = merge(sandwich_data_l2, sd_l2_m, all.x = T)
sandwich_data_l2$Sd_n = sandwich_data_l2$Sd/sqrt(reps)

png(paste("plot/", filename, "_distance_l2se.png", sep = ""), 
    width = 600, height = 400)
ggplot(sandwich_data_l2, aes(x = Quantile, Value, ymin=Value-2*Sd_n, 
                                ymax=Value+2*Sd_n, fill = Eps)) + 
  facet_wrap(~Eps) +
  geom_line(aes(color = Eps), size = 1) +
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(breaks= tau[seq(1, 20, 3)]) +
  ggtitle("L2 Distance to the Truth by Epsilon Allocation to the Main Quantiles") + 
  ylab("Mean Distance") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation") + ylim(c(0, 50))
dev.off()



