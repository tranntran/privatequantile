# This file can be used to plot any results from icds_distance_sandwich_koenker.R.

library(ggplot2)
library(data.table)
rm(list = ls())
filename = "sandwich_koenker_sd_2_e0.1"
load(paste("output/", filename, ".Rdata", sep = ""))

tau = c(seq(0.05, 0.95, 0.05), 0.99)
allocate_eps = seq(0.1, 0.9, 0.1)


png(paste("plot/", filename, "_distance_intercept.png", sep = ""), 
    width = 600, height = 400)
colnames(sandwich_koenker_int) = tau
sandwich_koenker_int = as.data.table(sandwich_koenker_int)
sandwich_koenker_int$Eps = allocate_eps
sandwich_koenker_int = melt(id.vars = 21, sandwich_koenker_int)
colnames(sandwich_koenker_int) = c("Eps", "Quantile", "Value")
sandwich_koenker_int$Eps = as.factor(sandwich_koenker_int$Eps)
sandwich_koenker_int$Quantile = as.numeric(as.character(sandwich_koenker_int$Quantile))
ggplot(sandwich_koenker_int, aes(x = Quantile, log(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau)  +
  ggtitle("Distance to the True Intercept by Epsilon Allocation to the Median") +
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation") + 
  ylim(c(1.25, 3.25))
dev.off()

png(paste("plot/", filename, "_distance_slope.png", sep = ""), 
    width = 600, height = 400)
colnames(sandwich_koenker_slope) = tau
sandwich_koenker_slope = as.data.table(sandwich_koenker_slope)
sandwich_koenker_slope$Eps = allocate_eps
sandwich_koenker_slope = melt(id.vars = 21, sandwich_koenker_slope)
colnames(sandwich_koenker_slope) = c("Eps", "Quantile", "Value")
sandwich_koenker_slope$Eps = as.factor(sandwich_koenker_slope$Eps)
sandwich_koenker_slope$Quantile = as.numeric(as.character(sandwich_koenker_slope$Quantile))
ggplot(sandwich_koenker_slope, aes(x = Quantile, log(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ggtitle("Distance to the True Slope by Epsilon Allocation to the Median") +
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation") +
  ylim(c(-2, 2))
dev.off()



colnames(sandwich_koenker_l2) = tau
sandwich_koenker_l2 = as.data.table(sandwich_koenker_l2)
sandwich_koenker_l2$Eps = allocate_eps
sandwich_koenker_l2 = melt(id.vars = c(21), sandwich_koenker_l2)
sandwich_koenker_l2 = sandwich_koenker_l2[order(sandwich_koenker_l2$Eps)]
colnames(sandwich_koenker_l2) = c("Eps", "Quantile", "Value")
sandwich_koenker_l2$Eps = as.factor(sandwich_koenker_l2$Eps)
sandwich_koenker_l2$Quantile = as.numeric(as.character(sandwich_koenker_l2$Quantile))

png(paste("plot/", filename, "_distance_meanl2.png", sep = ""), 
    width = 600, height = 400)
ggplot(sandwich_koenker_l2, aes(x = Quantile, log(Value))) +
  geom_line(aes(group = Eps, color = Eps, linetype = Eps),  size = 1) +
  scale_x_continuous(breaks= tau) +
  ggtitle("L2 Distance to the Truth by Epsilon Allocation to the Median") +
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation") + 
  ylim(c(1.25, 3))
dev.off()


colnames(sd_l2) = tau
sd_l2 = as.data.table(sd_l2)
sd_l2$Eps = allocate_eps
sd_l2_m = melt(sd_l2, id.vars = "Eps")
colnames(sd_l2_m) = c("Eps", "Quantile", "Sd")
sd_l2_m$Eps = as.factor(sd_l2_m$Eps)
sd_l2_m$Quantile = as.numeric(as.character(sd_l2_m$Quantile))
setkey(sandwich_koenker_l2, Eps, Quantile)
setkey(sd_l2_m, Eps, Quantile)
sandwich_koenker_l2 = merge(sandwich_koenker_l2, sd_l2_m, all.x = T)
sandwich_koenker_l2$Sd_n = sandwich_koenker_l2$Sd/sqrt(1000)

png(paste("plot/", filename, "_distance_l2se.png", sep = ""), 
    width = 600, height = 400)
ggplot(sandwich_koenker_l2, aes(x = Quantile, log(Value), ymin=log(Value-2*Sd_n), 
                             ymax=log(Value+2*Sd_n), fill = Eps)) + 
  facet_wrap(~Eps) +
  geom_line(aes(color = Eps), size = 1) +
  geom_ribbon(alpha = 0.2) +
  scale_x_continuous(breaks= tau[seq(1, 20, 3)]) +
  ggtitle("L2 Distance to the Truth by Epsilon Allocation to the Main Quantiles") + 
  ylab("log(Mean Distance)") +
  labs(color = "Epsilon \nAllocation", linetype = "Epsilon \nAllocation") +
  ylim(c(1.25, 3.25))
dev.off()



