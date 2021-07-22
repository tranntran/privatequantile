library(data.table)
ut_l = NULL
ut_li =NULL
for (i in 100:104){
  filename = paste0('./output/simulations_seed', i,'_e1_10reps.Rdata')
  load(filename)
  ut_l = rbind(ut_l, ut_logit)
  ut_li = rbind(ut_l, ut_logit_inter)
}

tab1 = apply(ut_logit, 1, quantile, na.rm = T)
colnames(tab1) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
tab1
rownames(tab1) = c("Min", "Q1", "Median", "Q3", "Max")
t(round(tab1, 4))

#
library(gridExtra)
png(paste("./plot/utility_seed100104.png", sep =""), height = 50*nrow(tab1),
    width = 100*ncol(tab1))
grid.table(t(round(tab1, 5)))
dev.off()

tab2 = apply(ut_logit_inter, 1, quantile)
colnames(tab2) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
rownames(tab2) = c("Min", "Q1", "Median", "Q3", "Max")
t(round(tab2,4))
png(paste("./plot/utility_inter_seed100104.png", sep =""), height = 50*nrow(tab2),
    width = 100*ncol(tab2))
grid.table(t(round(tab2, 5)))
dev.off()

tab3 = rbind(apply(ut_logit, 1, mean), apply(ut_logit_inter, 1, mean))
rownames(tab3) = c("Mean Utility", "Mean Utility w/ Interaction")
colnames(tab3) = c("Original", "Stepwise-Fixed Slope", "Stepwise-Varying Slope",
                   "Sandwich-Fixed Slope", "Sandwich-Varying Slope", "Non-Private")
t(round(tab3, 5))
png(paste("./plot/utility_mean_seed100104.png", sep =""), height = 30*nrow(tab3),
    width = 150*ncol(tab3))
grid.table(round(tab3, 5))
dev.off()