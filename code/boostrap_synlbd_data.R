# This code uses bootstrap on synthetic data to find the median and confidence
# intervals of specific utility metrics (such as gross employment level)

n = 10
data = matrix(sample(c(1:100, NA), 25*n, replace = TRUE), ncol = 25)
data = as.data.frame(data)
colnames(data) = paste("emp", c(1976:2000), sep = "_")


reps = 100
gross_emp_matrix = NULL

median_gross_emp = matrix(NA, nrow = reps, ncol = 25)  
for (i in 1:reps) {
  gross_emp_matrix = matrix(NA, nrow = 5, ncol = 25)  
  for (j in 1:5) {
    #read in data here
    samp_row = sample(c(1:n), n, replace = TRUE)
    samp_data = data[samp_row, ]
    gross_emp_matrix[j, ] = apply(samp_data, 2, sum, na.rm = TRUE)
  }
  median_gross_emp[i, ] = apply(gross_emp_matrix, 2, median)
}
apply(median_gross_emp, 2, quantile, c(0.025, 0.975))

bootstrap.median <- function(data, B = 1000, confidence = 0.95) {
  n <- length(data)
  median <- rep(0, B)
  for (i in 1:B) {
    median[i] <- median(sample(data, size = n, replace = T))
  }
  med.obs <- median(median)
  c.l <- round((1 - confidence) / 2 * B, 0)
  c.u <- round(B - (1 - confidence) / 2 * B, 0)
  l <- sort(median)[c.l]
  u <- sort(median)[c.u]
  cat(c.l / 1000 * 100, "-percentile:      ", l, "\n")
  cat("Median: ", med.obs, "\n")
  cat(c.u / 1000 * 100, "-percentile:      ", u, "\n")
  return(median)
}