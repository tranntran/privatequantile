# This code uses bootstrap on synthetic data to find the median and confidence
# intervals of specific utility metrics (such as gross employment level)
# 
# Main code idea: Each method has 5 different synthetic datasets
# For each method:
# 1. Find the actual median gross employment:
#   - Summing up the employment of one year for synthetic dataset 1, 2,...5
#   - Find the median gross employment for that year
#   - Repeat for 25 years
# 2. Find the confidence interval for median gross employment:
#   - Sample 5 new synthetic data version (1', 2',..., 5') from 
#     synthetic dataset 1, 2,...5
#   - Summing up the employment of one year for synthetic dataset 1', 2',...5'
#   - Find the median gross employment for that year
#   - Repeat for 25 years
#   After repeat the 4 steps above for 100 times, we have 100 median gross emp
#   and we can use this to get the confidence interval.
# We then repeat this for all the methods

tab_plot_median = NULL
for (k in 1:6) {
  n = 10
  data = list()
  for (i in 1:5) {
    data[[i]] = as.data.frame(matrix(sample(c(1:100, NA), 25*n, replace = TRUE), ncol = 25))
    colnames(data[[i]]) = paste("emp", c(1976:2000), sep = "_")
  }
  
  tab = apply(sapply(data, function(x) apply(x, 2, sum, na.rm = TRUE), simplify = TRUE), 1, median)
  
  reps = 100
  gross_emp_matrix = NULL
  
  median_gross_emp = matrix(NA, nrow = reps, ncol = 25)  
  for (i in 1:reps) {
    gross_emp_matrix = matrix(NA, nrow = 5, ncol = 25)  
    for (j in 1:5) {
      #read in data here
      samp_row = sample(c(1:n), n, replace = TRUE)
      samp_data = data[[j]][samp_row, ]
      gross_emp_matrix[j, ] = apply(samp_data, 2, sum, na.rm = TRUE)
    }
    median_gross_emp[i, ] = apply(gross_emp_matrix, 2, median)
  }
  confidence_int = apply(median_gross_emp, 2, quantile, c(0.025, 0.975))
  tab = as.data.frame(t(rbind(tab, confidence_int)))
  tab$year = as.numeric(substring(row.names(tab), 5, 8))
  tab$method = paste("Method", k)
  rownames(tab) = NULL
  colnames(tab) = c("value", "lower", "upper", "year", "method")
  tab_plot_median = rbind(tab_plot_median, tab)
}


# This code is similar as above but uses the mean instead of the median
tab_plot_mean = NULL
for (k in 1:6) {
  n = 10
  data = list()
  for (i in 1:5) {
    data[[i]] = as.data.frame(matrix(sample(c(1:100, NA), 25*n, replace = TRUE), ncol = 25))
    colnames(data[[i]]) = paste("emp", c(1976:2000), sep = "_")
  }
  
  tab = apply(sapply(data, function(x) apply(x, 2, sum, na.rm = TRUE), simplify = TRUE), 1, mean)
  
  reps = 100
  gross_emp_matrix = NULL
  
  mean_gross_emp = matrix(NA, nrow = reps, ncol = 25)  
  for (i in 1:reps) {
    gross_emp_matrix = matrix(NA, nrow = 5, ncol = 25)  
    for (j in 1:5) {
      #read in data here
      samp_row = sample(c(1:n), n, replace = TRUE)
      samp_data = data[[j]][samp_row, ]
      gross_emp_matrix[j, ] = apply(samp_data, 2, sum, na.rm = TRUE)
    }
    mean_gross_emp[i, ] = apply(gross_emp_matrix, 2, mean)
  }
  confidence_int = apply(mean_gross_emp, 2, quantile, c(0.025, 0.975))
  tab = as.data.frame(t(rbind(tab, confidence_int)))
  tab$year = as.numeric(substring(row.names(tab), 5, 8))
  tab$method = paste("Method", k)
  rownames(tab) = NULL
  colnames(tab) = c("value", "lower", "upper", "year", "method")
  tab_plot_mean = rbind(tab_plot_mean, tab)
}






















# bootstrap.median <- function(data, B = 1000, confidence = 0.95) {
#   n <- length(data)
#   median <- rep(0, B)
#   for (i in 1:B) {
#     median[i] <- median(sample(data, size = n, replace = T))
#   }
#   med.obs <- median(median)
#   c.l <- round((1 - confidence) / 2 * B, 0)
#   c.u <- round(B - (1 - confidence) / 2 * B, 0)
#   l <- sort(median)[c.l]
#   u <- sort(median)[c.u]
#   cat(c.l / 1000 * 100, "-percentile:      ", l, "\n")
#   cat("Median: ", med.obs, "\n")
#   cat(c.u / 1000 * 100, "-percentile:      ", u, "\n")
#   return(median)
# }