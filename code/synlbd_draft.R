library(data.table)
dat = cbind(sample(0:1, 10, replace = T),
            sample(1:10, 10, replace = T),
            sample(1:10, 10, replace = T))
colnames(dat) = c("mu", "emp_1977", "pay_1976")
dat = as.data.table(dat)


product_col = function(data, col1, col2){
  data = as.data.table(dat)
  
  if (length(col1) == 0 | length(col2) == 0){
    return(data)
  }
}