library(data.table)
dat = cbind(sample(0:1, 10, replace = T),
            sample(1:10, 10, replace = T),
            sample(1:10, 10, replace = T))
colnames(dat) = c("mu", "emp_1977", "pay_1976")
dat = as.data.table(dat)


product_col = function(data, col1, col2){
  data = as.data.table(data)
  
  if (length(col1) == 0 | length(col2) == 0){
    return(data)
  }
  vars = paste(col1, col2, sep = "*")
  data[, (vars) := lapply(.SD, function(x) x*data[[col1]]), .SDcols = col2]
  
  return(data)
}

dat1 = product_col(dat, col1, col2)

DT <- data.table(id=1:1000,year=round(runif(1000)*10), 
                 inc1 = runif(1000), inc2 = runif(1000), inc3 = runif(1000),    
                 deflator = rnorm(1000))
inc_cols = c('inc1', 'inc2','inc3')
DT[, (inc_cols) := lapply(.SD, function(x) 
  x * DT[['deflator']] ), .SDcols = inc_cols]