rm(list = ls())
n<-1000
beta1<-4
beta2<-5
eps<-0.1
tau<-0.5
lambda=0.01

x<-rexp(n,rate=lambda)
#x <- runif(n, -100, 100)
y<-beta1+beta2*x+rexp(n,rate=lambda)
#y<-beta1+beta2*x+rnorm(n)

#x = log(x)
#y = y/max(abs(y))
maxx<-max(x)
meanx<-mean(x)

loss<-function(b){
  #exp(
    -(eps*n / (4 * (1-tau)*maxx)) *norm(
      -tau * cbind(1, meanx) + sum(t(cbind(1, x))%*%((cbind(1,x)%*%b > y)))/n, type = "M") #-0.00001*t(b)%*%b
 #)
}

grid_n<-100
b_grid1<-seq(0,100,length=grid_n)
b_grid2<-seq(0,10,length=grid_n)

l_val<-matrix(nrow=grid_n,ncol=grid_n)

for(i in 1:grid_n){
  for(j in 1:grid_n){
    l_val[i,j] = loss(c(b_grid1[i],b_grid2[j]))
  }
}

image(b_grid1,b_grid2,l_val)

min(l_val)/max(l_val)
quantreg::rq(y~x, tau = 0.5)


