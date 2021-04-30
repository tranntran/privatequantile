
#library(mcmc)
library(MASS)
#library(gtools)
#library(rmutil)
library(tictoc)

set.seed(1000)

#n_vect=c(100,1000,10000,100000)#number of observations
#n_vect = c(1e1,1e2,1e3,1e4,1e5)#,1e6)#,1e7)
n_vect = 1e5
ep = 0.01
tau=1/2
m=1# number of predictors
truth = seq(-1,1-2/m,by=2/m)
beta = c(0,truth)
B=1# bound on L1 norm of beta (objective perturbation)
num_reps=100 #number of simulations we will average per  epsilon value
runs=1000#   for metrops.








norm2 = function(v){
    return(sqrt(sum(v^2)))
}

norm1 = function(v){
    return(sum(abs(v)))
}

DistanceToBeta = function(guess, truth){
  return(sqrt(sum( (guess-truth)^2)))
}



metrop1 = function(logD,init,nbatch,scale,blen){
###   blen doesnt do anything.
    dim = length(init)
    out = list(accept = 0,
    batch = matrix(rep(0,nbatch*dim),nrow=nbatch,ncol=dim))
    U = matrix(runif(nbatch*dim,min=0,max=1),nrow=nbatch,ncol=dim)
    Prop = matrix(rnorm(nbatch*dim,m=0,s=scale),nrow=nbatch,ncol=dim)

    for(r in 1:nbatch){
        if(r==1){
            out$batch[1,] =init
            oldLogD = logD(init)
        }
        else
            out$batch[r,] = out$batch[r-1,]
        ###
        for(i in 1:dim){
            oldVector = out$batch[r,]
           
            newVector = oldVector
            newVector[i] = newVector[i] + Prop[r,i]
            newLogD = logD(newVector)

            TestValue = exp((newLogD - oldLogD))
            if(U[r,i]<TestValue){
                out$batch[r,] = newVector
                out$accept = out$accept + 1/(nbatch*dim)
                oldLogD = newLogD
            }
        }
    }
    return(out)
}


getScale = function(logA){
 scale = .1
 prevBelow=0
 prevAbove = Inf
 count=0

 init = rep(0,m+1)
 #init = beta_hat
 
 out = metrop1(logA,init=init,nbatch=1000,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
#out$accept
 # out = metrop(logD,init=t(tail(out$batch,n=1)),nbatch=1000,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
 while((out$accept<.1 | out$accept>.2)&(count<=10)){
   if(out$accept<.1){
    prevAbove = scale 
   }
   else if(out$accept>.2){
    prevBelow = scale 
   }
   
   if(prevAbove<Inf){
    scale = (1/2)*(prevBelow+prevAbove) 
   }
   else{
    scale = 2*scale 
   }
   
   out = metrop1(logA,init=init,nbatch=1000,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
   count = count+1
 }
 if(count==11){
  print('scale did not converge')
   return(0)
 }
 return(scale)
}



KNG = function(ep,tau,sumX,X,Y,scale){

    init = rep(0,m+1)

    logA = function(beta){
        #if(norm1(beta)<=B){
        left = cbind(Y,X)%*% c(1,-beta)
        lessEq = (left <=0)
        return(-(ep/2) *
               max(abs(-tau*sumX + t(X)%*%lessEq))/
               ((1-tau)*2*1)+
               (-1/(2))*(beta%*%beta)###   N(0,1) is base measure.
               )
         #}
        #else{
        #    return(-Inf) 
        #}
    }
                                        #figure out a good step size for Metropolic hastings for Exponential Mechanism
                                        #if(rep==1){
    #if(scale==0)
    #scale = 0#getScale(logA=logA)
    while(scale==0){
        scale=getScale(logA=logA)
    }
                                        #}
                                        #
    out = metrop1(logA,init=init,nbatch=runs,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
    beta_kng =t(tail(out$batch,n=1))
    return(list(beta_kng,scale))
}





ExpMech = function(ep,tau,X,Y,scale){
    de = 2*max(tau,1-tau)*(1+B)# delta for exponential mechanism
    init=rep(0,m+1)
                                        #log-likelihood for the exponential mechanism
    logA = function(beta){
        
        if(norm1(beta)<=B){
            XB = X%*%beta
            lessEq = (Y<=XB)
            grEq = !lessEq
            weight = (tau-1)*lessEq  +  tau*grEq
            
            
            return(-ep/(2*de)*(t(Y-XB)%*%weight))
        }
        else{
            return(-Inf) 
        }
    }
    while(scale==0)
        scale = getScale(logA)
    
    out = metrop1(logA,init=init,nbatch=runs,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
    beta_exp =t(tail(out$batch,n=1))
    return(list(beta_exp,scale))
}


QuantileRegression = function(tau,sumX,X,Y){

    obj = function(beta){
        XB = X%*%beta
        lessEq = (Y<=XB)
        grEq = !lessEq
        weight = (tau-1)*lessEq  +  tau*grEq
        return((t(Y-XB)%*%weight))
    }
    grad = function(beta){
        XB = X%*%beta
        lessEq = (Y<=XB)
        return(-tau*sumX + t(X)%*%lessEq)
    }
    
    min = optim(par = rep(0,m+1), fn=obj,gr=grad,method = "L-BFGS-B")
    return(min$par)
}





############################################################
############################################################
###   SIMULATIONS
############################################################
############################################################



RegressionVect =matrix(rep(0,length(n_vect)*num_reps),nrow=length(n_vect),ncol=num_reps)
KngVect = RegressionVect
ExpVect = RegressionVect
#ObjPertVect = RegressionVect

tic("total")
for(index in 1:length(n_vect)){
    n=n_vect[index]
    print(n)
    scale_exp=0
    scale_kng=0
    tic("loop")
    for(rep in 1:num_reps){
        print(rep)
        uni = runif(n=n*m,min=-1,max=1)
        X = matrix(uni, nrow=n,ncol=m)# all entries between -1 and 1

        W = rnorm(n,mean=0, sd=1)
        Y = X%*%truth + W# the response
                                        # largest magnitude of y
        R = max(abs(Y))# in the future, should DP this.
        Y = Y/R #   make between -1 and 1.
        
        X = cbind(rep(1,n),X)
        sumX = apply(X=X,2,FUN=sum)
        #XtX = t(X)%*%X
        #XtY =t(X)%*%Y
        #YtY=t(Y)%*%Y

        #beta_hat = QuantileRegression(tau=tau,sumX=sumX,X=X,Y=Y)
        #out_exp = ExpMech(ep=ep,tau=tau,X=X,Y=Y,scale=scale_exp)
        #scale_exp = out_exp[[2]]
        #beta_exp = out_exp[[1]]

        out_kng = KNG(ep=ep,tau=tau,sumX=sumX,X=X,Y=Y,scale=scale_kng)
        scale_kng = out_kng[[2]]
        beta_kng= out_kng[[1]]

        #beta_obj = ObjPert(lambda=2*(m+1),xi=4*(1+B),ep=ep,XtX=XtX,XtY=XtY,YtY)


        #RegressionVect[index,rep] =DistanceToBeta(beta,R*beta_hat)
        #ExpVect[index,rep] = DistanceToBeta(beta,R*beta_exp)
        KngVect[index,rep] = DistanceToBeta(beta,R*beta_kng)
        #ObjPertVect[index,rep] = DistanceToBeta(beta,R*beta_obj)
    }
    toc()
}
toc()

RegressionMean = apply(X=RegressionVect,1,FUN=mean)
ExpMean        = apply(X=ExpVect,1,FUN=mean)
#ObjPertMean    = apply(X=ObjPertVect,1,FUN=mean)
KngMean        = apply(X=KngVect,1,FUN=mean)

RegressionSd   = apply(X=RegressionVect,1,FUN=sd)
ExpSd          = apply(X=ExpVect,1,FUN=sd)
#ObjPertSd      = apply(X=ObjPertVect,1,FUN=sd)
KngSd          = apply(X=KngVect,1,FUN=sd)


KngDist_1 = apply(KngVect, 1, quantile)
KngSd_1 = apply(X=KngVect,1,FUN=sd)

KngDist_0.1 = apply(KngVect, 1, quantile)
KngSd_0.1 = apply(X=KngVect,1,FUN=sd)

KngDist_0.01 = apply(KngVect, 1, quantile)
KngSd_0.01 = apply(X=KngVect,1,FUN=sd)

thick=4
font=1.5
pdf("KNGplotQm1.pdf",width=7,height=7)
par(mar=c(5.1,5.1,4.1,2.1))
plot(log(RegressionMean)/log(10),type="l",ylim = c(-2.5,.5),xaxt="n",xlab="n",lwd=thick,ylab=expression("log_10( Average L2 Distance to Truth )"),cex.lab=font,cex.axis=font,cex.sub=font)
lines(log(ExpMean)/log(10),col="red",lty=2,lwd=thick)
#lines(log(ObjPertMean)/log(10),col="blue",lty=3,lwd=thick)
lines(log(KngMean)/log(10),col="orange",lty=4,lwd=thick)
#abline(a=0,b=0)
axis(side=1,1:length(n_vect),n_vect,cex.axis=font)
legend("topright",c("NonPrivate","ExpMech","KNG"),
       col=c("black","red","orange"),
       lty=c(1,2,4),lwd=thick,cex=font)
dev.off()






summary(
c(log(KngMean+KngSd/sqrt(num_reps))/log(10)-log(KngMean-KngSd/sqrt(num_reps))/log(10),
  #
#log(ObjPertMean+ObjPertSd/sqrt(100))/log(10)-log(ObjPertMean-ObjPertSd/sqrt(100))/log(10),
#
log(RegressionMean+RegressionSd/sqrt(num_reps))/log(10)-log(RegressionMean-RegressionSd/sqrt(num_reps))/log(10),
#
log(ExpMean+ExpSd/sqrt(num_reps))/log(10)-log(ExpMean-ExpSd/sqrt(num_reps))/log(10)#
))
 
