rm(list=ls())
#experiment for n=100

setwd("/Users/darenwang/Dropbox/SEPP/code/scene2/")
library(penalized)
source("functions.R")
n=100

set.seed(n)
#lambda is lasso penalization
lambda=500
#M is dimension 
M=50
#s is sparsity
s=30
#penalize size of partition
gam=n/50

gridsize=5
#grid size =5 is expensive but a lot more accurate

#large factor gives exact recovery
factor=0.12


delta=n/2-1
delta2=1.5*n



intercept=1/2
threshold=4
#define transition matrices
A1=matrix(0,M,M)
diag(A1[,-1])=1
diag(A1)=1
diag(A1[-1,])=-1
A1=A1*factor
A1[(s+1):M,(s+1):M]=0

A2=matrix(0,M,M)
diag(A2[,-1])=1
diag(A2)=-1
diag(A2[-1,])=1
A2=A2*factor
A2[(s+1):M,(s+1):M]=0

A3=matrix(0,M,M)
diag(A3[,-1])=1
diag(A3)=1
diag(A3[-1,])=-1
A3=A3*factor
A3[(s+1):M,(s+1):M]=0

# end of transition matrix definition

#rec is the record of 50 repition
rec=vector("list",length=50)

for ( rr in 1:50){
  
  print(rr)
data1=gen.seep.data(intercept,M,A1,threshold,n,vzero=c())
data2=gen.seep.data(intercept,M,A2,threshold,n,vzero=data1[,n])
data3=gen.seep.data(intercept,M,A3,threshold,n,vzero=data2[,n])

data=cbind(data1,data2,data3)



 
#########end of generate data


Best.value=rep(0, 3*n)
partition=rep(0,3*n)
start_time <- Sys.time()


#start of DP
for( r in 1:(3*n)){
  Best.value[r]=Inf
  if(r %%gridsize==0){
    #print(r)
    b=sapply(1:r,function(l) 
      inner.function.dp(r,l,gridsize,gam,Best.value,data,intercept,M,lambda,threshold,delta,delta2)     
            )
        if (min(b) <Best.value[r]){Best.value[r]=min(b)
        partition[r]=which.min(b)-1} 
      }
      
    }
end_time <- Sys.time()
 print(paste("time=",end_time-start_time,"mins"))


#end of DP


#pull out the change points
cc=3*n
estimate.change=cc
while(cc[1]!=0){
  estimate.change=c(partition[cc],estimate.change)
  cc=partition[cc]}
print(estimate.change)

rec[[rr]]=estimate.change
}


save.image(file="n=100.RData")


dist=c()
for ( rr in 1:50){
  dist=c(dist, hausdorff.distance(rec[[rr]][c(-1,-length(rec[[rr]]))], c(n,2*n)))
  
  
}
mean(dist/450)
