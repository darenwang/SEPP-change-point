rm(list=ls())
#experiment for n=100

setwd("/Users/darenwang/Dropbox/SEPP/code/scene3/")
library(penalized)
source("functions.R")
n=150
set.seed(50)


#iteration is the number of iterated experiments
iteration=50

#penalize size of partition

gridsize=10
#grid size =5 is expensive but a lot more accurate


delta=50
#the lower bound of the estimated change point spacing: if delta is too small, then restricted eigenvalue condition breaks
#and numerically unstable
delta2=2*n
#maximum of the estimated change point spacing: reduce computation cost


intercept=1/5
threshold=4
#define transition matrices


#rec is the record of 50 repition
rec=vector("list",length=50)
rec.list=vector("list",length=4)
M.seq=c(30,40,50,60)

for ( bb in 1:4){
  print(paste("bb=",bb))
  M=M.seq[bb]
  lambda=600
  #lambda= 200*log(M)
  gam=15*log(M)
  #gam=15*log(M)
  rec.list[[bb]]=vector("list", length=iteration)
  
  A0=matrix(0, nrow=M, ncol=M-3)
  v1=c(c(-3,4.5, 5,-3),rep(0,M-4))
  v1=v1/(2*norm(v1,type="2"))
  #A1=v1%*%t(v1)
  
  v2=c(c(2.5,-1.5, -0.5,1.5),rep(0,M-4))
  v2=v2/(2*norm(v2,type="2"))
  #A2=v2%*%t(v2)
  sum((v1-v2)^2)
  v3=c(c(-3,2, -5,-3),rep(0,M-4))
  v3=v3/(2*norm(v3,type="2"))
  #A3=v3%*%t(v3)
  sum((v2-v3)^2)
  A1=cbind(v1,v2,v3,A0)
  A2=cbind(v2,v3,v1,A0)
  
  A3=cbind(v3, v2,v1,A0)
  # end of transition matrix definition
for ( rr in 1:iteration){
  print(paste("rr=",rr))
  

data1=gen.seep.data(intercept,M,A1,threshold,n,vzero=c())
data2=gen.seep.data(intercept,M,A2,threshold,n,vzero=data1[,n])
data3=gen.seep.data(intercept,M,A3,threshold,n,vzero=data2[,n])

data=cbind(data1,data2,data3)



# 
# compute.residuals.from.data(data[,1:n],intercept,M,lambda/n^(0.5),threshold ,delta,delta2=3*n)+
#   compute.residuals.from.data(data[,1:n+n],intercept,M,lambda/n^(0.5),threshold ,delta,delta2=3*n)-
#   compute.residuals.from.data( cbind(data1,data2),intercept,M,lambda/(2*n)^(0.5),threshold ,delta,delta2=3*n)
# 
# 
# compute.residuals.from.data(data[,1:n+n],intercept,M,lambda/n^(0.5),threshold ,delta,delta2=3*n)+
#   compute.residuals.from.data(data[,(2*n+1):(3*n)],intercept,M,lambda/n^(0.5),threshold ,delta,delta2=3*n)-
#   compute.residuals.from.data(data[,(n+1):(3*n)],intercept,M,lambda/(2*n)^(0.5),threshold ,delta,delta2=3*n)
# 
# 
# compute.residuals.from.data(data[,1:(n/2)],intercept,M,lambda/(n/2)^(0.5),threshold ,delta,delta2=3*n)+
#   compute.residuals.from.data(data[,(n/2+1):(n)],intercept,M,lambda/(n/2)^(0.5),threshold ,delta,delta2=3*n)-
#   compute.residuals.from.data(data[,1:n],intercept,M,lambda/(n)^(0.5),threshold ,delta,delta2=3*n)
# 
# 
# B=matrix.estimation(data[,1:n],intercept,M,lambda/n^(0.5),threshold)
# norm(B-A1,type="2")
# norm(A1,type="2")
 
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

rec.list[[bb]][[rr]]=estimate.change
}

}
save.image(file="M.RData")


dist1=matrix(0,nrow=4,ncol=50)
for( bb in 1:4){
for ( rr in 1:50){
  dist1[bb,rr] =hausdorff.distance(rec.list[[bb]][[rr]][c(-1,-length(rec.list[[bb]][[rr]]))], c(n,2*n))
  
  
}
}
rowMeans(dist1/(3*n))





