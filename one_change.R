#SEPP two change  code
rm(list=ls())
set.seed(50)
#setwd("/Users/darenw/Dropbox/SEPP/code/")
setwd("/Users/darenwang/Dropbox/SEPP/code/scene1")
source("functions.R")
library(penalized)
#M is the dimension
M=30
#n is the spacing between change point
n=300
#lambda is the lasso parameter
lambda=1500
#thresholding makes the process stable
threshold=6
#intercept of the model. Assume to be known as in the existing literature
intercept=1/2

#change point estimation 
estimate.change=matrix(0,nrow=100,ncol=6)

for ( bb in 1:6){

  
#define two   
v1=(2*(seq(1,30,1)%%2)-1)*(0.05+0.05*bb)
v2=-v1
AA=matrix(0,nrow = 30,ncol=28)

A=cbind(v1,v2,AA)
A2=cbind(v2,v1,AA)

#50 experiments with
for ( aa in 1:50){
  print(aa)
#when n is super large, it works fine


data1=gen.seep.data(intercept,M,A,threshold,n,vzero=c())
data2=gen.seep.data(intercept,M,A2,threshold,n,vzero=data1[,n])

data=cbind(data1,data2)
rec.rr=rep(0,150)

for ( j in 1:150 ){
  i=2*j+150
  #print(i)
  data.temp.1=data[,1:i]
  data.temp.2=data[,(i+1):(2*n)]
  lambda1= lambda/(i)^(1/2)
  lambda2= lambda/(2*n-i)^(1/2)
  
  B1=matrix.estimation(data.temp.1,intercept,M,lambda=lambda1,threshold)
  rr1=compute.residuals(data.temp.1,intercept,M,lambda=lambda1,B1,threshold)
  B2=matrix.estimation(data.temp.2,intercept,M,lambda=lambda2,threshold)
  rr2=compute.residuals(data.temp.2,intercept,M,lambda=lambda2,B2,threshold)
  rec.rr[j]=rr1+rr2
   #print(rec.rr[j])
}
 #rec.rr1=rec.rr
 #data=data.1
#plot(rec.rr)
estimate.change[aa,bb]=which.min(rec.rr)*2+150
print(estimate.change[aa,bb])
}

}
colMeans(abs(estimate.change-300   )/600)
