#generate stable SEPP data

gen.seep.data= function(intercept,M,A,threshold,n,vzero){

#X is the data matrix with horizontal axis being time
  if(length(vzero)==0){
    vthreshold=rpois(M, intercept)}else{
vthreshold=vzero
vthreshold[which(vzero>threshold)]=threshold
                                  }
  X=matrix(0,ncol=n,nrow=M)
  
X[,1]= rpois(M, lambda=exp(intercept+A%*%as.matrix(vthreshold)))

for ( t in 2:n){
  X.temp = X[,t-1]
  X.temp[which(X[,t-1]>threshold)]=threshold
  
    X[,t]=rpois(M,lambda= exp(intercept+A%*%X.temp))

  
}
return(X)}

#compuate matrix estimation on given interval(dataset)
matrix.estimation=function(data,intercept,M,lambda,threshold){
  estimate=matrix(0,nrow=M,ncol=M)
  n.temp=ncol(data)
  data.x=data
  data.x[which(data.x >threshold)]=threshold
  
  for ( m in 1:M){
    pen <- penalized(data[m,2:n.temp]/exp(intercept), 
                     penalized = t(data.x[,1:(n.temp-1)]), unpenalized = ~0,
                     lambda1 = lambda,lambda2=10,model=c("poisson"),trace=F,maxiter=500)
   estimate[m,]= 
      coefficients(pen, "all")
     
    
  }
  return(estimate)
}

#compute residuals


residual.fun=function(data,data.x,intercept,m,B,n.temp){

  return( sum( sapply(2:n.temp, function(x)   exp(intercept+B[m,]%*%data.x[,x-1])-
                        data[m,x]*(intercept+B[m,]%*%data.x[,x-1])
  ))
    
    
  )
}


compute.residuals=function(data,intercept,M,lambda,B,threshold){
  data.x=data
  data.x[which(data.x >threshold)]=threshold
  n.temp=ncol(data)-1
  
  return(sum(sapply(1:M, function(m) residual.fun (data,data.x,intercept,m,B,n.temp))))
}

compute.residuals.from.data=function(data.temp,intercept,M,lambda,threshold ,delta,delta2){
  if( ncol(as.matrix(data.temp))>delta-1 && ncol(as.matrix(data.temp))<delta2+1){
     B.temp= matrix.estimation(data.temp,intercept,M,lambda,threshold)
     
     result=compute.residuals(data.temp,intercept,M,lambda,B.temp,threshold)}else{
       result=Inf
     }
    
    return(result)
}



#### dp main function
inner.function.dp=function(r,l,gridsize,gam,Best.value,data,intercept,M,lambda,threshold,delta,delta2){
 
 result=Inf
  if (l%% gridsize==1 ){
    cc= ifelse(l==1, -gam, Best.value[l-1] )
    result=  cc+gam+ 
      compute.residuals.from.data(data[,l:r],intercept,M,lambda/(r-l+1)^(1/2),threshold,delta,delta2)
  }
 return(result)
  }



