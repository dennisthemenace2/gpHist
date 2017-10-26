##estimate hyperparameters...

downhillsimplex = function(fc,nParams,bp=NULL,lower=0,upper=1,it=1000,tol=1e-10, control= list(alpha=1,gamma=2, beta=0.5,sigma = 0.5) ){
  
  # Make sure GP is right.
  if(! is.function(fc) ) {
    print("downhillsimplex(): fc has to be of type function")
    return(NaN);
  }
  if(nParams<=0){
    print("downhillsimplex(): nParams has to be greater than 1")
    return(NaN);
  }
  
  if(length(lower)< nParams){
    diff = nParams -length(lower)
    lower = c(lower,rep(lower[length(lower)],diff))
  }
  
  if(length(upper)< nParams){
    diff = nParams -length(upper)
    upper = c(upper,rep(upper[length(upper)],diff))
  }
  
  if(is.null(bp)){
    #bp = matrix(runif( (N+1)*N,0,5),nrow=N+1,ncol=N )
    bp = matrix(0,nrow=nParams,ncol=nParams+1 )
    for(i in 1:nParams){
      bp[i,] = runif( (nParams+1),lower[i],upper[i])
    }
  }else{
    if(ncol(bp) != nParams+1){
      print('solution points wrong dimension! your input will be replaced')
      bp = matrix(runif( (nParams+1)*nParams,0,5),nrow=nParams+1,ncol=nParams )
    }
    if(nrow(bp)>nParams){
      print('solution points dosnt need function value.')
      bp = bp[,1:nParams]
    }
  }
  
  alpha  = control$alpha #1 # =1 (Reflexion)
  gamma =control$gamma #2 #= 2 (Expansion)
  beta = control$beta#0.5 #1/2 (Kontraktion)
  sigma =control$sigma # 0.5 #Komprimierung
  
  
  # void Cdownhillsimplex(SEXP func,SEXP sys,int* nParams,double* bp,double* lower,double* upper,int *resPtr, double alpha= 1,double gamma=2,double beta=0.5,double sigma=0.5,unsigned int it=1000,double tol=1e-10){
  resPointer =1;
  ##  print(bp)
  output =.C("Cdownhillsimplex",
             func = fc,
             frame = sys.frame(),
             nParams  = as.integer(nParams),
             result = as.double(bp),
             lower = as.double(lower),
             upper = as.double(upper),
             resPointer = as.integer(resPointer),
             alpha=as.double(alpha),
             gamma= as.double(gamma),
             beta=as.double(beta),
             sigma = as.double(sigma),
             it=as.integer(it),
             tol = as.double(tol)
             
  )
  start = (output$resPointer-1)*nParams+1
  end  = start +nParams-1;
#  print(output$result)
  res = matrix(output$result[start:end],nrow=1);
  return(res);
}
