##estimate hyperparameters...

##make sure datatransform function does preserves orders of X's, this way the order has only be computed once.
estimateHyperParameters=function(X, Y, paramLower=0,paramUpper=1,datatransform=NULL,nParams=0,it = 100,tol=0.001){

  GP = NULL  
  dhfc = function(p){
    P = length(p)
    
    if(!is.null(datatransform) ){
      X_trans = datatransform(X,p[2:P])
      if(is.null(GP)){
        ret = gpHist(X_trans ,Y,p[1]);
      }else{
     #   ret = gpHist(X_trans ,Y,p[1],orders=GP$orders ,alphas = GP$alpha);
        ret = gpHist(X_trans ,Y,p[1],orders=GP$orders);
      }
    }else{
      if(is.null(GP)){
        ret = gpHist(X ,Y,p[1]);
      }else{
        ret = gpHist(X ,Y,p[1],orders=GP$orders);
      #  ret = gpHist(X ,Y,p[1],orders=GP$orders,alphas = GP$alpha);
      }
    }
    GP<<-ret

    ret$logmarginal
  }
  resmat = downhillsimplex(dhfc,nParams+1,lower=paramLower,upper=paramUpper,it = it,tol=tol)
  
}


estimateHyperParametersZZZ=function(X, Y, paramLower=0,paramUpper=1,shift =100,it = 100,tol=0.001){
  # Make sure GP is right.
  
  nParams = ncol(X) +1
  
  if(length(lower)< nParams){
    diff = nParams -length(lower)
    lower = c(lower,rep(lower[length(lower)],diff))
  }
  
  if(length(upper)< nParams){
    diff = nParams -length(upper)
    upper = c(upper,rep(upper[length(upper)],diff))
  }
  
  bp = matrix(0,nrow=nParams,ncol=nParams+1 )
  for(i in 1:nParams){
    bp[i,] = runif( (nParams+1),lower[i],upper[i])
  }

  
  
  # Make sure A and B are matrices.
  if ( ( ! is.matrix(X) ) || ( ! is.matrix(Y) ) ) {
    print("gpHist(): input must be matrices!")
    return(NaN);
  }
  
  orders = matrix(,ncol=ncol(X),nrow=nrow(X))
  for(i in 1:ncol(X)){
    orders[,i  ]= order(X[,i])   
  }

  alpha  = 1 # =1 (Reflexion)
  gamma = 2 #= 2 (Expansion)
  beta = 0.5 #1/2 (Kontraktion)
  sigma = 0.5 #Komprimierung
  
  
  # void Cdownhillsimplex(SEXP func,SEXP sys,int* nParams,double* bp,double* lower,double* upper,int *resPtr, double alpha= 1,double gamma=2,double beta=0.5,double sigma=0.5,unsigned int it=1000,double tol=1e-10){
  resPointer =1;
  ##  print(bp)
  output =.C("CdownhillsimplexIN",
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
             tol = as.double(tol),
             X = as.double(X),
             Y = as.double(Y),
             orders =as.double(orders),
             nrow = nrow(X),
             ncol= ncol(X)
             
  )
  start = (output$resPointer-1)*nParams+1
  end  = start +nParams-1;
  #  print(output$result)
  res = matrix(output$result[start:end],nrow=1);
  print(res)
  return(res);
}