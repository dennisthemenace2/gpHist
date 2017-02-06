##estimate hyperparameters...


Rdownhillsimplex = function(fc,N,bp=NULL,lower=0,upper=1,it=1000,tol=1e-10){
  
  
  chkBound=function(p,lower,upper){
    toolow =p<lower
    if(any(toolow) ){
      p[toolow] =lower[toolow]
    }
    toohigh =p>upper
    if(any(toohigh) ){
      p[toohigh] = upper[toohigh]
    }
    p
  }
  #######
  alpha  = 1 # =1 (Reflexion)
  gamma = 2 #= 2 (Expansion)
  beta = 0.5 #1/2 (Kontraktion)
  sigma = 0.5 #Komprimierung
  
  if(length(lower)<N){
    lower = rep(lower,N)
  }
  if(length(upper)<N){
    upper = rep(upper,N)
  }
  
  if(is.null(bp)){
    #bp = matrix(runif( (N+1)*N,0,5),nrow=N+1,ncol=N )
    bp = matrix(0,nrow=N+1,ncol=N )
    for(i in 1:N){
      bp[,i] = runif( (N+1),lower[i],upper[i])
    }
  }else{
    if(nrow(bp) != N+1){
      print('solution points wrong dimension! your input will be replaced')
      bp = matrix(runif( (N+1)*N,0,5),nrow=N+1,ncol=N )
    }
    if(ncol(bp)>N){
      print('solution points dosnt need function value.')
      bp = bp[,1:N]
    }
  }
  
  res=c()
  for(k in 1:(N+1) ){
    res =c(res, fc(bp[k,]) )
  }
  
  bp = cbind(bp,res)
  
  for(i in 1:it){
    
    ##plot
    #   points(bp[1,],bp[2,],col='red')
    #
    
    ord = order(bp[,N+1])
    bp = bp[ ord,] ##order

  #  print('bp');
  #  print(bp);
    
    if(N>1){
      m = colMeans(bp[1:N,1:N])
    }else{
      m = sum(bp[,1])
    }
    
 #   print('\n\n\nm');
#    print(m);
    #  print(bp[1,N+1]-bp[N+1,N+1] )
    if(abs(bp[1,N+1]-bp[N+1,N+1]) < tol){
      cat(c('converged after:',i,' iterations. value:', abs(bp[1,N+1]-bp[N+1,N+1]) ,'\n' ) )
      break;
    }
    
    r =(1+alpha)*m - alpha*bp[N+1,1:N]
 #  print('querry r')
#    print(r)
    r = chkBound(r,lower,upper )
    
    np = fc(r)
    if(np<bp[1,N+1]){ ##better
      e =(1+gamma)*m -gamma*bp[N+1,1:N]
 #        print('querry e')
#         print(e)
      ###check bounderys
      e = chkBound(e,lower,upper )
      
      npe = fc(e)
      if(npe<np){ ##e better
        bp[N+1,1:N] = e
        bp[N+1,N+1] = npe
      }else{
        bp[N+1,1:N] = r
        bp[N+1,N+1] = np
      }
   #   print('bp')
  #    print(bp)
      next; ##jump back
    }
    
    if(np<bp[N,N+1]){
      bp[N+1,1:N] = r
      bp[N+1,N+1] = np
      
 #     print('bp')
#      print(bp)
      next;
    }
    
#    print('np')
#    print(np)
    if(np < bp[N+1,N+1]){ ##np better
  #    print('np is better')
      h = r
  #    print("np\n");
  #    print(h);
    }else{
    #  print('n+1 is better')
      h =  bp[N+1,1:N]
      np = bp[N+1,N+1]
    }
    cp = beta*m+(1-beta)*h
  #  print('querry cp before chkbound')
  #  print(cp)
    cp = chkBound(cp,lower,upper )
  #  print('querry cp')
  #  print(cp)
    nc = fc(cp)
    if(nc< bp[N+1,N+1]){ ##c is better
      bp[N+1,1:(N+1) ] = c(cp,nc)
#      print('c is better bp')
 #     print(bp)
      next;
    }
    for(k in 2:(N+1) ){ ##shrink
      bp[k,1:N] = bp[1,1:N] + sigma*(bp[k,1:N] - bp[1,1:N] )
      ##eval
      bp[k,1:N] = chkBound(bp[k,1:N],lower,upper )
      
      bp[k,N+1] = fc(bp[k,1:N])
      #   print('shrink')
      #    print(bp[k,1:N])
      #    print(bp[k,N+1])
    }
   
    
  #  print('bp last')
  #  print(bp)
     
  }
  #  cat(c('iterations:',i,'\n') )
  ###just return parameters
  ret = matrix(bp[1,1:(ncol(bp)-1)],nrow=1)
}

downhillsimplex = function(fc,nParams,bp=NULL,lower=0,upper=1,it=1000,tol=1e-10){
  
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
  
  alpha  = 1 # =1 (Reflexion)
  gamma = 2 #= 2 (Expansion)
  beta = 0.5 #1/2 (Kontraktion)
  sigma = 0.5 #Komprimierung
  
  
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