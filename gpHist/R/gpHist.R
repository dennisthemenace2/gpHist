gpHist <-
function(X,Y,sigma,orders=NULL,alphas = NULL,k=1) {
  
    # Make sure A and B are matrices.
    if ( ( ! is.matrix(X) ) || ( ! is.matrix(Y) ) ) {
        print("gpHist(): input must be matrices!")
        return(NaN);
    }
  
    if(is.null(orders)){ ##create order but it would be better to have better strucutre
      orders = matrix(,ncol=ncol(X),nrow=nrow(X))
      for(i in 1:ncol(X)){
        orders[,i  ]= order(X[,i])   
      }
    }else{
      if( nrow(X)!= nrow(orders) || ncol(X) != ncol(orders) ){
        print("gpHist(): dimension of orders not correct!")
        return(NaN);
      }
    }
  
  nrowA      = nrow(X)
  if(is.null(alphas) ){
    multResult = matrix(rep(0,nrowA),nrow=nrowA,ncol=1)  
  }else{
    if(nrow(alphas)==nrow(X)){
      multResult = alphas;
    }else{
      print("gpHist(): alphas dimension not right!")
      return(NaN);
    }
  }
  if(nrow(Y)!= nrow(X)){
    print("gpHist(): number of examples of X and results Y are not matching!")
    return(NaN);
  }
  ##check k real quick
  if(k<=0){
    print("gpHist(): number of estimated eigenvalues can not be negative or zero!")
    return(NaN);
  }
  if(k>nrow(X) ){
    cat(c("gpHist(): maximal number of eigenvalues in K is ", nrow(X),'!'))
    return(NaN);
  }
  
   
   logmarginal = 0; 
   lambda = matrix(0,nrow=k, ncol=1); 
   eigenvec = matrix(0,nrow=nrow(X), ncol=k); 
   
     output =.C("CgpHist",
         result = as.double(multResult),
         mat1 = as.double(X),
         numRows  = as.integer(nrowA),
         numCols  = as.integer(ncol(X)),
         mat2 = as.double(Y),
         numRows2  = as.integer(nrow(Y)),
         numCols2  = as.integer(ncol(Y)),
         sigma = as.double(sigma),
         orders = as.double(orders),
         logmarginal= as.double(logmarginal),
         lambda = as.double(lambda),
         vector = as.double(eigenvec),
         k = as.integer(k)
     )
   
    ret = list('logmarginal'=output$logmarginal,'orders'=orders,'alpha'=matrix(output$result,nrow=nrowA),'sigma'=sigma,'lambda'=output$lambda,'vectors'=matrix(output$vector,nrow=nrow(X), ncol=k))

#    ret$newVec = matrix(output$vector,nrow=nrow(X), ncol=k);
#    ret$newlambda = output$lambda;
 #   print(output$lambda)
    
    ##debugging stuff
    if(F){
    # maybe do this in c and. you also have to recalculate the logmarginal.
    if(k>1){
      if(k>nrow(X)){
        cat(c('Matrix can not have ',k,'eigenvalues. Using ',nrow(X),' instead\n'))
        k =nrow(X)
      }
      b = matrix(runif(nrow(X),0.001,1) ,nrow=nrow(X),ncol=1) ##create enough space..
      alphas = matrix(0 ,nrow=k,ncol=1)
      betas = matrix(0 ,nrow=k-1,ncol=1)
  #    V = matrix(0,nrow=k,ncol=k);
      output =.C("Clanczos",
                 mat1 = as.double(X),
                 numRows  = as.integer(nrow(X)),
                 numCols  = as.integer(ncol(X)),
                 b = as.double(b),
                 k = as.integer(k),
                 alphas = as.double(alphas),
                 betas = as.double(betas),
                 orders = as.double(orders),
                 sigma = as.double(sigma) )
      ##craete matrix.. and call eigen...
    #  TT= diag(as.numeric(output$alphas) )
    #  print('output$alphas')
    #  print(output$alphas)
    #  for(i in 1:length(output$betas)){
    #    TT[i,i+1] = output$betas[i]
    #    TT[i+1,i] = output$betas[i]
    #  }
    #  eigenv = eigen(TT)
    #  ret$lambda=eigenv$values
 
    #  ret$vectors=matrix(output$V,nrow=k,ncol=k) %*% eigenv$vectors
      ######try sturm
      
      
      ###check c implementation real quick
      #void getEigenSturm(double *alpha,int *nalpha,double*beta,int* k,double *result){
      result =matrix(0, nrow=k,ncol=1)
      output =.C("CgetEigenSturm",
                 alpha = as.double(output$alphas),
                 nalpha  = as.integer(length(output$alphas)),
                 beta  = as.double( output$betas),
                 k = as.integer(k),
                 result = as.double(result))  
      
      print('C sturm eigen')
      print(output$result);
      
      eigenvals = output$result
      #now call inverse iterations
      # void CinverseInteration( double* A,int* nrow,int*ncol,double *lambda,double *result,double *orders,double *sigma){
      eigenvecs =matrix(0,nrow=nrow(X),ncol=length(eigenvals))
      for( i in 1:length(eigenvals)){ ##for each eigenvalue
        res = matrix(0, nrow=nrow(X),ncol=1)
        output =.C("CinverseInteration",
                   A = as.double(X),
                   nrow  = as.integer(nrow(X)),
                   ncol  = as.integer(ncol(X)),
                   lambda = as.double(eigenvals[i]),
                   result = as.double(res),
                   orders = as.double(orders),
                   sigma = as.double(sigma))
        eigenvecs[,i]=output$result;
      #  print(output$result)
      }
      
      ret$lambda = eigenvals
      ret$vectors = eigenvecs
      ####recalculate log marginal...
      N = nrow(X)
      mu1 = N*sigma + sum(X)
      
      beta= ret$lambda[1]
      mu2 =  sum( ret$lambda^2)

      t_bar = (beta*mu1 - mu2) /(beta*N - mu1)
      t_bar2 = t_bar^(2)
      beta2 = beta^2
      t1 = matrix(c(log(beta),log(t_bar)),nrow=1,ncol=2)
      t2 =(1/(beta*t_bar2 - t_bar*beta2 ))* matrix(c(t_bar2,-beta2,-t_bar,beta ),nrow=2,ncol=2 ) ## write inverse directly
      t3 = matrix(c(mu1,mu2),ncol=1,nrow =2)
      
      logdetM =  t1%*%t2%*%t3
 
      ret$approxmarginal= 0.5*t(Y)%*%ret$alpha + logdetM/2 + N/2*log(2*pi)
      
    }
    
    }
    
    return(ret)

}