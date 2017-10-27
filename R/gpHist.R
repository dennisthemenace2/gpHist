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
   
    return(ret)

}
