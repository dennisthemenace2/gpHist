gpHist <-
function(X,Y,sigma,orders=NULL,alphas = NULL) {
  
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

   
   logmarginal = 0; 
   lambda = 0; 
   
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
         lambda = as.double(lambda)
     )
   
    ret = list('logmarginal'=output$logmarginal,'orders'=orders,'alpha'=matrix(output$result,nrow=nrowA),'sigma'=sigma,'lambda'=output$lambda)
    
    return(ret)

}


## TEST ONLY
gpfMult <- function(X,Y,sigma,orders=NULL) {
    
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
    }
    
    
    nrowA      = nrow(X)
    # ncolB      = ncol(Y)
    multResult = matrix(rep(0,nrowA),nrow=nrowA,ncol=1)
    
    output =.C("CgpfMult",
               result = as.double(multResult),
               mat1 = as.double(X),
               numRows  = as.integer(nrow(X)),
               numCols  = as.integer(ncol(X)),
               mat2 = as.double(Y),
               numRows2  = as.integer(nrow(Y)),
               numCols2  = as.integer(ncol(Y)),
               sigma = as.double(sigma),
               orders = as.double(orders)
    )
    
  
    multResult = matrix(output$result,nrow=nrowA)	
    return(multResult)
    
  }



