
gpHistVariance= function( GP,X,x_pred){
 # Make sure GP is right.
  if(!is.list(GP) || is.null(GP$orders) || is.null(GP$alpha) ) {
    print("gpHistVariance(): Object appears not to be a GP!")
    return(NaN);
  }
  if ( ( ! is.matrix(X) ) || ( ! is.matrix(x_pred) ) ) {
    print("gpHistVariance(): input must be matrices!")
    return(NaN);
  }
  if(ncol(X)!= ncol(x_pred)){
    print("gpHistVariance(): GP data and prediction data does not match!")
    return(NaN);
  }
  if( nrow(GP$alpha)!= nrow(X) ) {
    print("gpHistVariance(): GP object and data X dimension missmatch!")
    return(NaN);
  }
  nsamples = nrow(x_pred)
  multResult = matrix(rep(0,nsamples),nrow=nsamples,ncol=1)
  
 ##check number of eigenvaules real quick
  if(length(GP$lambda)>1){ #Fine Approximation of
    output =.C("CgpHistVarianceFine",
               result = as.double(multResult),
               numRows  = as.integer(nrow(X)),
               numCols  = as.integer(ncol(X)),
               numRows2  = as.integer(nrow(x_pred)),
               numCols2  = as.integer(ncol(x_pred)),
               X = as.double(X),
               pred = as.double(x_pred),
               lambda = as.double(GP$lambda),
               nlambda = as.integer(length(GP$lambda)),
               vectors = as.double(GP$vectors),
               sigma = as.double(GP$sigma),
               orders = as.double(GP$orders) 
    )
  }else{ ## Coarse Approximation
    output =.C("CgpHistVarianceCoarse",
               result = as.double(multResult),
               mat1 = as.double(X),
               numRows  = as.integer(nrow(X)),
               numCols  = as.integer(ncol(X)),
               mat2 = as.double(x_pred),
               numRows2  = as.integer(nrow(x_pred)),
               numCols2  = as.integer(ncol(x_pred)),
               lambda = as.double(GP$lambda),
               sigma = as.double(GP$sigma),
               orders = as.double(GP$orders) 
    )
  }
  
  ret = output$result;
  return(ret)
  
}