##gpHistPredict

gpHistPredict =function(GP,X,x_pred){
  
  # Make sure GP is right.
  if(!is.list(GP) || is.null(GP$orders) || is.null(GP$alpha) ) {
    print("gpHistPredict(): Object appears not to be a GP!")
    return(NaN);
  }
  if ( ( ! is.matrix(X) ) || ( ! is.matrix(x_pred) ) ) {
    print("gpHistPredict(): input must be matrices!")
    return(NaN);
  }
  if(ncol(X)!= ncol(x_pred)){
    print("gpHistPredict(): GP data and prediction data does not match!")
    return(NaN);
  }
  if( nrow(GP$alpha)!= nrow(X) ) {
    print("gpHistPredict(): GP object and data X dimension missmatch!")
    return(NaN);
  }
  nsamples = nrow(x_pred)
  multResult = matrix(rep(0,nsamples),nrow=nsamples,ncol=1)
  
  output =.C("CgpHistPredict",
             result = as.double(multResult),
             mat1 = as.double(X),
             numRows  = as.integer(nrow(X)),
             numCols  = as.integer(ncol(X)),
             mat2 = as.double(x_pred),
             numRows2  = as.integer(nrow(x_pred)),
             numCols2  = as.integer(ncol(x_pred)),
             mat3 = as.double(GP$alpha),
             orders = as.double(GP$orders) 
  )
  
  ret = output$result;
  return(ret)
}