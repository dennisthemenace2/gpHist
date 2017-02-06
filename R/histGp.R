

rm(list = ls())

options(error = recover)
options(warnings = recover)

require(MASS)

fastMultiplyOne=function(X,Y){
  res = matrix(0,nrow=nrow(X),ncol=1)
  NR = nrow(X)
  
  prev = X[1:NR] * Y
 # res[NR] = sum(prev)

#  for(i in (NR-1):1){
#    int = sum(prev[1:i])
#    res[i] = int +sum(X[i]*Y[ (i+1):NR] )
#  }
  int = 0
  for(i in 1:(NR-1) ){
    res[i] = int +X[i]*sum(Y[ i:NR] )
    int = int + prev[i]
  }
  res[NR] = int+prev[NR]
  
  res
}
###only works for standard min kernel.....
fastMultiply=function(X,Y,orders, sigma,weights){
  res = matrix(0,nrow=nrow(X),ncol=1)

  for(i in 1:ncol(X)){
    ord = orders[[i]]#order(X[,i] )
    Xt= as.matrix(X[ord,i] )
    Yt= as.matrix(Y[ord])
    tres = fastMultiplyOne(Xt,Yt)#+ (sigma/ncol(X)) *Yt
    res[ord]= res[ord] +  weights[i]* tres 
  }
 # print(res)
  res = res + sigma*Y ##add sigma
#  print(res)
  
  res
}


#####conjugate gradient descnet only USE with kernel matrices! and must be orderd!
conjgradFP = function(A,b,x,orders,sigma,weights,max=1000){
  
  r=b-fastMultiply(A,x,orders,sigma,weights) #A%*%x;
  p=r;
  rsold=t(r)%*%r;
  
  for (i in 1:max ){
    Ap= fastMultiply(A,p,orders,sigma,weights) #A%*%p;
    alpha=as.numeric(rsold/(t(p) %*%Ap)) ;
    x=x+alpha*p;
    
    if(i %% 50 == 0){
    #  print(i)
      r=b-fastMultiply(A,x,orders,sigma,weights)#A%*%x;
      
    }else{
      r=r-alpha*Ap;
    }
    
    
    rsnew=t(r)%*%r;
    if(sum(sqrt(rsnew)) <1e-12){
    #  cat(c('break:', i, ' err:',sum(sqrt(rsnew)),'\n') )
      break;
    }
    p=r+as.numeric(rsnew/rsold)*p;
    rsold=rsnew;
  }
  
  x
}



downhillsimplex = function(fc,N,bp=NULL,lower=0,upper=1,it=1000,tol=1e-10){
  
  
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
  
  if(N==1){
    print('this will not work for 1 parameter only')
  }
  if(length(lower)<N){
    lower = rep(lower,N)downhillsimplex
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
    if(N>1){
      m = colMeans(bp[1:N,1:N])
    }else{
      m = sum(bp[,1])
    }
    
  #  print(bp[1,N+1]-bp[N+1,N+1] )
    if(abs(bp[1,N+1]-bp[N+1,N+1]) < tol){
      cat(c('converged after:',i,' iterations. value:', abs(bp[1,N+1]-bp[N+1,N+1]) ,'\n' ) )
      break;
    }
    
    r =(1+alpha)*m - alpha*bp[N+1,1:N]
   # print('querry r')
  #  print(r)
    r = chkBound(r,lower,upper )
    
    np = fc(r)
    if(np<bp[1,N+1]){ ##better
      e =(1+gamma)*m -gamma*bp[N+1,1:N]
   #   print('querry e')
  #    print(e)
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
      next; ##jump back
    }
    
    if(np<bp[N,N+1]){
      bp[N+1,1:N] = r
      bp[N+1,N+1] = np
      next;
    }
    
    if(np < bp[N+1,N+1]){ ##np better
      h = np
    }else{
      h =  bp[N+1,1:N]
      np = bp[N+1,N+1]
    }
    cp = beta*m+(1-beta)*h
   # print('querry cp')
  #  print(cp)
    cp = chkBound(cp,lower,upper )
    nc = fc(cp)
    if(nc< bp[N+1,N+1]){ ##c is better
      bp[N+1,1:(N+1) ] = c(cp,nc)
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
    
  }
#  cat(c('iterations:',i,'\n') )
  if(i == it){
    cat(c('potentially not converged. Diff:',abs(bp[1,N+1]-bp[N+1,N+1]) ,'\n') )
  }
  bp
}




X = as.matrix(1:5)
Y = as.matrix(X^2*-0.1 +4)
N = nrow(X)

X = as.matrix(1:5)
X = cbind(X,c(3,1,9,10,0))
Y = as.matrix(X[,1]^2*-0.1 +4 + X[,2])

N = nrow(X)


##shuffle
sh = sample(length(X))
X = as.matrix(X[sh])
Y = as.matrix(Y[sh])

x_pred = as.matrix(1.5)


x_pred = as.matrix(c(1.5,1.5) )
x_pred = t(x_pred)


minkern = function(x,y){
  res=c()
  for(i in 1:length(x)){
    res = c(res,min(x[i],y[i]))
  } 
  sum(res)
}

paramMinkern = function(x,y,w){
  res=c()
  for(i in 1:length(x)){
    res = c(res, w[i]*min( x[i] ,y[i]) )
  } 
  sum(res)
}

gpVariance= function( GP,X,x_pred){
  N= nrow(X)

  if(GP$approx){
  
    ##calc sigma
    k_dist_fast = c()
   #ksKks=c()
    for(i in 1:nrow(x_pred) ){
      ##each dim
      C =0
   # ksKks[i] =0
      for(d in 1:ncol(x_pred)){
        ord = GP$orders[[d]]
        tmp = X[ord,d]
        
        idx= which(tmp <x_pred[i,d]  ) 
        if(length(idx)>0){
          idx = max(idx)
          A = sum( (GP$weights[d]*tmp[1:idx])^2)
        #  ksKks[i] = ksKks[i] + sum( (GP$weights[d]*tmp[1:idx])^2 * 1/GP$lambdas$values[1:idx] )
     #     ksKks[i] = ksKks[i] + sum( ( matrix(GP$weights[d]*tmp[1:idx],nrow=1)%*%rt$u[1:idx,1:idx] )^2 * 1/GP$lambdas$values[1:idx] )
        }else{
          A=0
          idx=0
        #  A= A + sum( (GP$weights[d]*tmp[1])^2)
        #  C= C + sum( (GP$weights[d]*tmp[N])^2)
          
          #idx = 1
         # next;
        }
        B = (N-(idx) )*( (GP$weights[d]*x_pred[i,d] ) ^2)
        
        #ksKks[i] = ksKks[i] + sum( (GP$weights[d]*x_pred[i,d])^2 * 1/GP$lambdas$values[(idx+1):N ] )
      #  ksKks[i] = ksKks[i] + sum( (matrix(rep(GP$weights[d]*x_pred[i,d],N-idx) ,nrow=1)%*% rt$u[(idx+1):N,(idx+1):N ]) ^2 * 1/GP$lambdas$values[(idx+1):N ] )
        
      #   (1/GP$lambdas$values[1])
        
        C = C+ (A+B)
        
      }
      k_dist_fast = c(k_dist_fast ,C)
    }
    
    ###
    #K_starStar = outer(1:(nrow(x_pred)), 1:nrow(x_pred), FUN =Vectorize( function(i,j) paramMinkern(x_pred[i,],x_pred[j,],GP$weights) ) )
    #K = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) paramMinkern(X[i,],X[j,],GP$weights) ) ) ##with params
    # K = K +diag(N)*GP$sigma
    #
    #  k_star =  outer(1:(N), 1:nrow(x_pred), FUN =Vectorize( function(i,j) minkern(X[i,],x_pred[j,],GP$weights) ) )
    
    #L =chol(K)
    #sL = ginv(L)#solve(L)
    #v = t(sL) %*% k_star
    #
    
    #V =K_starStar - t(v)%*%v
    
    
    sum2 =  (1/GP$lambdas$values[1]) *(k_dist_fast )
    #sum2 = 0
    #for(i in 1:length(GP$lambdas$values)){
    #  sum2 = sum2+ (1/GP$lambdas$values[i]) *(k_dist_fast )
    #}
    
  #  K_starStarFast = rowSums(x_pred%*%GP$weights) ### this is only the diagonal
   # transXp = GP$gx(x_pred[,1], GP$px[1] )
  #  if(ncol(x_pred) >1){
  #    for(i in 2:ncol(x_pred)){
  #      transXp = cbind( transXp,GP$gx(x_pred[,i], GP$px[i] ) )
  #    }
  #  }
    
    # k_star_dist = sum(k_star^2)
    #
  
    K_starStarFast = rowSums( x_pred %*%GP$weights) ### this is only the diagonal
    
    var = K_starStarFast -sum2 + GP$sigma
    if(any(var<0) ){
      print('variance is negative fix this')
    #  stop()
      var = abs(var)
    }
  }else{
    
    K_starStar = outer(1:(nrow(x_pred)), 1:nrow(x_pred), FUN =Vectorize( function(i,j) GP$kernel(x_pred[i,],x_pred[j,],GP$weights) ) )
    k_star =  outer(1:(N), 1:nrow(x_pred), FUN =Vectorize( function(i,j) GP$kernel(X[i,],x_pred[j,],GP$weights) ) )
    
    v = t(GP$sL) %*% k_star
    var= diag(K_starStar - t(v)%*%v)
    if(any(var<0)){
      print('variance is negative fix this')
      stop()
      var = abs(var)
    }
    
  }
  
  var
}


gpHistPredict =function(GP,X,x_pred){

  NRA = nrow(GP$alpha)
  
  if(GP$approx){
  
    pred = c()
    for(i in 1:nrow(x_pred) ){
      ##each dim
      C=0
      for(d in 1:ncol(x_pred)){
        
        ord = GP$orders[[d]]#order(X[,d])
        tmp = X[ord,d]
        
        idx= which(tmp < x_pred[i,d]) 
        if(length(idx)>0){
          idx = max(idx) 
          if(idx==NRA){###largest
            C = C+ GP$weights[d]* (GP$alpha[ord][1:idx] %*% tmp[1:idx])
         #   C = C+ sum(GP$alpha[ord][1:idx]) %*% x_pred[i,d]
            next;
          }else{
            A = GP$alpha[ord][1:idx] %*% tmp[1:idx]
          }
        }else{ ###point even further left....
          A=0
          idx = 0;
          #C= C+ GP$weights[d]* GP$alpha[ord][1] ### ???
          ##can i predict most left value ??
        #  next;
        }
        B = sum(GP$alpha[ord][(idx+1):NRA]) * x_pred[i,d]
        C = C+ GP$weights[d]*(A+B)
      }
      if(is.na(C)){
        print('prediction is na')
        stop()
      }
      pred = c(pred,C)
    }
  }else{
    k_star =  outer(1:(nrow(X)), 1:nrow(x_pred), FUN =Vectorize( function(i,j) GP$kernel(X[i,],x_pred[j,],GP$weights) ) )
    pred = t(k_star)  %*% GP$alpha
    
  }
  pred
}


DEBUG = TRUE
x_pred = matrix(0)

gx=function(X,P){
  for(i in 1:nrow(X)){
    X[i,] = abs(X[i,])^P
  }
  X
}

gexp = function(x,p){
  (exp(p*abs(x))) / (exp(x))
}

gx = function(X,P){
  for(i in 1:nrow(X)){
    X[i,] = (exp(P*abs(X[i,]))) / (exp(X[i,]))
  }
  X
}

gx = function(X,P){
  for(i in 1:nrow(X)){
    X[i,] = (exp(P*abs(X[i,]))-1) / (exp(P)-1)
  }
  X
}

gx = function(X,P){
  for(i in 1:nrow(X)){
    X[i,] = (P+X[i,])
  }
  X
}

####make some hist efficent gp
gphist=function(X,Y,sigma,minkern,weights,orders=NULL,alpha_prev=NULL,approx=T){
  N = nrow(X)
  P = ncol(X)
  ###just tranfom g
  
  if(is.null(alpha_prev) ){
    alpha_prev = matrix(0,nrow=nrow(Y))
  }else{
    if(nrow(alpha_prev)!=nrow(Y)){
    #  print('alpha dimension missmatch')
      alpha_prev = rbind(alpha_prev,matrix(rep(0,nrow(Y)-nrow(alpha_prev) ) ) )
    #  alpha_prev = rbind(alpha_prev,matrix(rep(alpha_prev[length(alpha_prev)],nrow(Y)-nrow(alpha_prev) ) ) )
    }
  }
  
 # cat(c('samples:',N, ' features:',P ,'\n' ) )
  # N = length(X)
  noise = sigma *diag(N)
  
  if(!is.matrix(weights)){ ##check for matrix
    weights = matrix(weights)
  }
  
  ###order
#  ord =order(rowSums(X) )
#  X= as.matrix(X[ord,])
#  Y= as.matrix(Y[ord])
  
  ##create order list
  if(is.null(orders)){
    orders = list()
    for(i in 1:ncol(X)){
      orders[[i]] = order(X[,i])  
    }
  }
  
 ##debug
  if(DEBUG){
    N = nrow(X)
    K = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(X[i,],X[j,],weights) ) ) ##with params
    K = K +noise
    
    ###test kernel
    if( any( abs(K%*%Y- fastMultiply(X,Y,orders,sigma,weights))>1e-4 )  ){
      print('error multiplication missmatch')
     # stop()
    }
    
  }
  
  ####for testing
  if(F){
    ord =order(X[,1] )
    Xt= as.matrix(X[ord,1])
    K1 = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(Xt[i],Xt[j],weights) ) )
    
    ord =order(X[,2] )
    Xt= as.matrix(X[ord,2])
    K2 = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(Xt[i],Xt[j],weights) ) )
  }
  
  ####
  if(DEBUG){
    K_inv = ginv(K)#solve(K)
    alpha = K_inv %*% Y
  }
  
  ###calc alpha directly 
  if(approx){
    alpha_direct = conjgradFP(X,Y,alpha_prev,orders,sigma,weights)  # conjgrad(K,Y,matrix(0,nrow=nrow(Y)) )
  }else{
    N = nrow(X)
    K = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(X[i,],X[j,],weights) ) ) ##with params
    K = K +noise
    
    K_inv = ginv(K)#solve(K)
    alpha_direct = K_inv %*% Y
    
    
  }
  
  if(DEBUG){
    if(any(abs(alpha_direct-alpha)>1e-3 )){
      print('error alpha missmatch')
      stop()
    }
  
    #minkern=function(x,y,w){
    #  res=c()
    #  for(i in 1:length(x)){
    #    res = c(res,w[i]* max(x[i],y[i]))
    #  } 
    #  sum(res)
    #}
    if(ncol(x_pred)<ncol(X) ){
      x_pred = t(matrix(rep(0,ncol(X))) )
    }
    
    k_star =  outer(1:(N), 1:nrow(x_pred), FUN =Vectorize( function(i,j) minkern(X[i,],x_pred[j,],weights) ) )
    ##
    
    ##
    
    predictions = t(k_star)  %*% alpha
    
    GPt = list('alpha'=alpha_direct,'orders'=orders,'weights'=weights,'sigma'=sigma,'approx'=T)
    pred_test =  gpHistPredict(GPt,X,x_pred)
    if(abs(predictions - pred_test)>0.0001){
      print('prediction missmatch')
      stop()
    }
   
  }
  ### find index

  if(DEBUG){
  #  GP = list('orders'=orders,'alpha'=alpha,'weights'=weights,'sigma'= sigma)
  #  pred = gpHistPredict(GP = GP,X= X , x_pred = x_pred)
    
#
  #  if(any(abs(pred - predictions)>1e-9 )) {
  #    print('error prediction differ')
  #    stop()
  #  }
  
    L =chol(K)
    sL = ginv(L)#solve(L)
    v = t(sL) %*% k_star
    K_starStar = outer(1:(nrow(x_pred)), 1:nrow(x_pred), FUN =Vectorize( function(i,j) minkern(x_pred[i,],x_pred[j,],weights) ) )
    
    V =outer(1:(nrow(x_pred)), 1:nrow(x_pred), FUN =Vectorize( function(i,j) minkern(x_pred[i,],x_pred[j,],weights) ) ) - t(v)%*%v
   
  }
  
  ####
  # for each dimension 
  ####
  
  #X_trans = gx(X[,1],px[1])
  #if(ncol(X)>1){
  #  for(i in 1:ncol(X)){
  #    X_trans = gx(X[,i],px[i])
  #  }
  #}
  if(approx){
    mu1 = N*sigma + sum(X %*%weights)  #### 
    
    if(DEBUG){
      if(abs(sum(diag(K)) -mu1)>1e-6){
        print('error calculating mu1')
        stop()
      }
    }
    ####get eigen values !
    
    b = matrix(runif( N,0.5,1)) ##mabe do some better stuff here to get eigenvalues
    TT= lanczosKERN(X,b,orders,sigma,weights)
    
    lambdas = eigen(TT, only.values = TRUE)
    beta= lambdas$values[1]
    mu2 =  sum(lambdas$values[1]^2)
   
     
    if(DEBUG){
      ###check real quick
      if(abs(prod(lambdas$values)-det(K))>1e-10 ){
        #print('error eigen values not matching')
        if(abs(lambdas$values[1]-eigen(K)$values[1])>1e-4 ){
          print('first eigenvalue not matching')
          #stop()      
        }
        
      }
    }
    t_bar = (beta*mu1 - mu2) /(beta*N - mu1)
    if(is.na(t_bar)){
      print('t_bar is na !')
      stop()
    }
    t_bar2 = t_bar^(2)
    
    beta2 = beta^2
    
    t1 = matrix(c(log(beta),log(t_bar)),nrow=1,ncol=2)
    #t2 = matrix(c(beta,beta^2,t_bar,t_bar^(2) ),nrow=2,ncol=2 )
    
    t2 =(1/(beta*t_bar2 - t_bar*beta2 ))* matrix(c(t_bar2,-beta2,-t_bar,beta ),nrow=2,ncol=2 ) ## write inverse directly
    
    t3 = matrix(c(mu1,mu2),ncol=1,nrow =2)
      
  #  logdetM =  t1%*%solve(t2) %*%t3
    logdetM =  t1%*%t2%*%t3
    print(logdetM)
  }else{
    L = chol(K)
    sL = ginv(L)#solve(L)
    
    logdetM =  sum(log(diag(L)))*2
    ###for compatibility reasons
    lambdas  = eigen(K)
  }
  
  
  ##calc sigma
#  sum1 = 0
#  vecsum= 0
 #   if(length(lambdas$values)>1){
#      for(i in length(lambdas$values):2){
#        sum1 =sum1 +(1/lambdas$values[i]) * ( X%*%weights)
#        vecsum = vecsum +  X%*%weights
#      }
#    }
  
 # k_star2 = k_star^2# %*% K
#  k_dist = sum(k_star2) 
  ## cal k_dist
  
 # k_dist_fast = c()
#  for(i in 1:nrow(x_pred) ){
    ##each dim
#    C =0
#    for(d in 1:ncol(x_pred)){
#      ord = orders[[d]]
#      tmp = X[ord,d]
#      idx= which(tmp <x_pred[i,d] ) 
#      if(length(idx)>0){
#        idx = max(idx)
#        A = sum( (weights[d]*tmp[1:idx])^2)
#      }else{
#        A=0
#        idx=0
#      }
#      B = (N-(idx) )*( (weights[d]*x_pred[i,d]) ^2)
#      C = C+ (A+B)
#    }
#    k_dist_fast = c(k_dist_fast ,C)
#  }
  
  #if(k_dist_fast> k_dist){
  #  print('error')
  #  stop()
  #}
  
  ###
 # sum2 =  (1/lambdas$values[1]) *(k_dist_fast)
#  K_starStarFast = rowSums(x_pred%*%weights) ### this is only the diagonal
  
#  if(DEBUG){
#    if(any(K_starStarFast!=diag(K_starStar))){
#      print('kstar missmatch')
#      stop()
#    }
#  }
  ###
 # var = K_starStarFast - sum1 -sum2 + sigma
 
 # var2 = gpVariance(X,orders,weights,lambdas, sigma,x_pred)
#  if(any(var!=var2) ){
#    print('var should be equal')
#    stop()
#  }
 
  
  ####real var
  #realVar = K_starStar - t(k_star) %*%K_inv %*%k_star
  
  #####
  
   
 # logmarginal  = -0.5*t(Y) %*%alpha - sum(log(diag(L))) #- nrow(X)/2*log(2*pi)
 # logmarginal  = 0.5*t(Y) %*%alpha + sum(log(diag(L))) #- nrow(X)/2*log(2*pi)
  
#  logmarginal2  = t(Y) %*%alpha + sum(log(diag(L)))
  
  approxmarginal= 0.5*t(Y)%*%alpha_direct + logdetM/2 + nrow(X)/2*log(2*pi)
  #print(' nrow(X)/2*log(2*pi)')
  #print( nrow(X)/2*log(2*pi))
  print('t(Y)%*%alpha_direct')
  print(t(Y)%*%alpha_direct)
  
  print("approxmarginal")
print(approxmarginal)
#  cat(c('marginal approx dif:',approxmarginal - logmarginal , '\n'))
  ##get log det bound
  
  ##optimze ?
  if(!approx || DEBUG){
    ret = list('logmarginal'=approxmarginal,'orders'=orders,'alpha'=alpha_direct,'weights'=weights,'lambdas'=lambdas,'sigma'=sigma,'approx'=approx,'L'=L,'sL'=sL,'kernel'=minkern)
  }else{
    ret = list('logmarginal'=approxmarginal,'orders'=orders,'alpha'=alpha_direct,'weights'=weights,'lambdas'=lambdas,'sigma'=sigma,'approx'=approx)
  }
  ret
}

gpHistUpdata = function(X,Y,X_old,Y_old,gpList){
  ###find entries for new values
  for(d in 1:ncol(X)){
    for(i in 1:nrow(X)){gphist(X ,Y,p[1],paramMinkern,c(p[2]));
       idx = which(X[i,d]> X_old[gpList$orders[[d]],d] )  
       if(length(idx)>0){ ##found
         idx = max(idx) +1
         if(idx>nrow(X)){
           gpList$orders[[d]] = c(gpList$orders[[d]],idx)  
         }else{
           incx = which(gpList$orders[[d]]>=idx ) 
           gpList$orders[[incx]] =  gpList$orders[[incx]] +1
           gpList$orders[[d]] = c(gpList$orders[[d]],idx)  
         }
       }else{ ##is first
         gpList$orders[[d]] = gpList$orders[[d]] +1
         gpList$orders[[d]] = c(gpList$orders[[d]],1)
       }
    }
  }
  
}


####ONLY FOR KERNEL
lanczosKERN= function(A,b,orders,sigma,weights){
  
  v_prev = matrix(0,nrow=nrow(A), ncol=1)
  v= b /as.numeric(sqrt(t(b)%*%b) )
  
  beta = matrix(0,nrow=1, ncol=nrow(A))
  alpha = matrix(0,nrow=1, ncol=nrow(A))
  
  for(i in 1:(nrow(A)-1 ) ){
    w = fastMultiply(A,v,orders,sigma,weights)#A%*%v
    alpha[,i] =t(w)%*%v
    w = w-alpha[,i]*v - beta[,i]*v_prev
    v_prev = v
    beta[,i+1]= sqrt(t(w)%*%w)
    v=w/beta[,i+1]
  }
  
  w= fastMultiply(A,v,orders,sigma,weights)#A%*%v
  alpha[,nrow(A)] =t(w)%*%v
  
  TT= diag(as.numeric(alpha) )
  for(i in 2:ncol(beta)){
    TT[i,i-1] = beta[,i]
    TT[i-1,i] = beta[,i]
  }
  
  TT
}

###test stufff

mapfc = function(p){
  if(any(p<0) ){
    print('negative.. all need to be positive')
    print(p)
    return(99999)
    #stop()
  }  
  X= matrix(c(1:10))
  set.seed(12345)
  Y = matrix(c(X^2*0.4+3+rnorm(10,0,0.01) ) )
  #gpist(X,Y,sigma,minkern,weights, x_pred)
  
  ret = gphist(X ,Y,p[1],paramMinkern,c(p[2]) ,matrix(1) );
  ret$logmarginal
  points(1,ret$predictions[1],col='red')
  ret$logmarginal
}

resmat = downhillsimplex(mapfc,2,it = 1000,tol=0.01)

#res =  optim(c(2,2), mapfc, method = c("Nelder-Mead") )



mapfc2 = function(p){
  if(any(p<0) ){
    print('negative.. all need to be positive')
    print(p)
    return(99999)
    #stop()
  }  
  X= matrix(c(1:10))
  X = cbind(X,c(10:1))
  set.seed(12345)
  Y = matrix(c(X[,2]+X[,1]^2*0.4+3+rnorm(nrow(X),0,0.1) ) )
  #gpist(X,Y,sigma,minkern,weights, x_pred)
  x_pred = matrix(c(1:5))
  x_pred = cbind(x_pred, matrix(c(1:5) ) )
  
  ret = gphist(X ,Y,p[1],paramMinkern,c(p[2],p[3]) ,x_pred );
  #ret$logmarginal
  #points(1,ret$predictions[1],col='red')
  ret$logmarginal
}

resmat = downhillsimplex(mapfc2,3,it = 100,tol=0.01)

res =  optim(c(2,2,1), mapfc2, method = c("Nelder-Mead") )

X1 = matrix(X[1:10,],ncol=2)
Y1 = matrix(Y[1:10])

sigma = resmat[1,1]
weights = resmat[1,2:3]

gpres = gphist(X,Y,sigma,paramMinkern,weights, x_pred,alpha_prev=NULL);
  
X=rbind(X,c(12,13))
Y = rbind(Y,c(13+12^2*0.4+3+rnorm(1,0,0.1) ) )

gpres_full = gphist(X,Y,sigma,paramMinkern,weights, x_pred,alpha_prev=NULL);

gpres_part = gphist(X,Y,sigma,paramMinkern,weights, x_pred,alpha_prev=gpres$alpha);

###
y_max = max(Y)
kappa = 2.576
eps = 0.0
acq = "ucb"
ndim = 2
npoints = 10
####

testBO =function(){
  ##create data
  #X= matrix(c(10:20))
  X= matrix(c(0.001,0.2))
  X= matrix(c(seq(0.0,1,0.1)))
  Y = matrix(dnorm(X,0.5,0.1))
 # X = cbind(X,c(10:1))
  set.seed(12345)
  getY = function(X){
    Y= matrix((c(X-0.5)^2 ) )
    Y
  }
  getY = function(X){
    Y= matrix( sin( 2*pi*X*2 ) )
    Y
  }
  
  
#  Y = matrix(c(X[,2]+X[,1]^2*0.4+3+rnorm(nrow(X),0,0.1) ) )
  Y = getY(X)
#  x_pred = matrix(c(1:5))
#  x_pred = cbind(x_pred, matrix(c(1:5) ) )

  GP = gphist(X ,Y,1,paramMinkern,c(1));

    
  GP = NULL
  dhfc = function(p){
    if(any(p<0) ){
      print('negative.. all need to be positive')
      print(p)
      return(99999)
      #stop()
    }  
    X_trans = gx(X,p[3])
   # if(is.null(GP)){
  #    ret = gphist(X_trans ,Y,p[1],paramMinkern,c(p[2]));
  #    GP<<-ret
  #  }else{
   #   ret = gphist(X_trans ,Y,p[1],paramMinkern,c(p[2]), orders = GP$orders, alpha_prev = GP$alpha); 
   #   GP<<-ret
    #}
    ret = gphist(X_trans ,Y,p[1],paramMinkern,c(p[2]));
    GP<<-ret
    
    ret$logmarginal
  }
  resmat = downhillsimplex(dhfc,3,it = 100,tol=0.001)
  ##found parameters
  
  ###FOR DEBUG ONLY
  weights = GP$weights
  orders = GP$orders
  sigma = GP$sigma
  alpha_direct = GP$alpha
  alpha =  GP$alpha
  px = GP$px
  #####DEBUG ONLY
  
  ##plot some stuff
  
  for( i in 1:10){
    #plot(X,Y,xlim=c(0,20),ylim=c(-10,200))
    plot(X,Y)
    
    title(i)
    predPoints = matrix(seq(0,1,0.1))
    X_trans = gx(X,resmat[1,3])
    x_pred_trans = gx(predPoints , resmat[1,3])
    
    preds = gpHistPredict(GP=GP,X=X_trans,x_pred = x_pred_trans)
    #points(predPoints,preds,col='red')
    lines(predPoints,preds,col='red')
    vars =  gpVariance(GP=GP,X=X_trans,x_pred = x_pred_trans)
    
    lines(predPoints,preds-sqrt(vars),col='green' )
    lines(predPoints,preds+sqrt(vars),col='green' )
    
    
   # utmat = Utility_Max(10,X,GP,max(Y))
    utmat = Utility_Max(10,X_trans,GP,max(Y),gx = gx,px=resmat[1,3])
    
    ##which is the min...
    idx= which.min(utmat[,2])
    if(length(idx)<=0){
      print('error')
    }
    ##querry points
    newX = matrix(utmat[idx,1])
    newY = getY(newX )
    X= rbind(X, newX)
    Y= rbind(Y,newY )
    ##update gaussian process
    points(newX,newY,col='blue')
    
    px = resmat[1,3]
    X_trans = rbind(X_trans, gx(newX,resmat[1,3]) ) 
    
    GP = gphist(X_trans ,Y,GP$sigma,paramMinkern,GP$weights,alpha_prev = GP$alpha ); ##think about using old orders....
  ##retrain ?
  }
  
  #GP =  gphist(X,Y,sigma,paramMinkern,weights,alpha_prev=NULL);
}


Utility <- function(x_vec, X,GP,y_max, acq = "ucb", kappa, eps,gx=NULL,px=NULL) {
  # Gaussian Process Prediction
  #GP_Pred <- gpHistPredict(GP$alpha, GP$orders, X,GP$weights ,x_vec)#GPfit::predict.GP(object = GP, xnew = matrix(x_vec, nrow = 1))
#  print('hello')
  
  if(!is.null(gx)){
  #  print('teansforme')
    testpoint = matrix(gx(matrix(x_vec,nrow=1),px),nrow=1)  
  }else{
    testpoint = matrix(x_vec,nrow=1)  
  }

  ##for each row
#  print(X)
  GP_Mean <-gpHistPredict(GP ,X,testpoint ) # GP_Pred$Y_hat
  GP_MSE <- gpVariance( GP,X,testpoint)#GP_Pred$MSE %>% pmax(., 1e-9)
 # print(GP_Mean)
#  print(GP_MSE)
  if(GP_MSE<0){
    print('error GP_MSE < 0')
    print(x_vec)
    stop()
  }
  # Utility Function Type
  if (acq == "ucb") {
   # print('acq')
    value <- GP_Mean - kappa * sqrt(GP_MSE)
    #  GP_Mean + kappa * sqrt(GP_MSE)
    
  } else if (acq == "ei") {
    z <- (GP_Mean - y_max - eps) / sqrt(GP_MSE)
    value <- (GP_Mean - y_max - eps) * pnorm(z) + sqrt(GP_MSE) * dnorm(z)
  } else if (acq == "poi") {
    z <- (GP_Mean - y_max - eps) / sqrt(GP_MSE)
    value <- pnorm(z)
  }
  return(value)
}
###try to find next value...
Utility_Max <- function(npoints, X,GP,  y_max,lower=0,upper=1,acq=  "ucb", gx=NULL,px=NULL,kappa= 2.576, eps=0.0) {
  # Try Different Initial Values
  ndim = ncol(X)
  if(length(lower)==1){
    testmat <- matrix( runif(ndim*npoints , lower ,upper) ,nrow=npoints, ncol=ndim )
    if(ndim>1){
      lower = rep(lower,ndim )
      upper = rep(upper,ndim )
    }
  }else{
    testmat = matrix(runif(ndim , lower[1] ,upper[1]) ,nrow=npoints,ncol=1 )
    for(i in 2:length(lower)){
      testmat = cbind(testmat,runif(ndim , lower[i] ,upper[i])  )
    }
  }
    #Matrix_runif(100, #  lower = rep(0, length(DT_bounds[, Lower])),
                          #  upper = rep(1, length(DT_bounds[, Upper])))
  
  # Negative Utility Function Minimization
  #Mat_optim <- foreach(i = 1:nrow(testmat), .combine = "rbind") %do% {
  ##prepare set matrix
  remat = matrix(0,nrow=npoints,ncol=ndim+1)
  for(i in 1:nrow(testmat)){
    optim_result <- optim(par = testmat[i,] ,
                          fn = Utility,
                          X= X, GP = GP, acq = acq, y_max = y_max,gx=gx,px=px, kappa = kappa, eps = eps,
                          method = "L-BFGS-B",
                          lower = lower,
                          upper = upper,
                          control = list(maxit = 100,factr = 5e11))
    remat[i,]=  c(optim_result$par, optim_result$value)
  }

  #} %>%
  #  data.table(.) %>%
  #  setnames(., old = names(.), new = c(DT_bounds[, Parameter], "Negetive_Utility"))
  # Return Best set of Parameters
  #argmax <- as.numeric(Mat_optim[which.min(Negetive_Utility), DT_bounds[, Parameter], with = FALSE])
#  return(argmax)
  remat
}
############
fc = function(p){
  ret = 100*(p[2]-p[1]^2)^2+(1-p[1])^2
}

fc2 = function(p1,p2){
  ret = 100*(p2-p1^2)^2+(1-p1)^2
  ret = ret*-1
  list(Score = ret,Pred=ret*-1)
}


BayesianOptimization(fc2,bounds = list(p1 = c(0, 10),p2=c(0,10)), init_points = 20, n_iter = 20,
acq = "ucb", kappa = 2.576, eps = 0.0,verbose = TRUE)


fctest = function(p){
 # ret = (p^2)+2
  ret = sin(2*pi*p*0.1)-2 #+p +sin(2*pi*p*0.7) 
}

res =optimizeMe(fctest,npoints=20,ndim=1,lower=0,upper=10,0.0001,2,10,5,approx=F,kern =expkern,acq = 'ei')

gx = function(X,P){
  for(i in 1:nrow(X)){
    X[i,] = (110+X[i,])
  }
#  print(P)
  X
}



res =optimizeMe(fctest,npoints=20,ndim=1,lower=-10,upper=10,0.0001,0.1,20,5,approx=T,kern =paramMinkern,gx = gx,noGxParams = T)

res=optim(c(10,5),fc)
resmat = downhillsimplex(fc,2)
res =optimizeMe(fc,npoints=20,ndim=2,lower=0,upper=10,0.0001,2,5,1,approx=T,kern =paramMinkern)

expkern=function(xi,xj,l){
  if(length(l)>1){
    V =diag(as.numeric(l^2))
  }else{
    V=matrix(l^2)
  }
  
  x = (xi - xj)
  exp(-0.5 *  t(x)%*%V%*%(x) )
 # x = (xi - xj)^2
#  exp(-0.5 *  ( x/l^2) )
}

res =optimizeMe(fc,npoints=20,ndim=2,lower=0,upper=10,0.0001,2,5,1,approx=F,kern =expkern)


X= res$X
Y= res$Y
GP=res$GP
resmat = res$resmat
px = resmat[1,4:5]
X_trans = gx(X,px)
upper = 10
lower = 0
##### check some stuff

getParams = function(p){
  gp = gphist(X ,Y,p[1],expkern,p[2],approx=F ); ##think about using old orders....
  gp$logmarginal
}

res = optim(c(1,1),getParams,method = "L-BFGS-B",lower = 0.000001,upper = 100)

GP$approx = F  
p1 = gpHistPredict(GP,X_trans,matrix(c(40185910,3730411169),nrow=1))
GP$approx = T
p2 = gpHistPredict(GP,X_trans,matrix(c(40185910,3730411169),nrow=1))

GP$approx = F 
p3 = gpHistPredict(GP,X_trans,matrix(c(100,100),nrow=1))

p3 = gpHistPredict(GP,X_trans,matrix(c(102.7377,108.2059164),nrow=1))

GP$approx = F  
var1 = gpVariance(GP,X_trans,matrix(c(100,100),nrow=1))
GP$approx = T
var2 = gpVariance(GP,X_trans,matrix(c(100,100),nrow=1))
var3 = gpVariance(GP,X_trans,matrix(c(102.7377,108.2059164),nrow=1))


PLOT=TRUE
library(rgl)
require(akima)
require(MASS)
######optimise some crap
res =optimizeMe(getY,2,1,0,1,0.0001,10,10,1)

optimizeMe =function(getY,npoints,ndim,lower,upper,paramLower,paramUpper,it,update,kern,approx=T,gx=NULL,noGxParams=F,acq = 'ucb',kappa=2.576,eps=0.0){

  modelEvaluations = 0

  if(length(lower)==1){
    X <- matrix( runif(ndim*npoints , lower ,upper) ,nrow=npoints, ncol=ndim )
#    lower = rep(lower,ndim )
#    upper = rep(upper,ndim )
  }else{
    if(ndim != length(lower)){
      print('you have to specify a boundery for each dimension')
      return(0)
    }
    X = matrix(runif(npoints , lower[1] ,upper[1]) ,nrow=npoints,ncol=1 )
    for(i in 2:length(lower)){
      X = cbind(X,runif(npoints , lower[i] ,upper[i])  )
    }
  }
  ####remove double
  X = unique(X)

  Y = matrix(,ncol=1,nrow=nrow(X)) ##create Y
  for(i in 1:nrow(X)){
    Y[i] = getY(X[i,])
  }
  modelEvaluations = modelEvaluations +nrow(X)
  ###check parambounds
  
  if(is.null(gx) | noGxParams ){
    trueDim =ndim+1 ##for sigma
  }else{
    trueDim =ndim*2+1 ##for sigma and gx
  }
  
  if( (length(paramLower)<trueDim &length(paramLower)>1) |  (length(paramUpper)<trueDim &length(paramUpper)>1)){
    cat(c('You have to specify:',trueDim,'bounderies for the parameters! ',length(paramLower),' are defined.\n'))
    return(0)
  }
  
  if(length(paramLower)<trueDim){
    print('paramLower is changed!')
    paramLower = rep(paramLower,trueDim)
  }
  if(length(paramUpper)<trueDim){
    print('paramUpper is changed!')
    paramUpper = rep(paramUpper,trueDim)
  }
  
  ###
  
  GP = NULL
  dhfc = function(p){
    if(any(p<paramLower | p>paramUpper) ){
      print('left dimensions')
      print(p)
    #  print(paramLower)
    #  print(paramUpper)
    #  stop()
    #  return(99999)
      #stop()
    }  
    P = length(p)
    if(is.null(gx) ){
      ret = gphist(X ,Y,p[1],kern,matrix(p[ 2:(1+ncol(X))]),approx=approx);
    }else{
      if(noGxParams){
        X_trans = gx(X,0)
        ret = gphist(X_trans ,Y,p[1],kern,matrix(p[ 2:(1+ncol(X))]),approx=approx);
      }else{
        X_trans = gx(X,p[(P-ncol(X)+1):P])
        ret = gphist(X_trans ,Y,p[1],kern,matrix(p[ 2:(1+ncol(X))]),approx=approx);
      }
    }
  #  if(is.null(GP)){
      #  ret = gphist(X_trans ,Y,p[1],kern,matrix(p[ 2:(1+ncol(X))]),approx=approx);
        GP<<-ret
  #  }else{
    # ret = gphist(X_trans ,Y,p[1:ncol(X)],paramMinkern,matrix(p[(ncol(X)+1):((ncol(X)*2))]), orders = GP$orders, alpha_prev = GP$alpha);
     #ret = gphist(X_trans ,Y,p[1],paramMinkern,matrix(p[ 2:(1+ncol(X))] ),  alpha_prev = GP$alpha); 
    # GP<<-ret
    #}
    
    ret$logmarginal
  }
  if(is.null(gx)){
    resmat = downhillsimplex(dhfc,ndim+1,lower=paramLower,upper=paramUpper,it = 100,tol=0.001)
   
  }else{
    px=0;
    if(noGxParams){
      resmat = downhillsimplex(dhfc,ndim+1,lower=paramLower,upper=paramUpper,it = 100,tol=0.001)
    }else{
      resmat = downhillsimplex(dhfc,ndim*2+1,lower=paramLower,upper=paramUpper,it = 100,tol=0.001)
      px = resmat[1, (ncol(resmat)-ncol(X)):(ncol(resmat)-1) ]
      ##plot some stuff
    }
    X_trans = gx(X,px) ##maybe i can save this step
  }
  ##found parameters
 
  for( i in 1:it){
    #plot(X,Y,xlim=c(0,20),ylim=c(-10,200))
    if(PLOT & ndim==1){
      plot(X,Y,col='green',xlim=c(lower,upper))
      
      title(i)
      predPoints = matrix(seq(lower,upper,0.1))
  
      if(is.null(gx)){
        preds = gpHistPredict(GP=GP,X=X,x_pred = predPoints)
        vars =  gpVariance(GP=GP,X=X,x_pred = predPoints)
      }else{
        x_pred_trans = gx(predPoints , resmat[1,3])
        preds = gpHistPredict(GP=GP,X=X_trans,x_pred = x_pred_trans)
        vars =  gpVariance(GP=GP,X=X_trans,x_pred = x_pred_trans)
      }
    
      tp = getY(predPoints)
      lines(predPoints,tp,col='blue')
      ##plot stuff  
    
      lines(predPoints,preds,col='red')
      lines(predPoints,preds-sqrt(vars),col='black' )
      lines(predPoints,preds+sqrt(vars),col='black' )
      #####
    }else if(PLOT & ndim==2){
    #plot3d(X[,1],X[,2], Y, col="red", size=5) 
      grid = 0.5
      pd = (upper-lower)/ grid +1
      points = matrix(,nrow=pd^2, ncol=2)
      rowidx = 1
      for(xi in seq(lower,upper, grid)){
        for(yi in seq(lower,upper, grid)){
          points[rowidx,] = c(xi,yi)
          rowidx = rowidx + 1
        }
      }
      if(is.null(gx)){
        preds = gpHistPredict(GP=GP,X=X,x_pred = points)
        vars =  gpVariance(GP=GP,X=X,x_pred = points)
      }else{
        x_pred_trans = gx(points ,px) 
        preds = gpHistPredict(GP=GP,X=X_trans,x_pred = x_pred_trans)
        vars =  gpVariance(GP=GP,X=X_trans,x_pred = x_pred_trans)
      }
      
      ##plot stuff  
       
      tp=c()
      for(i in 1:nrow(points)){
        tp[i] = getY(points[i,])
      }
      plot3d(points[,1],points[,2], preds, col="red", size=5) 
      plot3d(points[,1],points[,2], tp, col="blue", size=5,add=T) 

      s=interp(points[,1], points[,2],preds)
      k=interp(points[,1], points[,2],tp)
      
      #surface3d(s$x,s$y,s$z)
      #plot3d(predict[,1], predict[,2],res[,1], col="red", size=5)
      
      #plot3d(allParams[,1], allParams[,2],allError[,1], col="red", size=5)
      persp3d(s$x,s$y,s$z, aspect=c(10, 1, 0.5),  col = "red", add=TRUE )
      persp3d(k$x,k$y,k$z, aspect=c(10, 1, 0.5),  col = "lightblue", add=TRUE )
      plot3d(X[,1],X[,2], Y, col="green", size=10,add=T) 
      
    
      v=interp(points[,1], points[,2],preds-sqrt(vars) )
      persp3d(v$x,v$y,v$z, aspect=c(10, 1, 0.5),  col = "yellow", add=TRUE )
      ###
      #plot3d(0,0,0,  col = "purple",size=20,add=TRUE )
      
    }
    
    # utmat = Utility_Max(10,X,GP,max(Y))
    if(is.null(gx)){
      utmat = Utility_Max(10,X,GP,max(Y),lower,upper,gx = NULL,px=NULL,acq = acq,kappa=kappa,eps=eps)
    }else{
      utmat = Utility_Max(10,X_trans,GP,max(Y),lower,upper,gx = gx,px=px,acq = acq, kappa=kappa,eps=eps) 
    }
    
    
    ##which is the min...
   for( w in 1:nrow(utmat)){
      idx= which.min(utmat[,ndim+1])
      if(length(idx)<=0){
        print('error')
        stop()
      }
      ##querry points
      newX = utmat[idx,1:ndim]
      if(any(X == newX) ){
      #  print('already in set, dont use this one !')
        utmat[idx,] =Inf
        idx = -1
        next
      }else{
        break;
      }
      
    }
    
    if(idx<0){
      print('no new point') #we can update now or leave
     # break;
      
      GP = NULL
      ##add some noise to resmat
      
      if(is.null(gx)){
        resmat = downhillsimplex(dhfc,ndim+1,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) 
      }else{
        if(noGxParams){
          resmat = downhillsimplex(dhfc,ndim+1,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) 
        }else{
          resmat = downhillsimplex(dhfc,ndim*2+1,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) 
          px = resmat[1, (ncol(resmat)-ncol(X)):(ncol(resmat)-1) ]
        }
        X_trans = gx(X,px)
      }
      next;
      
    #  stop()
    }
    newX = matrix(newX,nrow=1)
    newY = getY(newX ) ###get real Y
    modelEvaluations = modelEvaluations +1
    
    cat(c('new point:',newX,'\n'))
    X= rbind(X, newX)
    Y= rbind(Y,newY )
    ##update gaussian process
    if(PLOT& ndim==1){
      points(newX,newY,col='purple')
    }else if(PLOT &ndim==2){
      plot3d(newX[1],newX[2],newY, col = "purple",size=10,add=TRUE )
    }
    if(is.null(gx)){
      ###do nottin
    }else{
      X_trans = rbind(X_trans, gx(newX,px ) ) 
    }
    
    if(i %% update){
      if(is.null(gx)){
        GP = gphist(X ,Y,GP$sigma,kern,GP$weights,alpha_prev = GP$alpha , approx=approx ); ##think about using old orders....
      }else{
        GP = gphist(X_trans ,Y,GP$sigma,kern,GP$weights,alpha_prev = GP$alpha,approx=approx ); ##think about using old orders....
      }
    }else{
      GP = NULL
      ##add some noise to resmat
  #    resmat = resmat + rnorm(length(resmat),0,0.1 )
      
      R = ncol(resmat)-1
      resmat = resmat[,1:R]
      best = resmat[1,]
      for(z in 1:ncol(resmat)){
        if(length(paramLower)>1){
          pl = paramLower[z]
        }else{
          pl = paramLower
        }
        if(length(paramUpper)>1){
          pu = paramUpper[z]
        }else{
          pu = paramUpper
        }
        resmat[2:(R+1),z] = runif(R,pl,pu)
        #resmat[2:(ndim*2+2),z] = runif(ndim*2+1,paramLower[z],paramUpper[z])
       # resmat[2:(ndim*2+2),z] = runif(ndim*2+1,paramLower,paramUpper)
      }
     # resmat = downhillsimplex(dhfc,ndim*2+1,bp=resmat,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) ##update params
      
      
     # if(is.null(GP)){
    #    resmat = downhillsimplex(dhfc,ndim*3,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) ##update params
    #  }
      if(is.null(gx)){
        resmat = downhillsimplex(dhfc,ndim+1,bp=resmat,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) ##update params
        print(sum(abs(best - resmat[1,1:(ndim+1)])))
      }else{
        if(noGxParams){
          resmat = downhillsimplex(dhfc,ndim+1,bp=resmat,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) ##update params
        }else{
          resmat = downhillsimplex(dhfc,ndim*2+1,bp=resmat,lower=paramLower,upper=paramUpper,it = 100,tol=0.001) ##update params
          print(sum(abs(best - resmat[1,1:(ndim*2+1)])))
          px = resmat[1, (ncol(resmat)-ncol(X)):(ncol(resmat)-1) ]
        }
        X_trans = gx(X,px)  
      }
    }
    ##retrain ?
  }
  yminpos = which.min(Y)
  xminpos = X[yminpos,]
  
  cat(c('minimum:  x:',xminpos,' y:',Y[yminpos],'\n'  ) )
  
  #list('X'=xminpos,'Y'=Y[yminpos])
  list('X'=X,'Y'=Y,'Xpos'=xminpos,'Ymin'=Y[yminpos],'GP'=GP,'resmat'=resmat,'modelEvaluations'=modelEvaluations)
}

#
gs = function(V){
  n = nrow(V);
  k = ncol(V);
  U = matrix(0,nrow=n,ncol=k);
  
  U[,1] = V[,1]/sqrt( t(V[,1])%*%V[,1]);
  for(i in 2:k ){
     U[,i] = V[,i];
     for(j in 1:(i-1) ){
       U[,i] = U[,i] - ( t(U[,i])%*%U[,j] ) %*% U[,j];
     }
    U[,i] = U[,i]/sqrt(t(U[,i])%*% U[,i]);
  }
  U
}
##qr decompostion gramm schmidt
qrgs= function(A){
  t(gs(A)) %*%A
}
####


lanczos= function(A,b){
  
  v_prev = matrix(0,nrow=nrow(A), ncol=1)
  v= b /as.numeric(sqrt(t(b)%*%b) )
  
  beta = matrix(0,nrow=1, ncol=nrow(A))
  alpha = matrix(0,nrow=1, ncol=nrow(A))
  
  for(i in 1:(nrow(A)-1 ) ){
    w = A%*%v
    alpha[,i] =t(w)%*%v
    w = w-alpha[,i]*v - beta[,i]*v_prev
    v_prev = v
    beta[,i+1]= sqrt(t(w)%*%w)
    v=w/beta[,i+1]
  }
  
  w= A%*%v
  alpha[,nrow(A)] =t(w)%*%v
  
 # TT= diag(as.numeric(alpha) )
#  for(i in 2:ncol(beta)){
#    TT[i,i-1] = beta[,i]
#    TT[i-1,i] = beta[,i]
#  }
  
#  TT
  list('alpha'=alpha,'beta'=beta)
}


### compute tridiac eigenvalues
eigenL=function(alpha,beta){
  N = length(alpha)
  val = rep(0,N)
  for(i in 1:N){
    val[i]  = alpha[i]+ 2*sqrt(2*beta[i])*cos( (i*pi)/(N+1)  ) 
  }
  val 
}

## 
getTriEigen=function(alpha,beta){
#  if(length(alpha)>length(beta)){
#    beta = c(beta,beta[length(beta)] )
#  }
  sturm = function(alpha,beta,lambda){
    #n = nrow(M)
    #r = ncol(M)
    n = length(alpha)
   
    val = rep(1,n+1)
    
    val[2] = lambda- alpha[1] #alpha[1,1]
    for(i in 2:n){
      val[i+1] = (lambda- alpha[i])*val[i] - beta[i-1]^2*val[i-1]
    }
   
    val[2:(n+1)] 
  }

##estimate bound

  bmin = Inf
  bmax = -Inf

  for(i in 1:length(alpha)){
    #g = sum(abs(M[i,]))- M[i,i]
    if(i == 1){
      g = abs(beta[i])
    }else if(i == length(alpha)){
      g = abs(beta[i-1])
    }else{
      g =  abs(beta[i-1])+abs(beta[i])
    }
  
    mi = alpha[i] - g
    ma = alpha[i] + g
    if(mi<bmin){
      bmin = mi
    }
    if(ma>bmax){
      bmax = ma
    }
  }

  cat(c("min:",bmin,' max:',bmax,'\n') )
  tval = c(floor(bmin) :ceiling(bmax) )

  
  res= matrix(0,nrow=length(tval),ncol=length(alpha))
  k = 0
  for(i in tval){
    k= k+1
    val =sturm(alpha,beta,i)
    res[k,] = val
  }

  calcChanges = function(res){
    changes =sign(res)
    changes=cbind(rep(1,nrow(res)),changes)
    changes[changes==-1] =0
    int = changes[,2:ncol(changes)] - changes[,1:(ncol(changes)-1) ]
    if(nrow(res)>1){
      ck = rowSums(abs(int))
    }else{
      ck = sum(abs(int))  
    }
    ck
  }
  
  changes = calcChanges(res)
  
  res = cbind(tval,res)
  
  eigenvalues = c()
  ##find eigen values now
  for(i in 1:(nrow(res)-1) ){
    dif = abs(changes[i]-changes[i+1])
    if(dif >0){
      cat(c('found',dif,'eigenvalues in interval:',tval[i],'-',tval[i+1] ,'\n') )
      ##for each eigenvalue found
      changeLow =changes[i]
      changeHigh =changes[i+1]
      for( e in 1:dif ){
        ##half interval...
        lb =tval[i]
        ub = tval[i+1]
        
        for( w in 1:20){
          ns = (lb+ub)/2
          val =sturm(alpha,beta,ns)
          tc = calcChanges(matrix(val,nrow=1) )
          if(tc==changeLow){ ##set lower bound
            #res[i,1]= ns
            lb= ns
          }else if(tc==changeHigh){
            #res[i+1,1]= ns
            ub= ns
          }else{
          #  print('help with bound !') ##most prpably 2 eigenvalues in this bound
            ##we set upper boudn and hope for the best
            ub= ns
            ##and we set the next iteration lower bound and hope for the best
            tval[i] = ns
            changeHigh = tc
           # print('ss')
           # stop()
          }
          cat(c('eigenvalues in interval:',lb,'-',ub ,'\n') )
          if(abs(lb-ub)<0.001){
            break
          }
        }
        ##apend eigenvalue....
        eigenvalues= c(eigenvalues,(lb+ub)/2 )
        changeLow = changeHigh
        changeHigh = changeHigh-1
        
      }
    }
    
    
  }
  eigenvalues
}

res = getTriEigen(c(2,2,2),c(-1,-1))
###power

getEigen=function(A,b){
  last = 0
  for(i in 1:10000){
    tmp = A%*%b
    norm = as.numeric(sqrt(t(tmp)%*%tmp))
    b = tmp / norm
    if(abs(last-norm)<1e-10 ){
      print(i);
      print('conv\n')
      break;
    }
    last= norm
  }
  list('value'= norm,'vec'=b )
}

getEigenKern=function(A,b,orders,sigma,weights){
  last = 0
  for(i in 1:10000){
    tmp = fastMultiply(A,b,orders,sigma,weights) #A%*%p;
    print(tmp)
    norm = as.numeric(sqrt(t(tmp)%*%tmp))
    print(norm)
    b = tmp / norm
    print(b)
    if(abs(last-norm)<1e-10 ){
      print(i);
      print('conv\n')
      break;
    }
    last= norm
  }
  list('value'= norm,'vec'=b )
}
