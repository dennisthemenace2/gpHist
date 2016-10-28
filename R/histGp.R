

options(error = recover)

options(warnings = recover)


fastMultiplyOne=function(X,Y){
  res = matrix(0,nrow=nrow(X),ncol=1)
  NR = nrow(X)
  
  prev = X[1:NR] * Y
  res[NR] = sum(prev)

  for(i in (NR-1):1){
    int = sum(prev[1:i])
    res[i] = int +sum(X[i]*Y[ (i+1):NR] )
  }
  
  res
}
###only works for standard min kernel.....
fastMultiply=function(X,Y,orders, sigma,weights){
  res = matrix(0,nrow=nrow(X),ncol=1)

  for(i in 1:ncol(X)){
    ord = orders[[i]]#order(X[,i] )
    Xt= as.matrix(X[ord,i])
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
      print(i)
      r=b-fastMultiply(A,x,orders,sigma,weights)#A%*%x;
      
    }else{
      r=r-alpha*Ap;
    }
    
    
    rsnew=t(r)%*%r;
    if(sum(sqrt(rsnew)) <1e-10){
      cat(c('break:', i, ' err:',sum(sqrt(rsnew)),'\n') )
      break;
    }
    p=r+as.numeric(rsnew/rsold)*p;
    rsold=rsnew;
  }
  
  x
}



downhillsimplex = function(fc,N,it=1000,tol=1e-10){
  
  #######
  alpha  = 1 # =1 (Reflexion)
  gamma = 2 #= 2 (Expansion)
  beta = 0.5 #1/2 (Kontraktion)
  sigma = 0.5 #Komprimierung
  
  if(N==1){
    print('this will not work for 1 parameter only')
  }
  
  bp = matrix(runif( (N+1)*N,0,5),nrow=N+1,ncol=N )
  
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
    np = fc(r)
    if(np<bp[1,N+1]){ ##better
      e =(1+gamma)*m -gamma*bp[N+1,1:N]
   #   print('querry e')
  #    print(e)
      
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
    
    nc = fc(cp)
    if(nc< bp[N+1,N+1]){ ##c is better
      bp[N+1,1:(N+1) ] = c(cp,nc)
      next;
    }
    for(k in 2:(N+1) ){ ##shrink
      bp[k,1:N] = bp[1,1:N] + sigma*(bp[k,1:N] - bp[1,1:N] )
      ##eval
      bp[k,N+1] = fc(bp[k,1:N])
    }
    
  }
  cat(c('iterations:',i,'\n') )
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
    res = c(res,w[i]* min(x[i],y[i]))
  } 
  sum(res)
}


####make some hist efficent gp
gphist=function(X,Y,sigma,minkern,weights, x_pred){
  N = nrow(X)
  P = ncol(X)
  cat(c('samples:',N, ' features:',P ,'\n' ) )
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
  orders = list()
  for(i in 1:ncol(X)){
    orders[[i]] = order(X[,i])  
  }
  
  # K  = outer(X, X, tmpkernel) + noise
  
  #K = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(X[i,],X[j,]) ) )
  K = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(X[i,],X[j,],weights) ) ) ##with params
  
  ###test kernel
  if( any( abs(K%*%Y- fastMultiply(X,Y,orders,sigma,weights)>1e-10 ) ) ){
    print('error multiplication missmatch')
    stop()
  }
  
  K = K +noise
  
  ####for testing
  if(F){
    ord =order(X[,1] )
    Xt= as.matrix(X[ord,1])
    K1 = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(Xt[i],Xt[j]) ) )
    
    ord =order(X[,2] )
    Xt= as.matrix(X[ord,2])
    K2 = outer(1:(N), 1:(N), FUN =Vectorize( function(i,j) minkern(Xt[i],Xt[j]) ) )
  }
  
  ####
  
  K_inv = solve(K)
  # K_inv = conjInv(K,matrix(0,nrow=nrow(K),ncol=ncol(K) ) )
  #  sum(abs(solve(K) - K_inv) )
  alpha = K_inv %*% Y
  ##no stuff
  
  ###calc alpha directly 
  alpha_direct = conjgradFP(X,Y,matrix(0,nrow=nrow(Y)),orders,sigma,weights)  # conjgrad(K,Y,matrix(0,nrow=nrow(Y)) )
  
  if(any(abs(alpha_direct-alpha)>1e-6 )){
    print('error alpha missmatch')
    stop()
  }
  
  k_star =  outer(1:(N), 1:nrow(x_pred), FUN =Vectorize( function(i,j) minkern(X[i,],x_pred[j,],weights) ) )
  predictions = t(k_star)  %*% alpha
  ### find index
  
  
  pred = c()
  for(i in 1:nrow(x_pred) ){
    ##each dim
    C =0
    for(d in 1:ncol(x_pred)){
      ord = orders[[d]]#order(X[,d])
      tmp = X[ord,d]
      
      idx= which(tmp <x_pred[i,d]) 
      if(length(idx)>0){
        idx = max(idx) 
        A = alpha_direct[ord][1:idx] %*% tmp[1:idx]
      }else{
        A=0
        idx = 0;
      }
      B = sum(alpha_direct[ord][(idx+1):nrow(alpha_direct)]) * x_pred[i,d]
      C = C+ weights[d]*(A+B)
    }
    pred = c(pred,C)
  }
  
  if(any(abs(pred - predictions)>1e-9 )) {
    print('error prediction differ')
    stop()
  }
  
  L =chol(K)
  sL = solve(L)
  v = t(sL) %*% k_star
  # V =outer(x_pred, x_pred,kern) - t(v)%*%v
  K_starStar = outer(1:(nrow(x_pred)), 1:nrow(x_pred), FUN =Vectorize( function(i,j) minkern(x_pred[i,],x_pred[j,],weights) ) )
  
  V =outer(1:(nrow(x_pred)), 1:nrow(x_pred), FUN =Vectorize( function(i,j) minkern(x_pred[i,],x_pred[j,],weights) ) ) - t(v)%*%v
  

  mu1 = N*sigma+   sum(X%*%weights) #sum(diag(K)) #### 
  if(abs(sum(diag(K)) -mu1)>1e-9){
    print('error calculating mu1')
    stop()
  }
  b = matrix(runif( N,0.5,1))
  TT= lanczosKERN(X,b,orders,sigma,weights )
  
  lambdas = eigen(TT, only.values = TRUE)
  beta= lambdas$values[1]
  mu2 =  sum(lambdas$values[1]^2)
  #mu2 =  sum(lambdas$values^2)
  #mu2 =  sum(K*K)
  
  ###check real quick
  if(abs(prod(lambdas$values)-det(K))>1e-10 ){
    print('error eigen values not matching')
    if(abs(lambdas$values[1]-eigen(K)$values[1])>1e-8 ){
      print('first eigenvalue not matching')
      stop()      
    }
    
  }
  
 
  t_bar = (beta*mu1 - mu2) /(beta*N - mu1)
  t_bar2 = t_bar^(2)
  
  beta2 = beta^2
  
  t1 = matrix(c(log(beta),log(t_bar)),nrow=1,ncol=2)
  #t2 = matrix(c(beta,beta^2,t_bar,t_bar^(2) ),nrow=2,ncol=2 )
  
  t2 =(1/(beta*t_bar2 - t_bar*beta2 ))* matrix(c(t_bar2,-beta2,-t_bar,beta ),nrow=2,ncol=2 ) ## write inverse directly
  
  t3 = matrix(c(mu1,mu2),ncol=1,nrow =2)
    
#  logdetM =  t1%*%solve(t2) %*%t3
  logdetM =  t1%*%t2%*%t3
  
  
  ##calc sigma
  sum1 = 0
  vecsum= 0
  #  if(length(lambdas$values)>1){
  #    for(i in length(lambdas$values):2){
  #      sum1 =sum1 +(1/lambdas$values[i]) * sum( lambdas$vectors[,i]^2)
  #      vecsum = vecsum +  sum( lambdas$vectors[,i]^2)
  #    }
  #  }
  
  k_star2 = k_star^2# %*% K
  k_dist = sum(k_star2) 
  ## cal k_dist
  
  k_dist_fast = 0
  for(i in 1:nrow(x_pred) ){
    ##each dim
    C =0
    for(d in 1:ncol(x_pred)){
      ord = orders[[d]]
      tmp = X[ord,d]
      idx= which(tmp <x_pred[i,d] ) 
      if(length(idx)>0){
        idx = max(idx)
        A = sum( (weights[d]*tmp[1:idx])^2)
      }else{
        A=0
        idx=0
      }
      B = (N-(idx) )*( (weights[d]*x_pred[i,d]) ^2)
      C = C+ (A+B)
    }
    k_dist_fast = k_dist_fast +C
  }
  
  if(k_dist_fast> k_dist){
    print('error')
    stop()
  }
  
  ###
  sum2 =  (1/lambdas$values[1]) *(k_dist  - vecsum)
  
  K_starStarFast = rowSums(x_pred) ### this is only the diagonal
  var = K_starStarFast - sum1 -sum2 + sigma
  
  
 # logmarginal  = -0.5*t(Y) %*%alpha - sum(log(diag(L))) #- nrow(X)/2*log(2*pi)
  logmarginal  = 0.5*t(Y) %*%alpha + sum(log(diag(L))) #- nrow(X)/2*log(2*pi)
  
#  logmarginal2  = t(Y) %*%alpha + sum(log(diag(L)))
  
  approxmarginal= 0.5*t(Y)%*%alpha_direct + logdetM/2
  
  cat(c('marginal approx dif:',approxmarginal - logmarginal , '\n'))
  ##get log det bound
  
  ##optimze ?
 
  list('predictions'=pred,'var'=var,'logmarginal'=approxmarginal)
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
  Y = matrix(c(X[2]+X[1]^2*0.4+3+rnorm(10,0,0.1) ) )
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
