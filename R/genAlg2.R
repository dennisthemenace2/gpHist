

genAlgSingle <- setRefClass("genAlgSingle",
                      fields = list( generations = 'numeric', population = 'numeric',nparams = 'numeric' ,P='list',lowerBound='numeric',upperBound='numeric',fnc='function',crossrate='numeric', crossidx='numeric', mutationrate='numeric', mutationidx='numeric'),
                      methods = list(
                        pow=function(a,b){
                          a^b
                        },
                        initialize = function(fn,gen, pop,npar,lb, ub,crossrate, crossidx, mutationrate, mutationidx){
                          .self$generations <<- gen
                          .self$population <<- pop
                          .self$nparams <<- npar
                          .self$lowerBound <<- lb
                          .self$upperBound <<- ub
                          .self$fnc <<- fn;
                          .self$crossrate <<- crossrate
                          .self$crossidx <<- crossidx
                          .self$mutationrate <<- mutationrate
                          .self$mutationidx <<- mutationidx
                        },
                        evaluate = function(pop){
                          for(i in 1:length(pop)){
                            pop[[i]]$err = .self$fnc( pop[[i]]$params)
                          } 
                          pop
                        },
                        
                        initPopulation = function(){
                          .self$P= list()
                          
                          for(i in 1:.self$population){
                            np = list()
                            
                            np$params = array(0,dim=c(.self$nparams))
                            for(p in 1:.self$nparams){ ##slow do in one call
                              np$params = runif(1,.self$lowerBound[p],.self$upperBound[p])
                            }
                            
                            .self$P[[i]] = np 
                          }
                          
                        },
                        selectParents = function(){
                          ###this shuld be canged in the future
                          a = sample(1:population,population ,replace=T )
                          b = sample(1:population,population ,replace=T )
                          selected =rep(0,population)
                          
                          for( pop in 1:population ){
                            if( .self$P[[a[pop] ]]$err < .self$P[[b[pop] ]]$err ){
                              selected[pop] =  a[pop]
                            }else{
                              selected[pop] =  b[pop]
                            }
                          }
                          selected
                        },
                        reproduce = function(selected){
                          
                          childrens = list()
                          
                          while(length(childrens)<population){
                            idx = sample(1:length(selected),2,replace=T )
                            
                            p = .self$P[[ selected[idx[1]] ]]
                            q = .self$P[[ selected[idx[2]] ]]
                            
                            newparams = crossover( p$params, q$params  )
                            for(i in 1:2){
                              child = list()
                              child$params = newparams[[i]]
                              child$params = mutation(child$params)
                              childrens[[length(childrens)+1]] = child;
                            }
                            
                          }
                          
                          childrens
                          
                        },
                        crossover = function(p,q){
                          if(runif(1,0,1) < crossrate){ 
                            return (list(p,q) )
                          }
                          params1 = rep(0,length(p))
                          params2 = rep(0,length(p))
                          
                          y1=0
                          y2=0;
                          betaq =0;
                          for( i in 1:length(p)){
                            if( runif(1,0,1) <= 0.5 ){
                              if ( p[i] < q[i]) {
                                y1 = p[i];
                                y2 = q[i];
                              } else {
                                y1 = q[i];
                                y2 = p[i];
                              }
                              yl =  .self$lowerBound[i]; 
                              yu = .self$upperBound[i];
                              randv = runif(1,0,1);
                              beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                              alpha = 2.0 - pow(beta,-(.self$crossidx +1.0));
                              if (randv <= (1.0/alpha)){
                                betaq = pow((randv*alpha),(1.0/(.self$crossidx +1.0)));
                              }else{
                                betaq = pow((1.0/(2.0 - randv*alpha)),(1.0/(.self$crossidx +1.0)));
                              }
                              c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                              beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                              alpha = 2.0 - pow(beta,-(.self$crossidx +1.0));
                              if (randv <= (1.0/alpha)){
                                betaq = pow ((randv*alpha),(1.0/(.self$crossidx +1.0)));
                              }else{
                                betaq = pow ((1.0/(2.0 - randv*alpha)),(1.0/(.self$crossidx +1.0)));
                              }
                              c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                              
                              if (c1 < yl) c1=yl;
                              if (c2 < yl) c2=yl;
                              if (c1 > yu) c1=yu;
                              if (c2 > yu) c2=yu;
                              if ( runif(1,0,1) <= 0.5) {
                                params1[i] = c2;
                                params2[i] = c1;
                              } else  {
                                params1[i] = c1;
                                params2[i] = c2;
                              }
                              
                            }else{
                              params1[i] =  p[i];
                              params2[i] =  q[i];
                            }
                            
                            
                          }
                          
                          list(params1, params2)
                          
                        },
                        mutation =function(p){
                          val=0
                          xy=0
                          deltaq=0
                          
                          for( i in 1:length(p)){
                            
                            if( runif(1,0,1) < mutationrate){
                              y = p[i];
                              
                              yl = .self$lowerBound[i]; 
                              yu = .self$upperBound[i];
                              delta1 = (y-yl)/(yu-yl);
                              delta2 = (yu-y)/(yu-yl);
                              randv =  runif(1,0,1) ;
                              
                              mut_pow = 1.0/(mutationidx+1.0);
                              if(randv <= 0.5){
                                xy = 1.0-delta1;
                                val = 2.0*randv+(1.0-2.0*randv)*(pow(xy,(mutationidx+1.0)));
                                deltaq =  pow(val,mut_pow) - 1.0;
                              }else{
                                xy = 1.0-delta2;
                                val = 2.0*(1.0-randv)+2.0*(randv-0.5)*(pow(xy,(mutationidx+1.0)));
                                deltaq = 1.0 - (pow(val,mut_pow));
                              }
                              y = y + deltaq*(yu-yl);
                              if (y < yl){ 
                                y = yl;
                              }
                              if (y > yu){
                                y = yu;
                              }
                              p[i] = y;
                              
                            }
                            #else skip
                          }
                          
                          p
                        },
                        clearPopulation =function(){
                          ###remove half of the population...
                          errs = rep(0,length(P))
                          for(i in 1:length(P)){
                            errs[i] = .self$P[[i]]$err
                          }
                          
                          ord = order(errs)
                          .self$P = .self$P[-ord[(length(P)/2 +1):length(P)] ]
                          
                        },
                        run = function(){
                          
                          .self$initPopulation()
                          .self$P = evaluate(.self$P)
                          
                          selected = selectParents();
                          children = reproduce(selected)
                          children= evaluate(children)
                          
                          
                          for(i in 1:.self$generations){
                            .self$P = append( .self$P , children)
                            
                            selected = selectParents();
                            children = reproduce(selected)
                            children= evaluate(children)
                            clearPopulation()
                            
                          }
                          
                          .self$P = append( .self$P , children)
                          .self$P
                        }
                        
                      ))


fctest = function(p){
  # ret = (p^2)+2
  ret = sin(2*pi*p*0.1)-2 #+p +sin(2*pi*p*0.7) 
}


alg = genAlgSingle(fn=fctest,gen=100, pop=10,npar=1,lb=c(0), ub=c(10),crossrate=0.7, crossidx=5, mutationrate=0.1, mutationidx=5)

res = alg$run()

points = seq(0,10,0.1)
plot(points,fctest(points))

for(i in 1:length(res)){
  points(res[[i]]$params, res[[i]]$err,col='red')
}