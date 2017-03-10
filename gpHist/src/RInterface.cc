
#include "RMat.h"
#include "R.h"     // R functions
#include "Rmath.h" // Rmath

#include <R_ext/Utils.h>
#include <Rinternals.h>


#include <iostream> // For cout, cerr, endl, and flush
#include <algorithm>
#include <stdlib.h>
#include <float.h>

#include <pthread.h>


void Cpplanczos(RMat& Xmat, RMat& bmat,RMat& vMatOrders,int k,double*alphas,double*betas,double sigma);

// calculate sturm seq
void sturmSeq(double* alpha,int nalpha, double *beta,double lambda,double* result){
  
  result[0] = lambda- alpha[0];
  result[1] = (lambda- alpha[1])*result[0] - ( (*beta) * (*beta) ) ;
  for(unsigned int i=2;i<nalpha;++i){
    result[i] = (lambda- alpha[i])*result[i-1] - (beta[i-1]*beta[i-1]) *result[i-2];
  }
}
//estimate number of changes in sturm seq
int calcSeqChanges(RMat& m1){
  //first make signs..
  double *ptr = m1.getValuesPtr();
  int past =1;
  int result = 0;
  for(unsigned int i =0;i<m1.NumRows();++i){
    if(*ptr>=0){
      if(past==0){ // is not positive
        ++result;      
      }
      past = 1;
    }else{
      if(past==1){ 
        ++result;      
      }
      past = 0;
    }
    ++ptr;
  }
  
  return result;
}
// estimates the eigenvalues of a tri diagonal matrix using sturm seq.
void getEigenSturm(double *alpha,int nalpha,double*beta,int k,double *result){
  
  if(k==0 || k>nalpha){
    std::cout<<"invalid number of eigenvalues. k="<<k<<" nalpha="<< nalpha<< "\n";
    return;
  } 
  //estimate bounderies
  double bmin=DBL_MAX;
  double bmax=-DBL_MAX;
  
  for(unsigned int i=0;i<nalpha;++i){
    double g;
    if(i == 0){
      g = fabs(beta[i]);
    }else if(i == nalpha-1){
      g = fabs(beta[i-1]);
    }else{
      g =  fabs(beta[i-1])+fabs(beta[i]);
    }
    
    double min = alpha[i] - g;
    double max = alpha[i] + g;
    if(min<bmin){
      bmin = min;
    }
    if(max>bmax){
      bmax = max;
    }
  }// end for estimate boundaries
  
  std::cout<<"search in interval:"<<bmin<<"-"<<bmax<<std::endl;
  // void sturmSeq(double* alpha,int nalpha, double *beta,double lambda,double* result){
  RMat seqhigh(nalpha,1);
  sturmSeq(alpha,nalpha,beta,bmax,seqhigh.getValuesPtr()  );
 // seqhigh.Print("seqHigh\n");
  
  RMat seqlow(nalpha,1);
  sturmSeq(alpha,nalpha,beta,bmin,seqlow.getValuesPtr()  );
 // seqlow.Print("seqlow\n");
  
  int changeslow_first = calcSeqChanges(seqlow);
  int changeslow= changeslow_first;
  
  int changeshigh = calcSeqChanges(seqhigh);
  int dif = abs(changeshigh -changeslow );
  std::cout<<"found "<<dif<<"eigenvalues"<<std::endl;
  
  std::cout<<"changes high "<<changeshigh<<" chnages low "<< changeslow<<std::endl;
  
  for(unsigned int e=0; e<k; ++e ){
    double lb = bmin;
    double ub = bmax;
    std::cout<<"interval:"<<lb<<"-"<<ub<<std::endl;
    std::cout<<"changes high "<<changeshigh<<" chnages low "<< changeslow<<std::endl;
    
    for(unsigned int w=0;w<10000; ++w){
      double ns = (lb+ub)/2.0;
      sturmSeq(alpha,nalpha,beta,ns,seqhigh.getValuesPtr()  );
      int tc = calcSeqChanges(seqhigh);
      if(tc==changeslow){ //##set lower bound
        lb= ns;
      }else if(tc==changeshigh){
        ub= ns;
      }else{
        bmax =ns;
        lb = ns;
        changeslow = tc;
      }
      
    //  std::cout<<"eigenvalues in interval"<<lb<<"-"<<ub<<"\n";
      if(fabs(lb-ub)<0.0001){
        std::cout<<"iterations:"<<w<<"\n";
        break;
      }
    }
    result[e]= (lb+ub)/2.0;
    std::cout<<"eigenvalues adding:"<<result[e]<<"\n";
    changeshigh = changeslow;
    changeslow = changeslow_first;
    
  }
  
}

/*
void fastMultiplyBACK(const RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){
  
//  RMat prev(Y.NumRows ,1);
  result.SetZero();
  unsigned int NR = X.NumRows();
    
  for(unsigned int i=1;i<=X.NumCols();++i){
  
    unsigned int idx = 0;
    double matcount = 0;
    double Ysum;
    
    for(unsigned int k=1;k<NR;++k ){ // estimate prev
      idx = orders(k,i);
      
      Ysum= 0;
      for(unsigned int x=k;x<=NR;++x ){ // estimate Y.sum(k,NR); 
        Ysum += Y(orders(x,i) ,1);
      }
      
      result(idx,1) +=  matcount + X(idx, i)*Ysum;
      matcount = matcount + X(idx,i) * Y(idx,1);

    }
    idx = orders(NR,i);
    
    result(idx,1) += matcount + X(idx,i) * Y(idx,1);
      
  }
 // for(unsigned i =1;i<=Y.NumRows();++i){
  //  result(i,1)  += sigma*Y(i,1);
//  }
  result += sigma*Y ;
  
}*/

///// Chose of parallel or not
#define COMPILE_PARALLEL 0 // set to 0 for no parallel execution

#if COMPILE_PARALLEL == 0
  
  void fastMultiply( RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){
    //std::cout<<"fastMultiply single"<<std::endl;
    
    
    //  RMat prev(Y.NumRows ,1);
    result.SetZero();
    unsigned int NR = X.NumRows();
    
    double *result_ptr = result.getValuesPtr();
    
    RMat YsumMat(NR,1);
    double*  YsumMat_ptr = YsumMat.getValuesPtr();
    double* Y_ptr = Y.getValuesPtr();
     
    for(unsigned int i=1;i<=X.NumCols();++i){
      
      unsigned int idx = 0;
      double matcount = 0;
  
      
      double* order_prt = orders.getColPtr(i);
      double* idx_prt = order_prt;
      double* X_ptr = X.getColPtr(i);
      
      
      //calc Y_sum first
      double Ysum=0;
      for(int x=NR-1;x>=0;--x ){ // estimate Y.sum(k,NR); 
        Ysum += Y_ptr[(int) order_prt[x]-1 ]; //Y(orders(x,i) ,1);
        YsumMat_ptr[x] = Ysum;
      }
      //YsumMat.Print("Ysummat");
      
      for(unsigned int k=0;k<NR-1;++k ){ // estimate prev
        idx =  *(idx_prt++); //orders(k,i);
        idx -=1;
        
        result_ptr[idx] +=  matcount + X_ptr[idx]*YsumMat_ptr[k]; // Ysum
        matcount = matcount + X_ptr[idx] * Y_ptr[idx];
        
      }
      idx = *(idx_prt++);//orders(NR,i);
      idx -=1; 
    //  std::cout<<"Y(idx,1)"<<Y(idx,1)<<std::endl;
      
      result_ptr[idx] += matcount + X_ptr[idx] * YsumMat_ptr[NR-1];
      
    }
    result += sigma*Y ;
    
    //std::cout<<"fastMultiply end"<<std::endl;
  }

#elif COMPILE_PARALLEL ==1
   /// PARALLEL Execution 
  #define NUM_THREADS 4  // NUMBER OF THREADS THAT ARE CREATED
   
  struct thread_data_shared{
    RMat* X;
    RMat* Y;
    RMat* orders;
  };
  
  struct thread_data{
    unsigned  int  c;//start
    unsigned int C;//end row
    thread_data_shared* shared;
  };
  
/*  struct ret_thread_data{
  //  double *result;
};*/
  
  void *threadFn(void *arg){
    
    struct thread_data *data;
    data = (struct thread_data *) arg;
    unsigned int NR = data->shared->X->NumRows();
    
 //   ret_thread_data *ret = new ret_thread_data;
    double *result = new double[NR];
    memset(result,0,sizeof(double)*NR);
    ////////
    RMat YsumMat(NR,1);
    double*  YsumMat_ptr = YsumMat.getValuesPtr();
    double* Y_ptr = data->shared->Y->getValuesPtr();
    // do right collumns
    for(unsigned int i=data->c+1; i<= data->C ;++i){
      
      unsigned int idx = 0;
      double matcount = 0;
      
      double* order_prt = data->shared->orders->getColPtr(i);
      double* idx_prt = order_prt;
      double* X_ptr = data->shared->X->getColPtr(i);
      
      //calc Y_sum first
      double Ysum=0;
      for(int x=NR-1;x>=0;--x ){ // estimate Y.sum(k,NR); 
        Ysum += Y_ptr[(int) order_prt[x]-1 ]; //Y(orders(x,i) ,1);
        YsumMat_ptr[x] = Ysum;
      }
      //YsumMat.Print("Ysummat");
      
      for(unsigned int k=0;k<NR-1;++k ){ // estimate prev
        idx =  *(idx_prt++); //orders(k,i);
        idx -=1;
        
        result[idx] +=  matcount + X_ptr[idx]*YsumMat_ptr[k]; // Ysum
        matcount = matcount + X_ptr[idx] * Y_ptr[idx];
        
      }
      idx = *(idx_prt++);//orders(NR,i);
      idx -=1; 
      //  std::cout<<"Y(idx,1)"<<Y(idx,1)<<std::endl;
      
      result[idx] += matcount + X_ptr[idx] * YsumMat_ptr[NR-1];
    
    }
    
    
    
    
    //////////
    /*
    
    for(unsigned i =data->c+1; i<= data->C ;++i){
      ///
      unsigned int idx = 0;
      double matcount = 0;
      double Ysum;
      
      for(unsigned int k=1;k<NR;++k ){ // estimate prev
        idx = (*data->shared->orders) (k,i);
        
        Ysum= 0;
        for(unsigned int x=k;x<=NR;++x ){ // estimate Y.sum(k,NR); 
          Ysum += (*data->shared->Y)( (*data->shared->orders)(x,i) ,1);
        }
        
        ret->result[idx-1] +=  matcount + (*data->shared->X)(idx, i)*Ysum;
        matcount = matcount + (*data->shared->X)(idx,i) * (*data->shared->Y)(idx,1);
        
      }
      idx = (*data->shared->orders)(NR,i);
      
      ret->result[idx-1] += matcount + (*data->shared->X)(idx,i) * (*data->shared->Y)(idx,1);
      
    }*/
    
    pthread_exit( (void *)result);
  }

  
  void fastMultiply( RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){
  //  std::cout<<"fastMultiply parallel"<<std::endl;
    
    //  RMat prev(Y.NumRows ,1);
    result.SetZero();
    unsigned int NR = X.NumRows();
    
    
    unsigned int numThreads = X.NumCols(); // create a thread for each dimension
    if(numThreads>NUM_THREADS){
      numThreads = NUM_THREADS;
    }
    unsigned int each = X.NumCols()/numThreads;
    
   // std::cout<<"create "<<numThreads<<" each can do:"<<each <<std::endl;
    
    thread_data *tdata = new thread_data[numThreads];
    thread_data_shared *shared = new thread_data_shared;
    
    //create shared data
    shared->X = &X;
    shared->Y = &Y;
    shared->orders = &orders;
    // shared data
    
    pthread_t *threads= new  pthread_t[numThreads];
    
    // set up data....
    for(unsigned int c=0;c<numThreads;++c){
      tdata[c].c=c*each;
      tdata[c].C=c*each + each;
      if(c== numThreads-1){
        tdata[c].C=X.NumCols();
      }
      tdata[c].shared = shared;
    }
    // set up threat attr.
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
    // start threads
    for(unsigned t=0; t<numThreads; ++t)  {            
      pthread_create(&threads[t],&attr, threadFn, (void *)&tdata[t]);
    }
  
    // wait for all threads
    for(unsigned t=0; t<numThreads; ++t)  {
      //ret_thread_data *ret;
      double *ret;
      pthread_join(threads[t], (void**) &ret);
  
      RMat resMat(ret,NR ,1);
      result +=resMat; // add dimensions
      
      //delete [] ret->result;
      delete []ret;
    }
    
    pthread_attr_destroy(&attr);
    
    delete [] threads;
    delete [] tdata;
    delete shared;
    
    result += sigma*Y ;
  
  }
#else
  
struct thread_data_shred{
  unsigned int initcnt;
  pthread_mutex_t initmutex;
  
  RMat* X;
  RMat* Y;
  RMat *result;
  RMat *orders;
  
  bool run;
  
};
  
struct thread_data_Mutex{
    unsigned  int  i;
    unsigned  int  c;
    unsigned int C;
  
    pthread_mutex_t xmutex;
    
    pthread_mutex_t pmutex;
    pthread_mutex_t tmutex;
    pthread_mutex_t tsmutex;

    
    thread_data_shred *shared;
    
};

void *clientFnMutex(void *arg){
    
    struct thread_data_Mutex *data;
    data = (struct thread_data_Mutex *) arg;
    //   cout<<"run thread"<<endl;
    
    //
    pthread_mutex_lock (&data->shared->initmutex);
    //   cout<<"get init lock"<<endl;
    pthread_mutex_lock (&data->pmutex);       
    data->shared->initcnt +=1; 
    pthread_mutex_unlock (&data->shared->initmutex);
    

    //cout<<"thread run code"<< data->c <<"-"<<data->C  <<endl;
    while(true){
      
      pthread_mutex_lock (&data->tmutex); // wait till i shall run
      // check abort
      if(!data->shared->run){ // leave
        break;
      }
      
   
      unsigned int NR = data->shared->X->NumRows();

      double *result = new double[NR];
      memset(result,0,sizeof(double)*NR);
      ////////
      RMat YsumMat(NR,1);
      double*  YsumMat_ptr = YsumMat.getValuesPtr();
      double* Y_ptr = data->shared->Y->getValuesPtr();
      // do right collumns
      for(unsigned int i=data->c+1; i<= data->C ;++i){
        
        unsigned int idx = 0;
        double matcount = 0;
        
        double* order_prt = data->shared->orders->getColPtr(i);
        double* idx_prt = order_prt;
        double* X_ptr = data->shared->X->getColPtr(i);
        
        //calc Y_sum first
        double Ysum=0;
        for(int x=NR-1;x>=0;--x ){ // estimate Y.sum(k,NR); 
          Ysum += Y_ptr[(int) order_prt[x]-1 ]; //Y(orders(x,i) ,1);
          YsumMat_ptr[x] = Ysum;
        }
        //YsumMat.Print("Ysummat");
        
        for(unsigned int k=0;k<NR-1;++k ){ // estimate prev
          idx =  *(idx_prt++); //orders(k,i);
          idx -=1;
          
          result[idx] +=  matcount + X_ptr[idx]*YsumMat_ptr[k]; // Ysum
          matcount = matcount + X_ptr[idx] * Y_ptr[idx];
          
        }
        idx = *(idx_prt++);//orders(NR,i);
        idx -=1; 
        //  std::cout<<"Y(idx,1)"<<Y(idx,1)<<std::endl;
        
        result[idx] += matcount + X_ptr[idx] * YsumMat_ptr[NR-1];
        
      }
      
      
      pthread_mutex_lock (&data->shared->initmutex);
      double *fres = data->shared->result->getValuesPtr();
      for(unsigned int i=0;i<NR;++i){
        fres[i] += result[i];
      }
      pthread_mutex_unlock (&data->shared->initmutex);
      
      delete [] result;
      //     cout<<"thread unlock stuff"<<i<<endl;
      
      pthread_mutex_unlock (&data->tmutex);
      pthread_mutex_lock(&data->xmutex);        
      pthread_mutex_unlock (&data->pmutex);   
      
      //   cout<<"thread sperre"<<i<<endl;
      
      pthread_mutex_lock (&data->tsmutex);
      pthread_mutex_lock (&data->pmutex);    
      pthread_mutex_unlock (&data->tsmutex);
      pthread_mutex_unlock(&data->xmutex); 
      
      // cout<<"thread did run:"<<i<<endl;
      
      
    }
    // cout<<"finished"<<endl;
    
    std::cout<<"ending thread"<<std::endl;
    pthread_exit(NULL);

}

thread_data_Mutex *tdata=0;  
pthread_t *threads=0;
pthread_attr_t attr;

#define NUM_THREADS 4
  
  
void CstartThreads(){
  if(threads != NULL){
    std::cout<<"Already initialized!"<<std::endl;
    return;
  }
  tdata = new thread_data_Mutex[NUM_THREADS];
  
  threads= new  pthread_t[NUM_THREADS];
  thread_data_shred *shared = new thread_data_shred;
  

    /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  
 // pthread_mutex_t initmutex;
  pthread_mutex_init(&shared->initmutex, NULL);
  
 // unsigned int initcnt = 0;
  
  shared->initcnt = 0;
  shared->run = true;
  //shared->initmutex = &initmutex;
  
  for(unsigned int c=0;c<NUM_THREADS;++c){
    tdata[c].shared = shared;
    
    
    pthread_mutex_init(&tdata[c].tmutex, NULL);
    pthread_mutex_init(&tdata[c].pmutex, NULL);
    pthread_mutex_init(&tdata[c].tsmutex, NULL);
    pthread_mutex_init(&tdata[c].xmutex, NULL);
    
    pthread_mutex_lock(&tdata[c].tmutex);
    
  }
  
  
  for(unsigned t=0; t<NUM_THREADS; ++t)  {
    pthread_create(&threads[t],&attr, clientFnMutex, (void *)&tdata[t]);
  }
  
  while(true){
    //  cout<<"waiting"<<endl;
    pthread_mutex_lock (&shared->initmutex);    
    if(shared->initcnt ==NUM_THREADS){
      pthread_mutex_unlock (&shared->initmutex);
      break;
    }
    pthread_mutex_unlock (&shared->initmutex);
  } // wait

  
  std::cout<<"init threads"<<std::endl;
}
void CstopThreads(){
  if(threads == NULL){
    std::cout<<"no threads initilized!"<<std::endl;
    return;
  }
  
  thread_data_shred *shared = tdata[0].shared;
  
  shared->run = false;
  
  // make them run
 // for(unsigned t=0; t<NUM_THREADS; ++t)  {
   // pthread_join(threads[t], (void**) &ret);
  //}
  for(unsigned int c=0;c<NUM_THREADS;c++){
    //  cout<<"lock sperre"<<i<<endl;
    pthread_mutex_unlock(&tdata[c].tmutex);
   // pthread_mutex_lock(&tdata[c].xmutex);  
  //  pthread_mutex_lock(&tdata[c].tsmutex); 
  }
    
  for(unsigned int c=0;c<NUM_THREADS;++c){
      pthread_mutex_destroy(&tdata[c].tmutex);
      pthread_mutex_destroy(&tdata[c].pmutex);
      pthread_mutex_destroy(&tdata[c].tsmutex);
      pthread_mutex_destroy(&tdata[c].xmutex);        
  }
    
  pthread_mutex_destroy(& shared->initmutex);
  
  pthread_attr_destroy(&attr);
  
  delete [] threads;
  delete shared;
    
  delete [] tdata;
  threads= NULL;
  tdata= NULL;
}
void fastMultiply( RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){
  
 
  
  //for(unsigned int c=0;c<NUM_THREADS;++c){
  //  tdata[c].shared = shared;
//  }
  
  result.SetZero();
  unsigned int NR = X.NumRows();
  
  //set up the threads
  thread_data_shred* shared = tdata[0].shared;
  shared->X = &X;
  shared->Y = &Y;
  shared->orders = &orders;
  shared->result = &result;
  
  
  unsigned int numThreads = X.NumCols(); // calculate the number of threads that we require.
  if(numThreads>NUM_THREADS){
    numThreads = NUM_THREADS;
  }
  unsigned int each = X.NumCols()/numThreads;

  // set up data....
  for(unsigned int c=0;c<numThreads;++c){
    tdata[c].c=c*each;
    tdata[c].C=c*each + each;
    if(c== numThreads-1){
      tdata[c].C=X.NumCols();
    }
    tdata[c].shared = shared;
  }
  
  // start the threads that we need.
  for(unsigned int c=0;c<numThreads;c++){
    //  cout<<"lock sperre"<<i<<endl;
    pthread_mutex_lock(&tdata[c].xmutex); 
    pthread_mutex_lock(&tdata[c].tsmutex); 
    //   cout<<"unlock tread"<<i<<endl;          
    pthread_mutex_unlock(&tdata[c].tmutex);
    pthread_mutex_unlock(&tdata[c].xmutex); 
  }
  
  //wait until they are all done...
  for(unsigned int c=0;c<numThreads;c++){         
    //    cout<<"get program lock"<<i<<endl;
    pthread_mutex_lock(&tdata[c].pmutex); 
    pthread_mutex_lock(&tdata[c].tmutex);
    //  cout<<"remove thread sperre"<<i<<endl;
    pthread_mutex_unlock(&tdata[c].pmutex);     
    pthread_mutex_unlock(&tdata[c].tsmutex);  
    //    cout<<"unlock prgramm sperre"<<i<<endl;
  }
  
  // add sigma
  result += sigma*Y ;
  
}
#endif


bool conjgradFP( RMat& A, RMat& b,RMat& x ,RMat& orders,double sigma, unsigned int max_iterations=1000 ) {

  RMat res( x.NumRows(),1 );

  fastMultiply(A,x,res,orders,sigma); //A%*%x;

//  std::cout<<"conjgradFP first "<<max_iterations <<std::endl;
  RMat r= b-res;// 
  RMat p=r;
  
  double rsold=r.SquareSum();//t(r)%*%r 
  double rsnew ;
    
  for (unsigned int i=1;i<=max_iterations;++i){
    fastMultiply(A,p,res,orders,sigma);// res =  A%*%p;
//double alpha= rsold / ( p.Transpose() *res);
//  std::cout<<alpha<<std::endl;
    double alpha= rsold / p.ScalarProd(res) ;

    x += alpha*p;
    
    if(i % 50 == 0){
      fastMultiply(A,x,res,orders,sigma);
      r=b-res; //#A%*%x;
  //    std::cout<<"update"<<std::endl;
  //    r.Print();
    }else{
      r -= alpha*res;
     // r.Print();
    }
    
    
    rsnew= r.SquareSum();//t(r)%*%r;
    if(sqrt(rsnew) <1e-10){
//#  cat(c('break:', i, ' err:',sum(sqrt(rsnew)),'\n') )
    //  std::cout<<"conjgradFP converged after "<<i<<std::endl;
      return true;
    }
    p*=(rsnew/rsold);
    p+=r;
    //p.Print();
    rsold=rsnew;
  }
    
  //std::cout<<"DIST:"<<  sqrt(rsnew)<<std::endl;
  return false;
}
// power method to get larges eigenvalue and vector
double getEigen( RMat& A,RMat &b,RMat& orders,double sigma){
  double last = 0;
  RMat res( b.NumRows(),1 );

  for(unsigned int i=0;i<2000;++i){ // or give it a max iterations parameter
    fastMultiply(A,b,res,orders,sigma); //A%*%b;
    double norm = sqrt(res.SquareSum()) ;
  //  res.Print();
   // std::cout<< norm<<std::endl;
    b = res;
    b /= norm;
  //  b.Print();
    if(fabs(last-norm)<1e-10 ){
      std::cout<<"converged:"<<i<<std::endl;  
        break;
    }
    last= norm;
  }
  return last;
}

//inverse iteration to get eigenvector for an eigenvalue
void inverseInteration( RMat& A,double lambda,RMat &result,RMat& orders,double sigma){
  unsigned int N = A.NumRows();
  result = 1.0/ double(N); // innit vector
 // result.Print("result after init:");
  double last = 0.0;
  RMat tmp( N,1);
  
  for(unsigned int i=0;i<1000;++i){
 //tmp = (A- I*lambda)%*%b
    //fastMultiply(A,result,tmp,orders,sigma-lambda); //A%*%b;
    conjgradFP( A, result, tmp , orders,sigma-lambda, 200 );
 //   tmp.Print("tmp:");
    double norm = sqrt(tmp.SquareSum());
    result = tmp;
    result /= norm;
 //   result.Print("result:");
    if(fabs(last-norm)<1e-10 ){
 //   print(i);
   //   print('conv\n')
      std::cout<<"converged:"<<i<<std::endl;        
        break;
    }
    last= norm;
  }

}



// two helper functions for predictions
double sumOrdered(RMat& alpha,RMat&orders,int start,int end,int col){
  double result = 0;
  double *orders_ptr = orders.getColPtr(col)+(start-1);
  double *alpha_ptr = alpha.getValuesPtr();
  
  for(int i= start;i<=end;++i){
    int idx = *(orders_ptr++);
    idx -=1; // or set pointer to invalid location and skip this ?
    result += alpha_ptr[idx];
  }
  return result;
}
  
double sumSquareOrdered(RMat& X,RMat&orders,int start,int end, int col){
  double result = 0;
  double *X_ptr = X.getColPtr(col);
  double *orders_ptr = orders.getColPtr(col)+(start-1);
  
  for(int i= start;i<=end;++i){
    int idx = *(orders_ptr++);
    idx -=1;
    //  std::cout<<"sumordidx"<<idx<<std::endl;
    double val = X_ptr[idx];
    result += val*val;
  }
  
  return result;
}
    
double sumOrdered(RMat& X,RMat& alpha,RMat&orders,int start,int end, int col){
  double result = 0;
  double *alpha_ptr = alpha.getValuesPtr();
  double *X_ptr = X.getColPtr(col);
  double *orders_ptr = orders.getColPtr(col)+(start-1);
  
  for(int i= start;i<=end;++i){
    int idx = *(orders_ptr++);
    idx -=1;
  //  std::cout<<"sumordidx"<<idx<<std::endl;
    result += X_ptr[idx]*alpha_ptr[idx];
  }
  
  return result;
}

int findIdx(RMat& X,RMat& orders,double value,int col){
 // try intersection...
  double *orders_ptr = orders.getColPtr(col); // pointer to column
  double *X_ptr = X.getColPtr(col); // pointer to column
  
  int min = 0;
  int max = X.NumRows()-1;
  int mid; 
//  std::cout<<"max"<<max<<std::endl;
//  std::cout<<"mid"<<mid<<std::endl;
  do{
    mid = (min+max)/2;
//    std::cout<<"mid"<<mid<<std::endl;
    
//    std::cout<<" X_ptr[(int) orders_ptr[mid] ] "<< X_ptr[(int) orders_ptr[mid] -1]<<std::endl;
    if(X_ptr[(int) orders_ptr[mid]-1 ] >value  ){ // go left
      max = mid;
//      std::cout<<"go left max "<<max<<std::endl;
    }else{ // go right should i consider equal 
      min = mid;
//      std::cout<<"go right min "<<min<<std::endl;
    }

  }while(max-min >1);

  if( value <= X_ptr[(int) orders_ptr[min]-1 ] ){
 //   std::cout<<"value smaller "<<std::endl;
    mid = min;//-1;
  }else{
    mid = max;
    if(value>X_ptr[(int) orders_ptr[max]-1 ]){
  //    std::cout<<"value bigger "<<std::endl;
      mid = max+1;
    }
  }
  
  return mid;
}

/*int findIdx(RMat& X,RMat& orders,double value,int col){
  int k =1;
  for(;k<=X.NumRows();++k){ // this is linear this suckz
    int idx = orders(k,col);
    if(value<X(idx,col)){
      break;
    }
  }
  
  return k-1;
}*/

// coarse approximation
void CppHistVarianceCoarse(double* result,int numRows,int numCols,int numRows2,int numCols2,double* mat1,double* mat2,double lambda,double sigma,double* orders){
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatPred(mat2,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
 // RMat  vMatAlpha(mat3,numRows,1); // not used remove from call
  RMat  vMatOrders(orders,numRows,numCols);
  
  for(int i=1; i<=numRows2;++i ){
      double A=0;
      double B=0;
      
    for(int d=1; d<= numCols2;++d){
      int idx = findIdx(vMatX,vMatOrders, vMatPred(i,d),d);
     // std::cout<<"idx="<<idx<<std::endl;
      
       if(idx>0){
           A = sumSquareOrdered(vMatX,vMatOrders,1,idx,d);
      //   std::cout<<"A="<<A<<std::endl;
        }else{
          A=0;
          idx=0;
        }
        
        B = (numRows-idx )* ( vMatPred(i,d)*vMatPred(i,d) );
        //std::cout<<"B="<<B<<std::endl;
     // C = C+ (A+B);
    }
    vMatResult(i,1) = A+B;
  }
//  vMatResult.Print("vmatresult");
  vMatResult *= (1.0/lambda);

  
 // vMatPred.RowSums().Print("rowssums") ;
  
  vMatResult = vMatPred.RowSums() -vMatResult;
 // vMatResult += sigma;

}

// finer approximation use all calculated lambdas
void CppHistVarianceFine(double* result,int numRows,int numCols,int numRows2,int numCols2,double* X,double* pred,double *lambda,int nlambda, double *vectors,double sigma,double* orders){
  RMat  vMatX(X,numRows,numCols);
  RMat  vMatPred(pred,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
 // RMat  vMatAlpha(alphas,numRows,1);
 // RMat  vMatOrders(orders,numRows,numCols); // not used...maybe use later on ?! 
  RMat  vMatVectors(vectors,numRows,nlambda);
  
  RMat  k_star(numRows,1); // number of alphas 
  
  vMatResult =   vMatPred.RowSums(); // init this one...
// vMatResult.Print("k_starstar");
  for(int s=1; s<=numRows2;++s ){ // each sample to predict
    //calc k_star
   //k_star =  outer(1:(N), 1:nrow(x_pred), FUN =Vectorize( function(i,j) GP$kernel(X[i,],x_pred[j,],GP$weights) ) )
    k_star.SetZero();
    //estimate k_star
    for(unsigned x=1;x<=numCols;++x){ // each dimension
      double *col = vMatX.getColPtr(x);
      double *starpointer =k_star.getValuesPtr();
      double value =  vMatPred(s,x); // check against this value
      
      for(unsigned y=1;y<=numRows;++y){ // each row in X
        
        if(value< *col){ // value smaller X
          *starpointer +=value;
        }else{
          *starpointer +=*col;
        }
        ++col; // next xvalue
        ++starpointer; // next k_star value
      }
    }
 
  //  k_star.Print("k_star");
    double secondTerm =0;
    double sumOfProjectionLengths = 0;
    
    for(unsigned int i=0; i<nlambda;++i ) {
      RMat  vector(vMatVectors.getColPtr(i+1), numRows,1); // get eigenvector
    //  vector.Print("eigenv:\n");
      double projectionLength =  vector.ScalarProd( k_star) ; //
      double projectionLength2= projectionLength*projectionLength;
    //  std::cout<<"projectionLength:"<<projectionLength<<std::endl;
      
      secondTerm += ( 1.0 / lambda[i] ) *projectionLength2;
      sumOfProjectionLengths +=  projectionLength2;
    }
   // std::cout<<"secondTerm:"<<secondTerm<<std::endl;
  //  std::cout<<"sumOfProjectionLengths:"<<sumOfProjectionLengths<<std::endl;
    
    // k_star2 = k_star^2# %*% K
    double k_dist = k_star.SquareSum(); 
 //   std::cout<<"k_dist:"<<k_dist<<std::endl;
   // std::cout<<"sumOfProjectionLengths:"<<sumOfProjectionLengths<<std::endl;
    
  //  std::cout<<"( 1.0 / lambda[nlambda-1] ):"<<( 1.0 / lambda[nlambda-1] )<<std::endl;
  //  secondTerm += ( 1.0 / lambda[nlambda-1] ) * ( k_dist - sumOfProjectionLengths );
    
    if(k_dist - sumOfProjectionLengths <0){
  //    std::cout<<"projection smaller zero"<<std::endl;
      std::cerr << "Attention: k_dist: " << k_dist << " - sumOfProjectionLengths: " <<  sumOfProjectionLengths << " is smaller than zero -possible ghost eigenvalues!" << std::endl;
    //  std::cout << "subtracting... "<< ( 1.0 / lambda[nlambda-1] ) * ( k_dist - sumOfProjectionLengths )<<"to the scnd term!";
    }
    //K_starStarFast = rowSums( x_pred %*%GP$weights)
      
   // std::cout<<"secondTerm:"<<secondTerm<<std::endl;
    //std::cout<<"k_star_star:"<<vMatResult(s,1)<<std::endl;
    vMatResult(s,1) -= secondTerm;
   // std::cout<<"result:"<<vMatResult(s,1)<<"\n"<<std::endl;
    
  }// end for each sample
  
}


void CppHistPredict(double *result,
                    int numRows, int numCols,int numRows2 ,int numCols2,
                    double *mat1, double *mat2,double* mat3,double *orders){
  
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatPred(mat2,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
  RMat  vMatAlpha(mat3,numRows,1);
  RMat  vMatOrders(orders,numRows,numCols);
  
  //NRA = nrow(GP$alpha)
    
  //  pred = c()
    for(int i=1; i<=numRows2;++i ){
      double C=0;
      double A=0;
      double B=0;
      
      for(int d=1; d<= numCols2;++d){
        // find idx....
        
        int idx = findIdx(vMatX,vMatOrders, vMatPred(i,d),d);
      //  std::cout<<"idx"<<idx<<std::endl;
        
        if(idx>0){
          if(idx==numRows){//###largest
            C = C+ sumOrdered(vMatX,vMatAlpha,vMatOrders,1,numRows,d);
        //    std::cout<<"C"<<C<<std::endl;
            continue;
          }else{
            A = sumOrdered(vMatX,vMatAlpha,vMatOrders,1,idx,d);
      //      std::cout<<"A"<<A<<std::endl;
          }
        }else{ //##point even further left....
          A=0;
        }
        B = sumOrdered(vMatAlpha,vMatOrders,idx+1,numRows,d) * vMatPred(i,d);
    //    std::cout<<"B"<<B<<std::endl;
        C = C+A+B;
      }
     
     vMatResult(i,1)=C;
    }
    
  
}

void CppHist(double *result,
                   int numRows, int numCols,int numRows2 ,int numCols2,
                   double *mat1, double *mat2,double sigma , double *orders, double* logmarginal,double* lambda, double *vector,int k)
{
    RMat  vMatX(mat1,numRows,numCols);
    RMat  vMatY(mat2,numRows2,numCols2);
    RMat  vMatAlphas(result,numRows,1);
    RMat  vMatOrders(orders,numRows,numCols);
   
  // if(!conjgradFP(vMatX, vMatY ,vMatAlphas,vMatOrders, sigma ) ){
  //   Rprintf("conjgradFP not converged\n");
  // }
   conjgradFP(vMatX, vMatY ,vMatAlphas,vMatOrders, sigma );
   //vMatAlphas.Print();
   //std::cout<<"after conjgradFP"  <<std::endl;
   
   double mu1 = numRows*sigma + vMatX.Sum();
   //std::cout<<"set mu1"  <<std::endl;
   // nedd to get larges eigenvalue
 //  RMat eigenVec(vector,numRows,1);
   //eigenVec = 1.0/(double)numRows;
   
  // eigenVec.Print("eigenvec");
   //eigenVec = 1.0;
   
  // std::cout<<"calc eigen"  <<std::endl;
  double mu2;
  double beta;
  //how many eigenvalues do we need ?
  if(k == 1){
    std::cout<<"only one eigen value used"  <<std::endl;
    RMat eigenVec(vector,numRows,1);
    eigenVec = 1.0/(double)numRows;
    beta = getEigen(vMatX,eigenVec,vMatOrders, sigma);
    *lambda = beta; 
    mu2 = beta *beta;
  }else{
    //mu2 =sum(lambdas$values[1]^2)
 //void Cpplanczos(RMat& Xmat, RMat& bmat,RMat& vMatOrders,int k,double*alphas,double*betas,double sigma){
      RMat bmat(numRows,1);
      bmat = 1.0/(double)numRows;
      
        // run lanczos for more iterations than eigenvalues requried, otherwise it will not necessarly 
      // give you the k largest eigenvalues.
      // reorthogorlize does the trick
     // int k_lancz = k;
    //  if(k_lancz>numRows){// too many
    //    k_lancz = k;
    //  }
    //  std::cout<<"k_lancz"<<k_lancz<<std::endl;
      RMat lalphas(k,1);
      RMat lbetas(k-1,1);
      
      Cpplanczos(vMatX, bmat,vMatOrders, k,lalphas.getValuesPtr(),lbetas.getValuesPtr(), sigma); // get tri diagonal
      
      RMat lambdaVec(lambda,k,1);
      //void getEigenSturm(double *alpha,int *nalpha,double*beta,int* k,double *result){
      getEigenSturm(lalphas.getValuesPtr(), lalphas.NumRows() ,lbetas.getValuesPtr(),k,lambdaVec.getValuesPtr() );
      beta = *lambda;
      //std::cout<<"beta:"<<beta<<std::endl;
      mu2 = lambdaVec.SquareSum();
    //  std::cout<<"mu2:"<<mu2<<std::endl;
      
      //RMat vectorVec(lambda,numRows,k);
      RMat eigenVec(vector,numRows,k);
      RMat newLambdaVec(numRows,1);
      for( unsigned int i=0;i< k ;++i ){
        //void CinverseInteration( double* A,int* nrow,int*ncol,double *lambda,double *result,double *orders,double *sigma)
        RMat resEigenVec( eigenVec.getColPtr(i+1), numRows,1);
        inverseInteration(vMatX,lambda[i], resEigenVec,vMatOrders, sigma);
        
        //update lambda
     //   std::cout<<"old lambda:"<< lambda[i]<<std::endl;
    //    fastMultiply( vMatX,resEigenVec,newLambdaVec,vMatOrders,sigma); //A*b
    //    lambda[i] = newLambdaVec.ScalarProd(resEigenVec);
    //    std::cout<<"new lambda:"<< lambda[i]<<std::endl;
        //inverseInteration( RMat& A,double lambda,RMat &result,RMat& orders,double sigma)
        
       // void inverseInteration( RMat& A,double lambda,RMat &result,RMat& orders,double sigma){
        //  CinverseInteration( double* A,int* nrow,int*ncol,double *lambda,double *result,double *orders,double *sigma)
      }
        
  }
   
  // std::cout<<"eigen calced"<<beta  <<std::endl;
//   std::cout<<beta<<std::endl;
 //  std::cout<<" "<<std::endl;
 //  eigenVec.Print();
   
// double mu2 =  beta *beta;// sum(lambdas$values[1]^2)
     
  // calculate logmarginal   
   double t_bar = (beta*mu1 - mu2) / (beta*numRows - mu1);
   double t_bar2 = t_bar*t_bar;
   double beta2 = beta*beta; ///
     
   RMat t1(1,2);
   t1(1,1) =log(beta);
   t1(1,2) =log(t_bar);
 
   RMat t2(2,2);
   t2(1,1) = t_bar2;
   t2(2,1) = -beta2;
   t2(1,2) = -t_bar;
   t2(2,2) = beta;  
   t2*=(1/(beta*t_bar2 - t_bar*beta2 ));
       
   RMat t3(2,1);
   t3(1,1)=mu1;
   t3(2,1)=mu2;
    
   //(t1 * t2).Print("t1*t2\n");
    //t3.Print("t3\n");
   //double logdetM =  (t1 * t2).ScalarProd(t3);
   
   double logdetM =  *(((t1 * t2)*t3).getValuesPtr());
  // std::cout<<"logdet:"<<logdetM<<std::endl;
   
   *logmarginal= 0.5*vMatAlphas.ScalarProd(vMatY)+ logdetM/2.0 + (numRows/2.0)*log(2.0*PI);
  // std::cout<<"   vMatY.ScalarProd( vMatAlphas ) :"<<   vMatY.ScalarProd( vMatAlphas ) <<std::endl;
   //std::cout<<"    numRows/2*log(2*PI):"<<   (numRows/2.0)*log(2.0*PI)<<std::endl;
   
  // std::cout<<"approxmarginal:"<<*logmarginal<<std::endl;


}

int compare ( const void *pa, const void *pb ){
  const double *a = (double*) pa;
  const double *b = (double*) pb;
  
  if(a[0]  > b[0] ) return 1;
  if(a[0] == b[0])  return 0;
  if(a[0]  < b[0])  return -1;
}

void sortIdx(RMat &vMatValues){

 // vMatValues.Print("before sort");
  qsort(vMatValues.getValuesPtr(), vMatValues.NumCols() ,2*sizeof(double),compare);
 // vMatValues.Print("after sort");
}

void chkBound(RMat& params,double* lower, double* upper){
  
  double *params_ptr = params.getValuesPtr();
  for(unsigned int i=0;i<params.NumRows();++i){
    if(params_ptr[i]<lower[i]){
      params_ptr[i] = lower[i];
    }else if(params_ptr[i]>upper[i]){
      params_ptr[i] = upper[i];
    } 
  }
  
}
// do just all stuff in c
/*double dummy2(RMat& X,RMat& Y,RMat& params,RMat &orders, RMat &result,RMat& X_trans){
  
  double lambda;
  double logmarginal;

  // coopy X to Xtrans.
  
  //memcpy(X_trans.getValuesPtr(),X.getValuesPtr(),sizeof(double)*X.NumRows()*X.NumCols() );
  double *Xp=X.getValuesPtr();
  double *Xtransp=X_trans.getValuesPtr();
  
  for(unsigned int i =0;i<X.NumRows()*X.NumCols();++i ){
    Xtransp[i] = Xp[i] + 100; // shift for n 
  }
  
  // transform data for me
  for(int i=1;i<=X.NumCols();++i){
    RMat point(X_trans.getColPtr(i),X_trans.NumRows(),1);
    point *=  params(i+1,1);
  }
    
//  void CgpHist(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *sigma,double *orders,double* logmarginal,double *lambda){
    //  std::cout<<"CgpHist"  <<std::endl;
  CppHist(result.getValuesPtr(),X.NumRows(),X.NumCols(),Y.NumRows(),Y.NumCols(),X_trans.getValuesPtr(),Y.getValuesPtr(),params(1,1),orders.getValuesPtr(),&logmarginal,&lambda);
  
  return(logmarginal);  
  
}*/

// this function calls the R transform function
double dummy(RMat& params,SEXP function_call,SEXP environment,double* input){
 
  //double
  //for( int i = 0;i<params.NumRows();++i){
  //  input[i] =;
  //}
  memcpy(input,params.getValuesPtr(),params.NumRows()*sizeof(double));
  
  SEXP s_fval;
  PROTECT(s_fval = eval(function_call, environment));
  if(length(s_fval) != 1) {
    error("Function should return error to minimize");
  }
  double *err = REAL(s_fval);
 // std::cout<<"res:"<< *err <<std::endl;
  double ret = *err;
  UNPROTECT(1); /* s_fval */

  return ret ;
}

//lanczos implementation, think about full reorthogornalisation......
void Cpplanczos(RMat& Xmat, RMat& bmat,RMat& vMatOrders,int k,double*alphas,double*betas,double sigma){
 
  RMat v =  bmat / sqrt(bmat.SquareSum() );
 //RMat v_prev(Xmat.NumRows(),1); // chnages htis
  RMat v_prev(Xmat.NumRows(),k); // chnages store v vectors.

  RMat w(Xmat.NumRows(),1);
  for(unsigned i = 0;i<k-1;++i){
    fastMultiply( Xmat,v, w,vMatOrders, sigma);
    alphas[i] = w.ScalarProd(v);
    w -= alphas[i]*v;
    if(i>0){
      RMat v_col_prev(v_prev.getColPtr(i),Xmat.NumRows(),1);
      w -= betas[i-1]*v_col_prev;
    }
    RMat vcol(v_prev.getColPtr(i+1),Xmat.NumRows(),1); // store v
    vcol = v;
    //v_prev = v;
    betas[i] =  sqrt(w.SquareSum());
    v= w;
    //reorthogonalize
    RMat v_prev_small(v_prev.getValuesPtr(),Xmat.NumRows(),i+1);
    //v-=  v_prev_small* (v_prev_small.Transpose()*v);
    v-=  v_prev_small* (v_prev_small.tMultiply(v));
    
    
    //std::cout<<"v_prev_small.Transpose()*v"<<std::endl;
  //  (v_prev_small.Transpose()*v).Print("TRASNPOSE\n");
   // std::cout<<"tMultiply"<<std::endl;
  //  v_prev_small.tMultiply(v).Print("tmultiply\n");
    
    std::cout<<"(v_prev_small.Transpose()*v)"<<sqrt((v_prev_small.Transpose()*v).SquareSum())<<std::endl;
   //if(sqrt((v_prev_small.Transpose()*v).SquareSum())>1e-25){
    //  std::cout<<"need to reorthogonalize again!"<<std::endl;
  //    v-=  v_prev_small* (v_prev_small.Transpose()*v);
  //  }
    v/=betas[i];
    
  }
  fastMultiply( Xmat,v, w,vMatOrders, sigma);
  alphas[k-1] = w.ScalarProd(v);
}

extern "C" {
    void CgpHist(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *sigma,double *orders,double* logmarginal,double *lambda,double* vector,int*k){
    //  std::cout<<"CgpHist"  <<std::endl;
     #if COMPILE_PARALLEL == 2
       if(threads==NULL){
         std::cout<<"CgpHist threads are not initialized. Call startThreads first!"  <<std::endl;
        
         std::cout<<"I call it for you but please call stopThreads() after you are finished!"  <<std::endl;
        
         CstartThreads();
         return;
       } 
     #endif
      CppHist(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,*sigma,orders,logmarginal,lambda, vector,*k);
    }

  void CgpHistPredict(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *mat3,double *orders){
  //  std::cout<<"CgpHistPredict"  <<std::endl;
    CppHistPredict(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,mat3,orders);
  }
  
// predict variance. coarse and fine are supported
  void CgpHistVarianceCoarse(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double*lambda,double* sigma, double *orders){
    CppHistVarianceCoarse(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,*lambda,*sigma,orders);
  }
  void CgpHistVarianceFine(double* result,int *numRows,int *numCols,int *numRows2,int *numCols2,double* X,double* pred,double *lambda,int *nlambda, double *vectors,double *sigma,double* orders){
    CppHistVarianceFine(result,*numRows,*numCols,*numRows2,*numCols2,X,pred,lambda, *nlambda,vectors,*sigma,orders);
  }
    

  // for debugging, set up matrices and call lanczos....
  void Clanczos(double *X, int *nrow,int *ncol, double *b, int* k,double*alphas,double*betas,double* orders,double *sigma){
    RMat Xmat(X,*nrow,*ncol);
    RMat bmat(b,*nrow,1);
    RMat vMatOrders(orders,*nrow,*ncol);
    Cpplanczos(Xmat, bmat,vMatOrders, *k,alphas,betas,*sigma);
  }
  
  //wrapper for debuggin. call eigen sturm
  void CgetEigenSturm(double *alpha,int *nalpha,double*beta,int* k,double *result){
    getEigenSturm(alpha,*nalpha,beta,*k,result);
  }
  
  //wrapper for debugging.. call inverse iteration
  void CinverseInteration( double* A,int* nrow,int*ncol,double *lambda,double *result,double *orders,double *sigma){
    RMat vMatA(A,*nrow,*ncol);
    RMat vMatorders(orders,*nrow,*ncol);
    RMat vMatResult(result,*nrow,1);
    
    inverseInteration( vMatA,*lambda,vMatResult, vMatorders,*sigma);
  }
 
 
  // downhill simplex with boundery check for estimation of hyperparameters
  void Cdownhillsimplex(SEXP func,SEXP sys,int* nParams,double* bp,double* lower,double* upper,int *resPtr, double *alpha,double *gamma,double *beta,double *sigma,unsigned int *it,double *tol){
    SEXP arg;
    SEXP environment;
    
    PROTECT(environment = sys);
    PROTECT(arg = allocVector(REALSXP, *nParams) );
 // memcpy(REAL(arg), 5, (*nParams) * sizeof(double));
    double *input = REAL(arg);
 
    //for( int i = 0;i<*nParams;++i){
    //  input[i] =i;
  //  }
    
    SEXP function_call;
    PROTECT(function_call = lang2(func, arg));
    
    //SEXP s_fval;
   // PROTECT(s_fval = eval(function_call, environment));
    
   /* if(length(s_fval) != 1) {
      error("Function should return error to minimize");
    }
    double *err = REAL(s_fval);
    std::cout<<"res:"<< *err <<std::endl;
    
    UNPROTECT(4); /* s_fval */
        
        
        /////////
        
        int Np1 = (*nParams)+1;
        RMat  vMatbp(bp,*nParams,Np1); //parameters but make downwords
        RMat  vMatValues( 2,Np1); //functiosn values..next row idx
        
        RMat  m(*nParams ,1); //mean over params
        RMat  r(*nParams ,1); //new parameters
        RMat  e(*nParams ,1); //new parameters
        
//  RMat  vMatIdx( (*nParams)+1,1); //functiosn values..
        
        double fcRes;
        for(int i=1;i<=Np1;++i){
          RMat point(vMatbp.getColPtr(i),*nParams ,1); // get vector of params
          fcRes = dummy(point,function_call,environment, input);
          vMatValues(1,i) = fcRes;
          //res =c(res, fc(bp[k,]) )
        }
        
        for(unsigned i=1;i<=vMatValues.NumCols();++i ){ // set idx...
          vMatValues(2,i)=i;
        }
        
        
        //vMatValues.Print("values");
       // std::cout<<"iterations:"<<*it<<std::endl;
        
        for(unsigned int i=0;i<*it;++i){
          // sort my function valiues
          
          sortIdx(vMatValues);
          // check convergenc
        // std::cout<<"convergence value"<< fabs(vMatValues(1,1)-vMatValues(1,Np1)) <<std::endl;
          
          
          if(fabs( vMatValues(1,1) -vMatValues(1,Np1) )  < *tol){
          //  cat(c('converged after:',i,' iterations. value:', abs(bp[1,N+1]-bp[N+1,N+1]) ,'\n' ) )
            std::cout<<"converged after:"<<i<<" iterations."<< fabs(vMatValues(1,1)-vMatValues(1,Np1))<<std::endl;
            break;
          }
          
        //  vMatbp.Print("vMat");
          //have to do stuff because not the last solution...
          int idx ;
          m.SetZero();
          for(int x=1;x<vMatbp.NumCols();++x){
            idx =  (int) vMatValues(2,x);
            RMat setup(vMatbp.getColPtr(idx), (*nParams) ,1);
            m += setup;
          }
          //m = vMatbp.RowSums();
          m /= (*nParams); // mean
          
      //    m.Print("\n\n\nm");
          
          idx =(int) vMatValues(2,Np1);
       //   m.Print("mean");
         // std::cout<<"calculated mean"<<std::endl;  
            
          RMat point(vMatbp.getColPtr(idx),*nParams ,1); // get vector of params
       //   point.Print("point...");
          r =(1.+(*alpha)) *m - point*(*alpha) ;//bp[N+1,1:N]
       //   r.Print("r");
          chkBound(r, lower, upper);
          
      //    std::cout<<"after chk bound r"<<std::endl;  
      //    r.Print("R\n");
          //np = fc(r) //querry point
          fcRes =dummy(r,function_call,environment, input);;
            
           if(fcRes<vMatValues(1,1)){ //##better
              e =(1+(*gamma)) *m -point* (*gamma);
          //   std::cout<<"after calc e"<<std::endl;  
         //    e.Print("e\n");
              chkBound(e,lower,upper );
               
               double fcRes2 = dummy(e,function_call,environment, input);;
               
               if(fcRes2<fcRes){ //e better
                 point = e;
                 vMatValues(1,Np1) = fcRes2 ;
               }else{
                 point = r;
                 vMatValues(1,Np1) = fcRes ;
               }
        //       vMatValues.Print("vMatValues\n");
         //      vMatbp.Print("vMatbp\n");
               
               continue;
           }
           
           
           if(fcRes<vMatValues(1,*nParams)){
            // idx =(int) vMatValues(2,*nParams);// set second worst solution
           //  point.setValuesPtr(vMatbp.getColPtr(idx));
             
             point = r;
             vMatValues(1,Np1) = fcRes ;
             
       //      vMatValues.Print("set r vMatValues\n");
       //      vMatbp.Print("vMatbp\n");
             
             continue;
           }
    //       std::cout<<"fc res"<< fcRes<<std::endl;
     //      std::cout<<"vMatValues(1,Np1) res"<< vMatValues(1,Np1)<<std::endl;
           
           
           if(vMatValues(1,Np1)<fcRes){ 
             e = (*beta) *m+(1.- *beta)*point;
      //       std::cout<<"n+1 is better"<<std::endl;
           }else{
             e = (*beta) *m+(1.- *beta)*r;
      //       std::cout<<"fc is better"<<std::endl;
        //     r.Print("my r");
           }
           
        //   e.Print("cp before ckbound");
           
           chkBound(e,lower,upper );
           fcRes = dummy(e,function_call,environment, input);
        //  e.Print("get cp query");
           if(fcRes< vMatValues(1,Np1)){
               point = e;
               vMatValues(1,Np1)= fcRes ;
           //    vMatValues.Print("vMatValues accteo cp");
           //    vMatbp.Print("vMatbp");
               
               continue;
           }
           
           // idx =(int) vMatValues(2,*nParams);// set second worst solution
           //  point.setValuesPtr(vMatbp.getColPtr(idx));
           
           idx = (int) vMatValues(2,1);
           RMat best(vMatbp.getColPtr(idx),*nParams ,1); 
           
           for(int k=2;k<= vMatbp.NumCols();++k){ 
             idx = (int) vMatValues(2,k);
             point.setValuesPtr(vMatbp.getColPtr(idx));
             
             e = best+ (*sigma) *( point -best);
             chkBound(e,lower,upper );
             point = e;
             fcRes = dummy(e,function_call,environment, input);
             vMatValues(1,k)= fcRes;
           }
           
        //   vMatValues.Print("vMatValues last");
       //    vMatbp.Print("vMatbp");
           

        }
        *resPtr =vMatValues(2,1) ;
        
    UNPROTECT(3);
   
    //Cdownhillsimplex SEXP
   // std::cout<<"OK"<<std::endl;
  }

  
  /*
  ////////// downhill all in c will be deleted....
void CdownhillsimplexIN(int* nParams,double* bp,double* lower,double* upper,int *resPtr, double *alpha,double *gamma,double *beta,double *sigma,unsigned int *it,double *tol, double *X, double *Y,double *orders,double *nrow, double *ncol){


  //double dummy2(RMat& X,RMat& Y,RMat& params,RMat &orders, RMat &result,RMat& X_trans){
  RMat Xmat(X,*nrow,*ncol);
  RMat Xtransmat(*nrow,*ncol);
  RMat Ymat(Y,*nrow,1);
  RMat result(*nrow,1);
  RMat ordersmat(orders,*nrow,*ncol);
  
  
  /////////
  
  int Np1 = (*nParams)+1;
  RMat  vMatbp(bp,*nParams,Np1); //parameters but make downwords
  RMat  vMatValues( 2,Np1); //functiosn values..next row idx
  
  RMat  m(*nParams ,1); //mean over params
  RMat  r(*nParams ,1); //new parameters
  RMat  e(*nParams ,1); //new parameters
  
  //  RMat  vMatIdx( (*nParams)+1,1); //functiosn values..
  
  double fcRes;
  for(int i=1;i<=Np1;++i){
    RMat point(vMatbp.getColPtr(i),*nParams ,1); // get vector of params
    fcRes = dummy2(Xmat,Ymat,point,ordersmat, result,Xtransmat);
    vMatValues(1,i) = fcRes;
    //res =c(res, fc(bp[k,]) )
  }
  
  for(unsigned i=1;i<=vMatValues.NumCols();++i ){ // set idx...
    vMatValues(2,i)=i;
  }
  
  
  //vMatValues.Print("values");
  // std::cout<<"iterations:"<<*it<<std::endl;
  
  for(unsigned int i=0;i<*it;++i){
    // sort my function valiues
    
    sortIdx(vMatValues);
    // check convergenc
    // std::cout<<"convergence value"<< fabs(vMatValues(1,1)-vMatValues(1,Np1)) <<std::endl;
    
    
    if(fabs( vMatValues(1,1) -vMatValues(1,Np1) )  < *tol){
      //  cat(c('converged after:',i,' iterations. value:', abs(bp[1,N+1]-bp[N+1,N+1]) ,'\n' ) )
      std::cout<<"converged after:"<<i<<" iterations."<< fabs(vMatValues(1,1)-vMatValues(1,Np1))<<std::endl;
      break;
    }
    
    //  vMatbp.Print("vMat");
    //have to do stuff because not the last solution...
    int idx ;
    m.SetZero();
    for(int x=1;x<vMatbp.NumCols();++x){
      idx =  (int) vMatValues(2,x);
      RMat setup(vMatbp.getColPtr(idx), (*nParams) ,1);
      m += setup;
    }
    //m = vMatbp.RowSums();
    m /= (*nParams); // mean
    
    //    m.Print("\n\n\nm");
    
    idx =(int) vMatValues(2,Np1);
    //   m.Print("mean");
    // std::cout<<"calculated mean"<<std::endl;  
    
    RMat point(vMatbp.getColPtr(idx),*nParams ,1); // get vector of params
    //   point.Print("point...");
    r =(1.+(*alpha)) *m - point*(*alpha) ;//bp[N+1,1:N]
    //   r.Print("r");
    chkBound(r, lower, upper);
    
    //    std::cout<<"after chk bound r"<<std::endl;  
    //    r.Print("R\n");
    //np = fc(r) //querry point
    fcRes =dummy2(Xmat,Ymat,r,ordersmat, result,Xtransmat);
    
    if(fcRes<vMatValues(1,1)){ //##better
      e =(1+(*gamma)) *m -point* (*gamma);
      //   std::cout<<"after calc e"<<std::endl;  
      //    e.Print("e\n");
      chkBound(e,lower,upper );
      
      double fcRes2 = dummy2(Xmat,Ymat,e,ordersmat, result,Xtransmat);
      
      if(fcRes2<fcRes){ //e better
        point = e;
        vMatValues(1,Np1) = fcRes2 ;
      }else{
        point = r;
        vMatValues(1,Np1) = fcRes ;
      }
      //       vMatValues.Print("vMatValues\n");
      //      vMatbp.Print("vMatbp\n");
      
      continue;
    }
    
    
    if(fcRes<vMatValues(1,*nParams)){
      // idx =(int) vMatValues(2,*nParams);// set second worst solution
      //  point.setValuesPtr(vMatbp.getColPtr(idx));
      
      point = r;
      vMatValues(1,Np1) = fcRes ;
      
      //      vMatValues.Print("set r vMatValues\n");
      //      vMatbp.Print("vMatbp\n");
      
      continue;
    }
    //       std::cout<<"fc res"<< fcRes<<std::endl;
    //      std::cout<<"vMatValues(1,Np1) res"<< vMatValues(1,Np1)<<std::endl;
    
    
    if(vMatValues(1,Np1)<fcRes){ 
      e = (*beta) *m+(1.- *beta)*point;
      //       std::cout<<"n+1 is better"<<std::endl;
    }else{
      e = (*beta) *m+(1.- *beta)*r;
      //       std::cout<<"fc is better"<<std::endl;
      //     r.Print("my r");
    }
    
    //   e.Print("cp before ckbound");
    
    chkBound(e,lower,upper );
    fcRes = dummy2(Xmat,Ymat,e,ordersmat, result,Xtransmat);
    //  e.Print("get cp query");
    if(fcRes< vMatValues(1,Np1)){
      point = e;
      vMatValues(1,Np1)= fcRes ;
      //    vMatValues.Print("vMatValues accteo cp");
      //    vMatbp.Print("vMatbp");
      
      continue;
    }
    
    // idx =(int) vMatValues(2,*nParams);// set second worst solution
    //  point.setValuesPtr(vMatbp.getColPtr(idx));
    
    idx = (int) vMatValues(2,1);
    RMat best(vMatbp.getColPtr(idx),*nParams ,1); 
    
    for(int k=2;k<= vMatbp.NumCols();++k){ 
      idx = (int) vMatValues(2,k);
      point.setValuesPtr(vMatbp.getColPtr(idx));
      
      e = best+ (*sigma) *( point -best);
      chkBound(e,lower,upper );
      point = e;
      fcRes = dummy2(Xmat,Ymat,e,ordersmat, result,Xtransmat);
      vMatValues(1,k)= fcRes;
    }
    
    
  }
  *resPtr =vMatValues(2,1) ;
  
}
  */
  ///////////// will be deleted
  
#if COMPILE_PARALLEL == 2
  void startThreads(){
    CstartThreads();
  }
  void stopThreads(){
    CstopThreads();
  }
      
#endif
}
