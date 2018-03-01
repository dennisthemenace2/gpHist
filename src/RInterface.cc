#include <R.h>     // R functions
#include <Rmath.h> // Rmath

#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "RMat.h"


#include <iostream> // For cout, cerr, endl, and flush
#include <algorithm>
#include <stdlib.h>
#include <float.h>

#include <pthread.h>



void Cpplanczos(RMat& Xmat, RMat& bmat,RMat& vMatOrders,unsigned int k,double*alphas,double*betas,double sigma);

// calculate sturm seq
void sturmSeq(double* alpha,unsigned int nalpha, double *beta,double lambda,double* result){
  
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
void getEigenSturm(double *alpha,unsigned int nalpha,double*beta,unsigned int k,double *result){
  
  if(k==0 || k>nalpha){//this should not happen anyway
    //std::cout<<"invalid number of eigenvalues. k="<<k<<" nalpha="<< nalpha<< "\n";
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
  
  RMat seqhigh(nalpha,1);
  sturmSeq(alpha,nalpha,beta,bmax,seqhigh.getValuesPtr()  );

  RMat seqlow(nalpha,1);
  sturmSeq(alpha,nalpha,beta,bmin,seqlow.getValuesPtr()  );

  int changeslow_first = calcSeqChanges(seqlow);
  int changeslow= changeslow_first;
  
  int changeshigh = calcSeqChanges(seqhigh);
 // int dif = abs(changeshigh -changeslow ); //monkas not used anymore

  for(unsigned int e=0; e<k; ++e ){
    double lb = bmin;
    double ub = bmax;

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
      
      if(fabs(lb-ub)<0.00001){
        break;
      }
    }
    result[e]= (lb+ub)/2.0;
    changeshigh = changeslow;
    changeslow = changeslow_first;
    
  }
  
}


///// Choose of parallel or not
#define COMPILE_PARALLEL 0 // set to 0 for no parallel execution

#if COMPILE_PARALLEL == 0

void fastMultiply( RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){

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

    for(unsigned int k=0;k<NR-1;++k ){ // estimate prev
      idx =  *(idx_prt++); //orders(k,i);
      idx -=1;
      
      result_ptr[idx] +=  matcount + X_ptr[idx]*YsumMat_ptr[k]; // Ysum
      matcount = matcount + X_ptr[idx] * Y_ptr[idx];
      
    }
    idx = *(idx_prt++);//orders(NR,i);
    idx -=1; 
    result_ptr[idx] += matcount + X_ptr[idx] * YsumMat_ptr[NR-1];
    
  }
  result += sigma*Y ;

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

    result[idx] += matcount + X_ptr[idx] * YsumMat_ptr[NR-1];
    
  }
  
  pthread_exit( (void *)result);
}


void fastMultiply( RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){
  result.SetZero();
  unsigned int NR = X.NumRows();
  
  unsigned int numThreads = X.NumCols(); // create a thread for each dimension
  if(numThreads>NUM_THREADS){
    numThreads = NUM_THREADS;
  }
  unsigned int each = X.NumCols()/numThreads;
  
  //create threads  
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
    double *ret;
    pthread_join(threads[t], (void**) &ret);
    
    RMat resMat(ret,NR ,1);
    result +=resMat; // add dimensions
    
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
  

  pthread_mutex_lock (&data->shared->initmutex);
  //   cout<<"get init lock"<<endl;
  pthread_mutex_lock (&data->pmutex);       
  data->shared->initcnt +=1; 
  pthread_mutex_unlock (&data->shared->initmutex);

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

      for(unsigned int k=0;k<NR-1;++k ){ // estimate prev
        idx =  *(idx_prt++); //orders(k,i);
        idx -=1;
        
        result[idx] +=  matcount + X_ptr[idx]*YsumMat_ptr[k]; // Ysum
        matcount = matcount + X_ptr[idx] * Y_ptr[idx];
        
      }
      idx = *(idx_prt++);//orders(NR,i);
      idx -=1; 
      
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

  pthread_exit(NULL);
}

thread_data_Mutex *tdata=0;  
pthread_t *threads=0;
pthread_attr_t attr;

#define NUM_THREADS 4


void CstartThreads(){
  if(threads != NULL){
   // std::cout<<"Already initialized!"<<std::endl;
    Rprintf("Already initialized!\n");
    return;
  }
  tdata = new thread_data_Mutex[NUM_THREADS];
  
  threads= new  pthread_t[NUM_THREADS];
  thread_data_shred *shared = new thread_data_shred;
  
  
  /* Initialize and set thread detached attribute */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  pthread_mutex_init(&shared->initmutex, NULL);

  shared->initcnt = 0;
  shared->run = true;

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

  //std::cout<<"init threads"<<std::endl;
}
void CstopThreads(){
  if(threads == NULL){
  //  std::cout<<"no threads initilized!"<<std::endl;
    Rprintf("no threads initilized!\n");
    return;
  }
  
  thread_data_shred *shared = tdata[0].shared;
  
  shared->run = false;
  
  // make them run
  for(unsigned int c=0;c<NUM_THREADS;c++){
    pthread_mutex_unlock(&tdata[c].tmutex);
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


bool conjgradFP( RMat& A, RMat& b,RMat& x ,RMat& orders,double sigma, unsigned int max_iterations=2000 ) {
  
  RMat res( x.NumRows(),1 );
  
  fastMultiply(A,x,res,orders,sigma); //A%*%x;
  

  RMat r= b-res;// 
  RMat p=r;
  
  double rsold=r.SquareSum();//t(r)%*%r 


  double rsnew ;
  
  for (unsigned int i=1;i<=max_iterations;++i){
    fastMultiply(A,p,res,orders,sigma);// res =  A%*%p;


    double alpha= rsold / p.ScalarProd(res) ; // the scalar product can overflow :(
    
    x += alpha*p;
    
    if(i % 50 == 0){
      fastMultiply(A,x,res,orders,sigma);
      r=b-res; //#A%*%x;
    }else{
      r -= alpha*res;
    }
    
    
    rsnew= r.SquareSum();//t(r)%*%r;

    if(sqrt(rsnew) <1e-10){
       // std::cout<<"conjgradFP converged after "<<i<<std::endl;
      return true;
    }
    p*=(rsnew/rsold);
    p+=r;
    rsold=rsnew;
  }

  return false;
}
// power method to get larges eigenvalue and vector
double getEigen( RMat& A,RMat &b,RMat& orders,double sigma){
  double last = 0;
  RMat res( b.NumRows(),1 );
  
  for(unsigned int i=0;i<2000;++i){ // or give it a max iterations parameter
    fastMultiply(A,b,res,orders,sigma); //A%*%b;
    double norm = sqrt(res.SquareSum()) ;
    b = res;
    b /= norm;
    if(fabs(last-norm)<1e-10 ){
	//std::cout<<"converged poweriteration:"<< i<<std::endl;
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
  RMat tmp( N,1);

//see
//http://www.netlib.org/utk/people/JackDongarra/etemplates/node96.html

  for(unsigned int i=0;i<2000;++i){

    conjgradFP( A, result, tmp , orders,sigma-lambda, 1000 );
   
    double th = result.ScalarProd(tmp);
    double len = sqrt (  (tmp - (result*th) ).SquareSum() );

    if(len< 1e-10*fabs(th)){
   //   std::cout << "converged "<<len<< " 1e-10*fabs(th)"<< 1e-10*fabs(th)<<std::endl<<std::flush; 
      result = tmp;
      result /= th;
      break;
    }

    double norm = sqrt(tmp.SquareSum());
    result = tmp;
    result /= norm;

  }
  
}

// two helper functions for predictions
double sumOrdered(RMat& alpha,RMat&orders,unsigned int start,unsigned int end,unsigned int col){
  double result = 0;
  double *orders_ptr = orders.getColPtr(col)+(start-1);
  double *alpha_ptr = alpha.getValuesPtr();
  
  for(unsigned int i= start;i<=end;++i){
    int idx = *(orders_ptr++);
    idx -=1; // or set pointer to invalid location and skip this ?
    result += alpha_ptr[idx];
  }
  return result;
}

double sumSquareOrdered(RMat& X,RMat&orders,unsigned int start,unsigned int end,unsigned int col){
  double result = 0;
  double *X_ptr = X.getColPtr(col);
  double *orders_ptr = orders.getColPtr(col)+(start-1);
  
  for(unsigned int i= start;i<=end;++i){
    unsigned int idx = *(orders_ptr++);
    idx -=1;
    double val = X_ptr[idx];
    result += val*val;
  }
  
  return result;
}

double sumOrdered(RMat& X,RMat& alpha,RMat&orders,unsigned int start,unsigned int end,unsigned int col){
  double result = 0;
  double *alpha_ptr = alpha.getValuesPtr();
  double *X_ptr = X.getColPtr(col);
  double *orders_ptr = orders.getColPtr(col)+(start-1);
  
  for(unsigned int i= start;i<=end;++i){
    unsigned int idx = *(orders_ptr++);
    idx -=1;
    result += X_ptr[idx]*alpha_ptr[idx];
  }
  
  return result;
}

unsigned int findIdx(RMat& X,RMat& orders,double value,unsigned int col){
  // try intersection...
  double *orders_ptr = orders.getColPtr(col); // pointer to column
  double *X_ptr = X.getColPtr(col); // pointer to column
  
  unsigned int min = 0;
  unsigned int max = X.NumRows()-1;
  unsigned int mid; 
  do{
    mid = (min+max)/2;
    if(X_ptr[(unsigned int) orders_ptr[mid]-1 ] >value  ){ // go left
      max = mid;
    }else{ // go right should i consider equal 
      min = mid;
    }
    
  }while(max-min >1);
  
  if( value <= X_ptr[(unsigned int) orders_ptr[min]-1 ] ){
    mid = min;//-1;
  }else{
    mid = max;
    if(value>X_ptr[(unsigned int) orders_ptr[max]-1 ]){
      mid = max+1;
    }
  }
  
  return mid;
}

// coarse approximation
void CppHistVarianceCoarse(double* result,unsigned int numRows,unsigned int numCols,unsigned int numRows2,unsigned int numCols2,double* mat1,double* mat2,double lambda,double sigma,double* orders){
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatPred(mat2,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
  RMat  vMatOrders(orders,numRows,numCols);
  
  for(unsigned int i=1; i<=numRows2;++i ){
    double A=0;
    double B=0;
    
    for(unsigned int d=1; d<= numCols2;++d){
      int idx = findIdx(vMatX,vMatOrders, vMatPred(i,d),d);

      if(idx>0){
        A = sumSquareOrdered(vMatX,vMatOrders,1,idx,d);
      }else{
        A=0;
        idx=0;
      }
      
      B = (numRows-idx )* ( vMatPred(i,d)*vMatPred(i,d) );
    }
    vMatResult(i,1) = A+B;
  }
  vMatResult *= (1.0/lambda);
  vMatResult = vMatPred.RowSums() -vMatResult;
}

// finer approximation use all calculated lambdas
void CppHistVarianceFine(double* result,unsigned int numRows,unsigned int numCols,unsigned int numRows2,unsigned int numCols2,double* X,double* pred,double *lambda,unsigned int nlambda, double *vectors,double sigma,double* orders){
  RMat  vMatX(X,numRows,numCols);
  RMat  vMatPred(pred,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
  RMat  vMatVectors(vectors,numRows,nlambda);
  
  RMat  k_star(numRows,1); // number of alphas 
  
  vMatResult =   vMatPred.RowSums(); // init this one...
  for(unsigned int s=1; s<=numRows2;++s ){ // each sample to predict
    //calc k_star
    k_star.SetZero();
    //estimate k_star
    for(unsigned int x=1;x<=numCols;++x){ // each dimension
      double *col = vMatX.getColPtr(x);
      double *starpointer =k_star.getValuesPtr();
      double value =  vMatPred(s,x); // check against this value
      
      for(unsigned int y=1;y<=numRows;++y){ // each row in X
        if(value< *col){ // value smaller X
          *starpointer +=value;
        }else{
          *starpointer +=*col;
        }
        ++col; // next xvalue
        ++starpointer; // next k_star value
      }
    }
    
    double secondTerm =0;
    double sumOfProjectionLengths = 0;
    
    for(unsigned int i=0; i<nlambda;++i ) {
      RMat  vector(vMatVectors.getColPtr(i+1), numRows,1); // get eigenvector

      double projectionLength =  vector.ScalarProd( k_star) ; //
      double projectionLength2= projectionLength*projectionLength;
      secondTerm += ( 1.0 / lambda[i] ) *projectionLength2;

      sumOfProjectionLengths +=  projectionLength2;
    }
    double k_dist = k_star.SquareSum(); 

    secondTerm += ( 1.0 / lambda[nlambda-1] ) * ( k_dist - sumOfProjectionLengths );
    
    if(k_dist - sumOfProjectionLengths <0){
      Rprintf("sumOfProjectionLengths doesnt add up, propably ghost eigenvalues\n");
    }

    vMatResult(s,1) -= secondTerm;

  }// end for each sample
  
}


void CppHistPredict(double *result,
                    unsigned int numRows, unsigned int numCols,unsigned int numRows2 ,unsigned int numCols2,
                    double *mat1, double *mat2,double* mat3,double *orders){
  
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatPred(mat2,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
  RMat  vMatAlpha(mat3,numRows,1);
  RMat  vMatOrders(orders,numRows,numCols);
  
  //  pred = c()
  for(unsigned int i=1; i<=numRows2;++i ){
    double C=0;
    double A=0;
    double B=0;
    
    for(unsigned int d=1; d<= numCols2;++d){
      // find idx....
      
      unsigned int idx = findIdx(vMatX,vMatOrders, vMatPred(i,d),d);
      
      if(idx>0){
        if(idx==numRows){//###largest
          C = C+ sumOrdered(vMatX,vMatAlpha,vMatOrders,1,numRows,d);
          continue;
        }else{
          A = sumOrdered(vMatX,vMatAlpha,vMatOrders,1,idx,d);
        }
      }else{ //##point even further left....
        A=0;
      }
      B = sumOrdered(vMatAlpha,vMatOrders,idx+1,numRows,d) * vMatPred(i,d);
      C = C+A+B;
    }
    
    vMatResult(i,1)=C;
  }
  
}

void CppHist(double *result,
             unsigned int numRows, unsigned int numCols,unsigned int numRows2 ,unsigned int numCols2,
             double *mat1, double *mat2,double sigma , double *orders, double* logmarginal,double* lambda, double *vector,unsigned int k)
{
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatY(mat2,numRows2,numCols2);
  RMat  vMatAlphas(result,numRows,1);
  RMat  vMatOrders(orders,numRows,numCols);
  

  conjgradFP(vMatX, vMatY ,vMatAlphas,vMatOrders, sigma );

    
  
  double mu1 = numRows*sigma + vMatX.Sum();
  
  double mu2;
  double beta;
  //how many eigenvalues do we need ?
  if(k == 1){
    RMat eigenVec(vector,numRows,1);
    eigenVec = 1.0/(double)numRows;
    beta = getEigen(vMatX,eigenVec,vMatOrders, sigma);
    *lambda = beta; 
    mu2 = beta *beta;
  }else{
    RMat bmat(numRows,1);
    bmat = 1.0/(double)numRows;
    
    // run lanczos for more iterations than eigenvalues requried, otherwise it will not necessarly 
    // give you the k largest eigenvalues.
    unsigned int k_lancz = ceil(k*1.5);// according to wiki 1.5
    if(k_lancz>numRows){// too many
      k_lancz = k;
    }
    RMat lalphas(k_lancz,1);
    RMat lbetas(k_lancz-1,1);
    

    Cpplanczos(vMatX, bmat,vMatOrders, k_lancz,lalphas.getValuesPtr(),lbetas.getValuesPtr(), sigma); // get tri diagonal
  
  
    RMat lambdaVec(lambda,k,1);


    getEigenSturm(lalphas.getValuesPtr(), lalphas.NumRows() ,lbetas.getValuesPtr(),k,lambdaVec.getValuesPtr() );

    beta = *lambda;

    mu2 = lambdaVec.SquareSum();


    RMat eigenVec(vector,numRows,k);
  //  RMat newLambdaVec(numRows,1);
    for( unsigned int i=0;i< k ;++i ){
      RMat resEigenVec( eigenVec.getColPtr(i+1), numRows,1);
      inverseInteration(vMatX,lambda[i], resEigenVec,vMatOrders, sigma);
    }
    
  }

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
  t2*=(1.0/(beta*t_bar2 - t_bar*beta2 ));
  
  RMat t3(2,1);
  t3(1,1)=mu1;
  t3(2,1)=mu2;
   
  double logdetM =  *(((t1 * t2)*t3).getValuesPtr());

  *logmarginal= 0.5*vMatAlphas.ScalarProd(vMatY)+ logdetM/2.0 + (numRows/2.0)*log(2.0*PI);
}

int compare ( const void *pa, const void *pb ){
  const double *a = (double*) pa;
  const double *b = (double*) pb;
  
  int ret = 0;	
  if(a[0]  > b[0] ){
	ret = 1;
  }else if(a[0]  < b[0]){
	 ret =  -1;
  }

  return ret;
}

void sortIdx(RMat &vMatValues){
  qsort(vMatValues.getValuesPtr(), vMatValues.NumCols() ,2*sizeof(double),compare);
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

// this function calls the R transform function
double dummy(RMat& params,SEXP function_call,SEXP environment,double* input){
  
  memcpy(input,params.getValuesPtr(),params.NumRows()*sizeof(double));
  
  SEXP s_fval;
  PROTECT(s_fval = eval(function_call, environment));
  if(Rf_length(s_fval) != 1) {
    error("Function should return error to minimize");
  }
  double *err = REAL(s_fval);

  double ret = *err;
  UNPROTECT(1); /* s_fval */

  return ret ;
}

//lanczos implementation, think about full reorthogornalisation......
void Cpplanczos(RMat& Xmat, RMat& bmat,RMat& vMatOrders,unsigned int k,double*alphas,double*betas,double sigma){
  
  RMat v =  bmat / sqrt(bmat.SquareSum() );
  //RMat v_prev(Xmat.NumRows(),1); // chnages htis
  RMat v_prev(Xmat.NumRows(),k+1); // chnages store v vectors.

  RMat vcol(v_prev.getColPtr(1),Xmat.NumRows(),1); // store v
  vcol = v;
  
  RMat w(Xmat.NumRows(),1);
  for(unsigned int i = 0;i<k-1;++i){
    fastMultiply( Xmat,v, w,vMatOrders, sigma);

    alphas[i] = w.ScalarProd(v); // calc allha
    w -= alphas[i]*v; // substract alpha v
    
    if(i>0){ // subtract last b*v
      RMat v_col_prev(v_prev.getColPtr(i),Xmat.NumRows(),1);
      w -= betas[i-1]*v_col_prev;
    }
 
    RMat vcol(v_prev.getColPtr(i+1),Xmat.NumRows(),1); // store v
    vcol = v;
    v= w;
    
    //reorthogonalize
    RMat v_prev_small(v_prev.getValuesPtr(),Xmat.NumRows(),i+2);
    //v-=  v_prev_small* (v_prev_small.Transpose()*v);
    RMat h1=v_prev_small.tMultiply(v);
    v-=  v_prev_small* h1;
    
   if(sqrt(h1.SquareSum())>1.490116e-07){
      v-=  v_prev_small* v_prev_small.tMultiply(v);
    }
    betas[i] =  sqrt(w.SquareSum());
    
    v/=betas[i];
    
  }
  fastMultiply( Xmat,v, w,vMatOrders, sigma);
  alphas[k-1] = w.ScalarProd(v);
}

extern "C" {
  void CgpHist(double *result,double *mat1,unsigned int *numRows,unsigned int *numCols, double *mat2,unsigned int *numRows2,unsigned int *numCols2,double *sigma,double *orders,double* logmarginal,double *lambda,double* vector,unsigned int *k){
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
  
  void CgpHistPredict(double *result,double *mat1,unsigned int *numRows,unsigned int *numCols, double *mat2,unsigned int *numRows2,unsigned int *numCols2,double *mat3,double *orders){
    CppHistPredict(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,mat3,orders);
  }
  
  // predict variance. coarse and fine are supported
  void CgpHistVarianceCoarse(double *result,double *mat1,unsigned int *numRows,unsigned int *numCols, double *mat2,unsigned int *numRows2,unsigned int *numCols2,double*lambda,double* sigma, double *orders){
    CppHistVarianceCoarse(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,*lambda,*sigma,orders);
  }
  void CgpHistVarianceFine(double* result,unsigned int *numRows,unsigned int *numCols,unsigned int *numRows2,unsigned int *numCols2,double* X,double* pred,double *lambda,unsigned int *nlambda, double *vectors,double *sigma,double* orders){
    CppHistVarianceFine(result,*numRows,*numCols,*numRows2,*numCols2,X,pred,lambda, *nlambda,vectors,*sigma,orders);
  }

  // downhill simplex with boundery check for estimation of hyperparameters
  void Cdownhillsimplex(SEXP func,SEXP sys,unsigned int* nParams,double* bp,double* lower,double* upper,unsigned int *resPtr, double *alpha,double *gamma,double *beta,double *sigma,unsigned int *it,double *tol){
    SEXP arg;
    SEXP environment;
    
    PROTECT(environment = sys);
    PROTECT(arg = allocVector(REALSXP, *nParams) );

    double *input = REAL(arg);

    SEXP function_call;
    PROTECT(function_call = lang2(func, arg));
    
    unsigned int Np1 = (*nParams)+1;
    RMat  vMatbp(bp,*nParams,Np1); //parameters but make downwords
    RMat  vMatValues( 2,Np1); //functiosn values..next row idx
    
    RMat  m(*nParams ,1); //mean over params
    RMat  r(*nParams ,1); //new parameters
    RMat  e(*nParams ,1); //new parameters

    double fcRes;
    for(unsigned int i=1;i<=Np1;++i){
      RMat point(vMatbp.getColPtr(i),*nParams ,1); // get vector of params
      fcRes = dummy(point,function_call,environment, input);
      vMatValues(1,i) = fcRes;
    }
    
    for(unsigned i=1;i<=vMatValues.NumCols();++i ){ // set idx...
      vMatValues(2,i)=i;
    }

    for(unsigned int i=0;i<*it;++i){
      // sort my function valiues
      
      sortIdx(vMatValues);
      // check convergenc

      if(fabs( vMatValues(1,1) -vMatValues(1,Np1) )  < *tol){
       // std::cout<<"converged after:"<<i<<" iterations."<< fabs(vMatValues(1,1)-vMatValues(1,Np1))<<std::endl;
        break;
      }
      
      //have to do stuff because not the last solution...
      int idx ;
      m.SetZero();
      for(unsigned int x=1;x<vMatbp.NumCols();++x){
        idx =  (int) vMatValues(2,x);
        RMat setup(vMatbp.getColPtr(idx), (*nParams) ,1);
        m += setup;
      }

      m /= (*nParams); // mean
      
      idx =(unsigned int) vMatValues(2,Np1);
      // std::cout<<"calculated mean"<<std::endl;  
      
      RMat point(vMatbp.getColPtr(idx),*nParams ,1); // get vector of params
      r =(1.+(*alpha)) *m - point*(*alpha) ;//bp[N+1,1:N]
      chkBound(r, lower, upper);
      
      fcRes =dummy(r,function_call,environment, input);;
      
      if(fcRes<vMatValues(1,1)){ //##better
        e =(1+(*gamma)) *m -point* (*gamma);
        chkBound(e,lower,upper );
        
        double fcRes2 = dummy(e,function_call,environment, input);;
        
        if(fcRes2<fcRes){ //e better
          point = e;
          vMatValues(1,Np1) = fcRes2 ;
        }else{
          point = r;
          vMatValues(1,Np1) = fcRes ;
        }

        continue;
      }

      if(fcRes<vMatValues(1,*nParams)){
        point = r;
        vMatValues(1,Np1) = fcRes ;
        continue;
      }

      if(vMatValues(1,Np1)<fcRes){ 
        e = (*beta) *m+(1.- *beta)*point;
        //       std::cout<<"n+1 is better"<<std::endl;
      }else{
        e = (*beta) *m+(1.- *beta)*r;
        //       std::cout<<"fc is better"<<std::endl;
      }
      
      //   e.Print("cp before ckbound");
      chkBound(e,lower,upper );
      fcRes = dummy(e,function_call,environment, input);
      if(fcRes< vMatValues(1,Np1)){
        point = e;
        vMatValues(1,Np1)= fcRes ;
        continue;
      }
      
      idx = (unsigned int) vMatValues(2,1);
      RMat best(vMatbp.getColPtr(idx),*nParams ,1); 
      
      for(unsigned int k=2;k<= vMatbp.NumCols();++k){ 
        idx = (unsigned int) vMatValues(2,k);
        point.setValuesPtr(vMatbp.getColPtr(idx));
        
        e = best+ (*sigma) *( point -best);
        chkBound(e,lower,upper );
        point = e;
        fcRes = dummy(e,function_call,environment, input);
        vMatValues(1,k)= fcRes;
      }

    }
    *resPtr =vMatValues(2,1) ;
    
    UNPROTECT(3);
    
}
  
#if COMPILE_PARALLEL == 2
  void startThreads(){
    CstartThreads();
  }
  void stopThreads(){
    CstopThreads();
  }
  
#endif
}

