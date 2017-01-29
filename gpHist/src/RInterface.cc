
#include "RMat.h"
#include "R.h"     // R functions
#include "Rmath.h" // Rmath

#include <R_ext/Utils.h>
#include <Rinternals.h>


#include <iostream> // For cout, cerr, endl, and flush
#include <algorithm>
#include <stdlib.h>
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
  
}
*/

void fastMultiply( RMat& X,RMat& Y, RMat& result, RMat& orders, double sigma){
  
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
}



bool conjgradFP( RMat& A, RMat& b,RMat& x ,RMat& orders,double sigma, unsigned int max_iterations=1000 ) {

  RMat res( x.NumRows(),1 );

  fastMultiply(A,x,res,orders,sigma); //A%*%x;

  RMat r= b-res;// 
  RMat p=r;
  
  double rsold=r.SquareSum();//t(r)%*%r 
    
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
    
    
    double rsnew= r.SquareSum();//t(r)%*%r;
    if(sqrt(rsnew) <1e-12){
//#  cat(c('break:', i, ' err:',sum(sqrt(rsnew)),'\n') )
      return true;
    }
    p*=(rsnew/rsold);
    p+=r;
    //p.Print();
    rsold=rsnew;
  }
    
  return false;
}
// power method
double getEigen( RMat& A,RMat &b,RMat& orders,double sigma){
  double last = 0;
  RMat res( b.NumRows(),1 );
  
  for(unsigned int i=0;i<10000;++i){ // or give it a max iterations parameter
    fastMultiply(A,b,res,orders,sigma); //A%*%b;
    double norm = sqrt(res.SquareSum()) ;
  //  res.Print();
   // std::cout<< norm<<std::endl;
    b = res;
    b /= norm;
  //  b.Print();
    if(fabs(last-norm)<1e-10 ){
        break;
    }
    last= norm;
  }
  return last;
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

void CppHistVariance(double* result,double numRows,double numCols,double numRows2,double numCols2,double* mat1,double* mat2,double *mat3,double lambda,double sigma,double* orders){
  RMat  vMatX(mat1,numRows,numCols);
  RMat  vMatPred(mat2,numRows2,numCols2);
  RMat  vMatResult(result,numRows2,1);
  RMat  vMatAlpha(mat3,numRows,1);
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
  vMatResult *= (1/lambda);

  
 // vMatPred.RowSums().Print("rowssums") ;
  
  vMatResult = vMatPred.RowSums() -vMatResult;
  vMatResult += sigma;

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
                   double *mat1, double *mat2,double sigma , double *orders, double* logmarginal,double* lambda)
{
    RMat  vMatX(mat1,numRows,numCols);
    RMat  vMatY(mat2,numRows2,numCols2);
    RMat  vMatAlphas(result,numRows,1);
    RMat  vMatOrders(orders,numRows,numCols);
    
   
   if(!conjgradFP(vMatX, vMatY ,vMatAlphas,vMatOrders, sigma ) ){
     Rprintf("conjgradFP not converged\n");
   }
   //vMatAlphas.Print();
   
   double mu1 = numRows*sigma + vMatX.Sum();
   // nedd to get larges eigenvalue
   RMat eigenVec(numRows,1);
   eigenVec = vMatY;//1.0/(*numRows);
   
   double beta = getEigen(vMatX,eigenVec,vMatOrders, sigma);
   *lambda = beta;
//   std::cout<<beta<<std::endl;
 //  std::cout<<" "<<std::endl;
 //  eigenVec.Print();
   
   double mu2 =  beta *beta;
     
   double t_bar = (beta*mu1 - mu2) / (beta*numRows - mu1);
   double t_bar2 = t_bar*t_bar;
   double beta2 = beta*beta; /// CHECK THIS 
     
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
    
   double logdetM =  (t1 * t2).ScalarProd(t3);
   
  // std::cout<<"logdet:"<<logdetM<<std::endl;
   
   *logmarginal= 0.5*vMatAlphas.ScalarProd(vMatY)+ logdetM/2 + (numRows/2.0)*log(2.0*PI);
  // std::cout<<"   vMatY.ScalarProd( vMatAlphas ) :"<<   vMatY.ScalarProd( vMatAlphas ) <<std::endl;
   //std::cout<<"    numRows/2*log(2*PI):"<<   (numRows/2.0)*log(2.0*PI)<<std::endl;
   
  // std::cout<<"approxmarginal:"<<*logmarginal<<std::endl;


}

void Cmtest(double *result,
            int *numRows, int *numCols,int *numRows2 ,int *numCols2,
            double *mat1, double *mat2,double *sigma , double *orders)
{
  RMat  vMatX(mat1,*numRows,*numCols);
  RMat  vMatY(mat2,*numRows2,*numCols2);
  RMat  vMatAlphas(result,*numRows,1);
  RMat  vMatOrders(orders,*numRows,*numCols);
  
  fastMultiply(vMatX, vMatY ,vMatAlphas,vMatOrders, *sigma );
  
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

extern "C" {
    void CgpHist(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *sigma,double *orders,double* logmarginal,double *lambda){
      CppHist(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,*sigma,orders,logmarginal,lambda);
    }

  void CgpHistPredict(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *mat3,double *orders){
    CppHistPredict(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,mat3,orders);
  }
  
  void CgpHistVariance(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *mat3,double*lambda,double* sigma, double *orders){
    CppHistVariance(result,*numRows,*numCols,*numRows2,*numCols2,mat1,mat2,mat3,*lambda,*sigma,orders);
  }

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
          
          m.Print("\n\n\nm");
          
          idx =(int) vMatValues(2,Np1);
       //   m.Print("mean");
         // std::cout<<"calculated mean"<<std::endl;  
            
          RMat point(vMatbp.getColPtr(idx),*nParams ,1); // get vector of params
       //   point.Print("point...");
          r =(1.+(*alpha)) *m - point*(*alpha) ;//bp[N+1,1:N]
          r.Print("r");
          chkBound(r, lower, upper);
          
      //    std::cout<<"after chk bound r"<<std::endl;  
      //    r.Print("R\n");
          //np = fc(r) //querry point
          fcRes =dummy(r,function_call,environment, input);;
            
           if(fcRes<vMatValues(1,1)){ //##better
              e =(1+(*gamma)) *m -point* (*gamma);
             std::cout<<"after calc e"<<std::endl;  
             e.Print("e\n");
              chkBound(e,lower,upper );
               
               double fcRes2 = dummy(e,function_call,environment, input);;
               
               if(fcRes2<fcRes){ //e better
                 point = e;
                 vMatValues(1,Np1) = fcRes2 ;
               }else{
                 point = r;
                 vMatValues(1,Np1) = fcRes ;
               }
               vMatValues.Print("vMatValues\n");
               vMatbp.Print("vMatbp\n");
               
               continue;
           }
           
           
           if(fcRes<vMatValues(1,*nParams)){
            // idx =(int) vMatValues(2,*nParams);// set second worst solution
           //  point.setValuesPtr(vMatbp.getColPtr(idx));
             
             point = r;
             vMatValues(1,Np1) = fcRes ;
             
             vMatValues.Print("set r vMatValues\n");
             vMatbp.Print("vMatbp\n");
             
             continue;
           }
           std::cout<<"fc res"<< fcRes<<std::endl;
           std::cout<<"vMatValues(1,Np1) res"<< vMatValues(1,Np1)<<std::endl;
           
           
           if(vMatValues(1,Np1)<fcRes){ 
             e = (*beta) *m+(1.- *beta)*point;
             std::cout<<"n+1 is better"<<std::endl;
           }else{
             e = (*beta) *m+(1.- *beta)*r;
             std::cout<<"fc is better"<<std::endl;
             r.Print("my r");
           }
           
           e.Print("cp before ckbound");
           
           chkBound(e,lower,upper );
           fcRes = dummy(e,function_call,environment, input);
           e.Print("get cp query");
           if(fcRes< vMatValues(1,Np1)){
               point = e;
               vMatValues(1,Np1)= fcRes ;
               vMatValues.Print("vMatValues accteo cp");
               vMatbp.Print("vMatbp");
               
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
           
           vMatValues.Print("vMatValues last");
           vMatbp.Print("vMatbp");
           

        }
        *resPtr =vMatValues(2,1) ;
        
    UNPROTECT(3);
   
    //Cdownhillsimplex SEXP
   // std::cout<<"OK"<<std::endl;
  }
  

     
       
  void CgpfMult(double *result,double *mat1,int *numRows,int *numCols, double *mat2,int *numRows2,int *numCols2,double *sigma,double *orders)
  {
    Cmtest(result,numRows,numCols,numRows2,numCols2,mat1,mat2,sigma,orders);  
  }
  
}
